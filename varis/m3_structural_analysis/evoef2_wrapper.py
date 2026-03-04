"""EvoEF2 Wrapper — ΔΔG stability prediction (MIT licensed, primary DDG tool).

Runs EvoEF2 to compute the change in folding free energy (ΔΔG) caused by
a missense mutation. Positive ΔΔG = destabilizing.

Workflow (all in a temp directory):
1. Copy PDB into tempdir
2. RepairStructure → repaired PDB
3. Write mutation file in EvoEF2 format
4. BuildMutant → mutant PDB
5. ComputeStability on wild-type and mutant
6. DDG = mutant_energy - wt_energy

Prerequisites:
    EvoEF2 requires its ``library/`` directory (parameter files) to be located
    next to the binary. When installing manually::

        git clone https://github.com/tommyhuangthu/EvoEF2.git /tmp/EvoEF2
        cd /tmp/EvoEF2 && g++ -O3 -ffast-math -o EvoEF2 src/*.cpp
        cp EvoEF2 /usr/local/bin/
        cp -r library /usr/local/bin/library

    Or install to a custom location and set ``EVOEF2_BINARY`` env var.

Populates: ddg_evoef2.
"""

import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from varis.config import EVOEF2_BINARY
from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

EVOEF2_TIMEOUT = 120  # seconds per subprocess call


def _find_binary() -> Optional[str]:
    """Locate the EvoEF2 binary.

    Checks EVOEF2_BINARY config first, then falls back to PATH lookup.

    Returns:
        Path to binary, or None if not found.
    """
    if EVOEF2_BINARY:
        if shutil.which(EVOEF2_BINARY):
            return EVOEF2_BINARY
        logger.warning("EVOEF2_BINARY=%s not found, falling back to PATH", EVOEF2_BINARY)
    for name in ("EvoEF2", "evoef2"):
        path = shutil.which(name)
        if path:
            return path
    return None


def _parse_total_energy(stdout: str) -> Optional[float]:
    """Parse total energy from EvoEF2 ComputeStability output.

    Looks for a line matching 'Total = <number>'.

    Args:
        stdout: Standard output from EvoEF2 ComputeStability.

    Returns:
        The total energy as a float, or None if not found.
    """
    match = re.search(r"Total\s*=\s*([-+]?\d*\.?\d+)", stdout)
    if match:
        return float(match.group(1))
    return None


def _build_mutation_string(ref_aa: str, chain: str, position: int, alt_aa: str) -> str:
    """Build EvoEF2 mutation string.

    Format: {ref_aa}{chain}{position}{alt_aa};

    Args:
        ref_aa: Single-letter reference amino acid.
        chain: PDB chain identifier.
        position: Residue position number.
        alt_aa: Single-letter alternate amino acid.

    Returns:
        Mutation string, e.g. 'RA1699W;'
    """
    return f"{ref_aa}{chain}{position}{alt_aa};"


def run_evoef2(variant_record: VariantRecord) -> VariantRecord:
    """Run EvoEF2 ΔΔG prediction.

    Pre-checks: binary exists, pdb_path set, ref/alt amino acids available.
    On failure, sets ddg_evoef2=None with appropriate reason code.

    Args:
        variant_record: The shared VariantRecord with pdb_path set.

    Returns:
        The variant_record with ddg_evoef2 populated (or None with reason).
    """
    # Pre-check: binary
    binary = _find_binary()
    if not binary:
        logger.info("EvoEF2 binary not found")
        variant_record.set_with_reason("ddg_evoef2", None, NullReason.TOOL_MISSING)
        return variant_record

    # Pre-check: structure
    pdb_path = variant_record.pdb_fixed_path or variant_record.pdb_path
    if not pdb_path:
        logger.info("EvoEF2 skipped: no structure available")
        variant_record.set_with_reason(
            "ddg_evoef2", None, NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        return variant_record

    # Pre-check: mutation details
    ref_aa = variant_record.ref_aa_single
    alt_aa = variant_record.alt_aa_single
    position = variant_record.residue_position
    if not all([ref_aa, alt_aa, position]):
        logger.info("EvoEF2 skipped: missing ref_aa_single/alt_aa_single/residue_position")
        variant_record.set_with_reason(
            "ddg_evoef2", None, NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        return variant_record

    chain = "A"  # Default chain for AlphaFold structures

    try:
        ddg = _compute_ddg(binary, pdb_path, ref_aa, chain, position, alt_aa)
        variant_record.set_with_reason("ddg_evoef2", round(ddg, 4))
        logger.info("EvoEF2 DDG = %.4f kcal/mol", ddg)
    except subprocess.TimeoutExpired:
        logger.warning("EvoEF2 timed out")
        variant_record.set_with_reason("ddg_evoef2", None, NullReason.TIMED_OUT)
    except Exception as e:
        logger.warning("EvoEF2 failed: %s", e)
        variant_record.set_with_reason("ddg_evoef2", None, NullReason.TOOL_CRASHED)

    return variant_record


def _compute_ddg(
    binary: str,
    pdb_path: str,
    ref_aa: str,
    chain: str,
    position: int,
    alt_aa: str,
) -> float:
    """Execute the full EvoEF2 DDG workflow in a temp directory.

    Args:
        binary: Path to EvoEF2 binary.
        pdb_path: Path to the input PDB file.
        ref_aa: Single-letter reference amino acid.
        chain: PDB chain ID.
        position: Residue position.
        alt_aa: Single-letter alternate amino acid.

    Returns:
        DDG value (mutant_energy - wt_energy).

    Raises:
        RuntimeError: If any EvoEF2 step fails.
        subprocess.TimeoutExpired: If any step exceeds timeout.
    """
    with tempfile.TemporaryDirectory(prefix="evoef2_") as tmpdir:
        tmpdir = Path(tmpdir)

        # Copy PDB into tempdir
        pdb_src = Path(pdb_path)
        pdb_name = pdb_src.name
        pdb_dest = tmpdir / pdb_name
        shutil.copy2(pdb_src, pdb_dest)

        stem = pdb_src.stem  # filename without extension

        # Step 1: RepairStructure
        result = subprocess.run(
            [binary, "--command=RepairStructure", f"--pdb={pdb_name}"],
            cwd=str(tmpdir),
            capture_output=True,
            text=True,
            timeout=EVOEF2_TIMEOUT,
        )
        if result.returncode != 0:
            raise RuntimeError(f"RepairStructure failed: {result.stderr[:200]}")

        repaired_pdb = f"{stem}_Repair.pdb"
        if not (tmpdir / repaired_pdb).exists():
            raise RuntimeError(f"Repaired PDB not found: {repaired_pdb}")

        # Step 2: Write mutation file
        mutation_str = _build_mutation_string(ref_aa, chain, position, alt_aa)
        mutation_file = tmpdir / "individual_list.txt"
        mutation_file.write_text(mutation_str)

        # Step 3: BuildMutant
        result = subprocess.run(
            [
                binary,
                "--command=BuildMutant",
                f"--pdb={repaired_pdb}",
                "--mutant_file=individual_list.txt",
            ],
            cwd=str(tmpdir),
            capture_output=True,
            text=True,
            timeout=EVOEF2_TIMEOUT,
        )
        if result.returncode != 0:
            raise RuntimeError(f"BuildMutant failed: {result.stderr[:200]}")

        mutant_pdb = f"{stem}_Repair_Model_0001.pdb"
        if not (tmpdir / mutant_pdb).exists():
            raise RuntimeError(f"Mutant PDB not found: {mutant_pdb}")

        # Step 4: ComputeStability on wild-type (repaired)
        result_wt = subprocess.run(
            [binary, "--command=ComputeStability", f"--pdb={repaired_pdb}"],
            cwd=str(tmpdir),
            capture_output=True,
            text=True,
            timeout=EVOEF2_TIMEOUT,
        )
        wt_energy = _parse_total_energy(result_wt.stdout)
        if wt_energy is None:
            raise RuntimeError(
                f"Could not parse WT energy from: {result_wt.stdout[:200]}"
            )

        # Step 5: ComputeStability on mutant
        result_mut = subprocess.run(
            [binary, "--command=ComputeStability", f"--pdb={mutant_pdb}"],
            cwd=str(tmpdir),
            capture_output=True,
            text=True,
            timeout=EVOEF2_TIMEOUT,
        )
        mut_energy = _parse_total_energy(result_mut.stdout)
        if mut_energy is None:
            raise RuntimeError(
                f"Could not parse mutant energy from: {result_mut.stdout[:200]}"
            )

        return mut_energy - wt_energy
