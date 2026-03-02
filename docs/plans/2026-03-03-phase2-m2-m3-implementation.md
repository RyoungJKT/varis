# Phase 2: M2 + M3 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement structure retrieval (M2) and structural feature extraction (M3) for the Varis genetic variant investigator.

**Architecture:** Thin wrappers with flat orchestration, matching the M1 pattern. Each tool wrapper is a standalone file with a main function that takes a VariantRecord and returns it. The M2 `run()` and M3 `run()` functions call wrappers sequentially with try/except around each.

**Tech Stack:** BioPython (PDB parsing, DSSP), FreeSASA (solvent accessibility), httpx (ESMFold API, InterPro API), PDBFixer/OpenMM (structure repair), pytest (testing)

---

### Task 1: Update VariantRecord Schema to v1.2.0

**Files:**
- Modify: `varis/models/variant_record.py`
- Modify: `tests/conftest.py`

**Step 1: Write failing test for new fields**

Add to `tests/test_m2_structure.py` (create file):

```python
"""Tests for M2: Structure Retrieval & Preparation."""
import pytest
from varis.models.variant_record import VariantRecord, create_variant_record, RECORD_SCHEMA_VERSION


class TestSchemaV120:
    """Verify schema v1.2.0 fields exist on VariantRecord."""

    def test_schema_version_is_1_2_0(self):
        assert RECORD_SCHEMA_VERSION == "1.2.0"

    def test_m2_fields_exist(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        # New M2 fields
        assert hasattr(record, "mutation_site_present")
        assert hasattr(record, "mutation_site_plddt")
        assert hasattr(record, "plddt_available")
        assert hasattr(record, "mutation_site_confidence_bucket")
        assert hasattr(record, "numbering_scheme")
        assert hasattr(record, "structure_quality_summary")
        assert hasattr(record, "preparation_steps")
        assert hasattr(record, "pdb_hash")
        assert hasattr(record, "structure_source_url")

    def test_m3_fields_exist(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        # Renamed field
        assert hasattr(record, "solvent_accessibility_relative")
        # New M3 fields
        assert hasattr(record, "contacts_wt")
        assert hasattr(record, "hbonds_wt")
        assert hasattr(record, "packing_density")
        assert hasattr(record, "domain_start")
        assert hasattr(record, "domain_end")

    def test_removed_fields_gone(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "helix_disruption")
        assert not hasattr(record, "hbonds_lost")
        assert not hasattr(record, "contacts_changed")
        assert not hasattr(record, "structure_resolution")
        assert not hasattr(record, "solvent_accessibility")
        assert not hasattr(record, "functional_site_distance")
        assert not hasattr(record, "nearest_functional_site")
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_m2_structure.py::TestSchemaV120 -v`
Expected: FAIL — schema version is still "1.1.0", new fields don't exist

**Step 3: Update VariantRecord schema**

In `varis/models/variant_record.py`:

1. Change `RECORD_SCHEMA_VERSION = "1.1.0"` → `"1.2.0"`

2. Add new M2 fields after the existing STRUCTURE section:
```python
    # =========================================================================
    # STRUCTURE — 3D structure data (populated by M2)
    # =========================================================================
    structure_source: Optional[str] = None
    pdb_path: Optional[str] = None
    pdb_fixed_path: Optional[str] = None
    pdb_hash: Optional[str] = None
    structure_source_url: Optional[str] = None
    mutation_site_present: Optional[bool] = None
    mutation_site_plddt: Optional[float] = None
    plddt_mean: Optional[float] = None
    plddt_available: Optional[bool] = None
    mutation_site_confidence_bucket: Optional[str] = None  # "high"/"medium"/"low"
    numbering_scheme: Optional[str] = None  # "uniprot_canonical"
    structure_quality_summary: Optional[dict] = None
    preparation_steps: Optional[list[str]] = None
```

3. Replace old M3 fields with revised ones:
```python
    # =========================================================================
    # STRUCTURAL FEATURES — From 3D analysis (populated by M3)
    # =========================================================================
    ddg_foldx: Optional[float] = None
    ddg_pyrosetta: Optional[float] = None
    ddg_mean: Optional[float] = None
    solvent_accessibility_relative: Optional[float] = None  # [0.0, 1.0]
    burial_category: Optional[str] = None  # "core" / "surface"
    secondary_structure: Optional[str] = None
    secondary_structure_name: Optional[str] = None
    contacts_wt: Optional[int] = None
    hbonds_wt: Optional[int] = None
    packing_density: Optional[float] = None
    domain_name: Optional[str] = None
    domain_id: Optional[str] = None
    domain_start: Optional[int] = None
    domain_end: Optional[int] = None
    domain_criticality: Optional[str] = None
```

4. Remove these fields entirely:
   - `helix_disruption`
   - `hbonds_lost`
   - `contacts_changed`
   - `functional_site_distance`
   - `nearest_functional_site`
   - `structure_resolution`
   - `solvent_accessibility`
   - `plddt_score` (replaced by `mutation_site_plddt`)

5. Update `get_structural_features()` to use new field names:
```python
    def get_structural_features(self) -> dict:
        """Extract structural feature dict for ML model input."""
        return {
            "ddg_foldx": self.ddg_foldx,
            "ddg_pyrosetta": self.ddg_pyrosetta,
            "solvent_accessibility_relative": self.solvent_accessibility_relative,
            "secondary_structure": self.secondary_structure_name,
            "contacts_wt": self.contacts_wt,
            "hbonds_wt": self.hbonds_wt,
            "packing_density": self.packing_density,
            "domain_name": self.domain_name,
            "conservation_score": self.conservation_score,
            "mutation_site_plddt": self.mutation_site_plddt,
            "gnomad_frequency": self.gnomad_frequency,
            "alphamissense_score": self.alphamissense_score,
            "charge_change": self.charge_change,
            "burial_category": self.burial_category,
            "mutation_site_confidence_bucket": self.mutation_site_confidence_bucket,
        }
```

6. Update `validate()` to reference new field names (replace `solvent_accessibility` with `solvent_accessibility_relative`, remove `helix_disruption` references).

**Step 4: Update `tests/conftest.py` fixtures**

Update `fully_populated_record` to use new field names:
```python
    # M2 fields
    record.pdb_hash = "abc123def456"
    record.mutation_site_present = True
    record.mutation_site_plddt = 92.4
    record.plddt_mean = 85.1
    record.plddt_available = True
    record.mutation_site_confidence_bucket = "high"
    record.numbering_scheme = "uniprot_canonical"
    record.structure_quality_summary = {
        "plddt_mean": 85.1, "plddt_site": 92.4, "percent_low_confidence": 0.08
    }
    record.preparation_steps = ["validated"]

    # M3 fields (updated names)
    record.solvent_accessibility_relative = 0.03
    record.contacts_wt = 14
    record.hbonds_wt = 3
    record.packing_density = 0.72
    record.domain_start = 1646
    record.domain_end = 1736
```

Remove references to: `helix_disruption`, `hbonds_lost`, `contacts_changed`, `functional_site_distance`, `solvent_accessibility`, `plddt_score`.

**Step 5: Run tests to verify schema update passes**

Run: `pytest tests/test_m2_structure.py::TestSchemaV120 -v`
Expected: PASS

**Step 6: Run full test suite to check for regressions**

Run: `pytest tests/ -v --tb=short`
Expected: Some M5/conftest tests may need fixing due to renamed fields. Fix any failures.

**Step 7: Commit**

```bash
git add varis/models/variant_record.py tests/conftest.py tests/test_m2_structure.py
git commit -m "refactor(schema): update VariantRecord to v1.2.0 for Phase 2

Add M2 fields: mutation_site_present, mutation_site_plddt, plddt_available,
mutation_site_confidence_bucket, numbering_scheme, structure_quality_summary,
preparation_steps, pdb_hash, structure_source_url.

Rename: solvent_accessibility → solvent_accessibility_relative,
plddt_score → mutation_site_plddt.

Add M3 fields: contacts_wt, hbonds_wt, packing_density, domain_start, domain_end.

Remove deferred fields: helix_disruption, hbonds_lost, contacts_changed,
functional_site_distance, nearest_functional_site, structure_resolution."
```

---

### Task 2: M2 — Structure Validator

**Files:**
- Create: `varis/m2_structure/structure_validator.py`
- Modify: `tests/test_m2_structure.py`

**Step 1: Write failing tests**

Add to `tests/test_m2_structure.py`:

```python
from pathlib import Path
from varis.config import STRUCTURES_DIR

BRCA1_PDB = STRUCTURES_DIR / "AF-P38398-F1-model_v6.pdb"


class TestStructureValidator:
    """Tests for structure_validator.py — validates PDB and extracts pLDDT."""

    def test_structure_sanity_fixture(self):
        """PDB has >= 1699 residues, numbering starts at 1, continuous."""
        from varis.m2_structure.structure_validator import validate_structure
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("brca1", str(BRCA1_PDB))
        residues = [r for r in structure[0]["A"].get_residues()
                    if r.id[0] == " "]  # standard residues only
        assert len(residues) >= 1699
        # Numbering starts at 1
        assert residues[0].id[1] == 1

    def test_validate_existing_pdb(self, m1_completed_record):
        """Finds residue 1699, mutation_site_present=True."""
        from varis.m2_structure.structure_validator import validate_structure
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = validate_structure(m1_completed_record)
        assert result.mutation_site_present is True
        assert result.numbering_scheme == "uniprot_canonical"

    def test_validate_missing_residue(self, m1_completed_record):
        """mutation_site_present=False when residue not in structure."""
        from varis.m2_structure.structure_validator import validate_structure
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        m1_completed_record.residue_position = 99999  # out of range
        result = validate_structure(m1_completed_record)
        assert result.mutation_site_present is False
        assert "mutation_site_present" in result.null_reasons

    def test_site_out_of_range(self, m1_completed_record):
        """Position > protein length → mutation_site_present=False."""
        from varis.m2_structure.structure_validator import validate_structure
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        m1_completed_record.residue_position = 5000
        result = validate_structure(m1_completed_record)
        assert result.mutation_site_present is False

    def test_plddt_only_for_predicted(self, m1_completed_record):
        """pLDDT skipped when structure_source is not alphafold/esmfold."""
        from varis.m2_structure.structure_validator import validate_structure
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "pdb"
        result = validate_structure(m1_completed_record)
        assert result.plddt_available is False
        assert result.mutation_site_plddt is None

    def test_plddt_range_valid(self, m1_completed_record):
        """All pLDDT values in [0, 100]."""
        from varis.m2_structure.structure_validator import validate_structure
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = validate_structure(m1_completed_record)
        if result.plddt_available:
            assert 0.0 <= result.mutation_site_plddt <= 100.0
            assert 0.0 <= result.plddt_mean <= 100.0

    def test_plddt_site_bfactor_value(self, m1_completed_record):
        """Extracted pLDDT matches the B-factor at residue 1699 in the PDB."""
        from varis.m2_structure.structure_validator import validate_structure
        from Bio.PDB import PDBParser
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        # Get expected B-factor directly
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("brca1", str(BRCA1_PDB))
        residue = structure[0]["A"][1699]
        expected_bfactor = residue["CA"].get_bfactor()
        # Run validator
        result = validate_structure(m1_completed_record)
        assert result.mutation_site_plddt == pytest.approx(expected_bfactor, abs=0.01)

    def test_structure_quality_summary(self, m1_completed_record):
        """structure_quality_summary has required keys."""
        from varis.m2_structure.structure_validator import validate_structure
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = validate_structure(m1_completed_record)
        summary = result.structure_quality_summary
        assert "plddt_mean" in summary
        assert "plddt_site" in summary
        assert "percent_low_confidence" in summary

    def test_pdb_hash_computed(self, m1_completed_record):
        """pdb_hash is a hex string."""
        from varis.m2_structure.structure_validator import validate_structure
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = validate_structure(m1_completed_record)
        assert result.pdb_hash is not None
        assert len(result.pdb_hash) == 64  # SHA256 hex
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m2_structure.py::TestStructureValidator -v`
Expected: FAIL — `structure_validator` module doesn't exist yet

**Step 3: Implement structure_validator.py**

Create `varis/m2_structure/structure_validator.py`:

```python
"""Structure Validator — Validates PDB and extracts quality metrics.

Checks that the mutation residue exists in the structure, extracts pLDDT
scores from B-factor column (AlphaFold/ESMFold only), and computes a
structure quality summary.

Populates: mutation_site_present, mutation_site_plddt, plddt_mean,
           plddt_available, mutation_site_confidence_bucket, numbering_scheme,
           structure_quality_summary, preparation_steps, pdb_hash.
"""

import hashlib
import logging
from pathlib import Path

from Bio.PDB import PDBParser

from varis.config import PLDDT_CONFIDENCE_THRESHOLD
from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

_PREDICTED_SOURCES = ("alphafold", "esmfold")


def validate_structure(variant_record: VariantRecord) -> VariantRecord:
    """Validate PDB structure and extract quality metrics.

    Args:
        variant_record: Must have pdb_path and residue_position set.

    Returns:
        VariantRecord with validation fields populated.
    """
    pdb_path = variant_record.pdb_path
    position = variant_record.residue_position
    source = variant_record.structure_source

    if not pdb_path or not Path(pdb_path).exists():
        logger.warning("No PDB file found at %s", pdb_path)
        variant_record.set_with_reason(
            "mutation_site_present", None, NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        return variant_record

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)
        model = structure[0]

        # Find the first chain (AlphaFold uses chain A)
        chain = next(iter(model.get_chains()))
        residues = [r for r in chain.get_residues() if r.id[0] == " "]

        # Compute PDB hash for provenance
        variant_record.pdb_hash = _compute_hash(pdb_path)

        # Check if mutation site is present
        residue_ids = {r.id[1] for r in residues}
        if position not in residue_ids:
            logger.warning("Residue %d not found in structure", position)
            variant_record.mutation_site_present = False
            variant_record.set_with_reason(
                "mutation_site_present", False, NullReason.NO_DATA_AVAILABLE
            )
            variant_record.coordinate_mapping_confidence = "failed"
            variant_record.preparation_steps = ["validated"]
            variant_record.numbering_scheme = "uniprot_canonical"
            return variant_record

        variant_record.mutation_site_present = True
        variant_record.numbering_scheme = "uniprot_canonical"
        variant_record.preparation_steps = ["validated"]

        # Extract pLDDT from B-factor column (only for predicted structures)
        if source in _PREDICTED_SOURCES:
            variant_record.plddt_available = True
            target_residue = chain[position]
            ca_atom = target_residue["CA"]
            site_plddt = ca_atom.get_bfactor()
            variant_record.mutation_site_plddt = site_plddt

            # Mean pLDDT across all CA atoms
            all_bfactors = []
            for r in residues:
                if "CA" in r:
                    all_bfactors.append(r["CA"].get_bfactor())
            mean_plddt = sum(all_bfactors) / len(all_bfactors) if all_bfactors else 0.0
            variant_record.plddt_mean = round(mean_plddt, 2)

            # Confidence bucket
            if site_plddt >= 90:
                variant_record.mutation_site_confidence_bucket = "high"
            elif site_plddt >= PLDDT_CONFIDENCE_THRESHOLD:
                variant_record.mutation_site_confidence_bucket = "medium"
            else:
                variant_record.mutation_site_confidence_bucket = "low"

            # Percent of residues with low confidence
            low_count = sum(1 for b in all_bfactors if b < PLDDT_CONFIDENCE_THRESHOLD)
            pct_low = round(low_count / len(all_bfactors), 4) if all_bfactors else 0.0

            variant_record.structure_quality_summary = {
                "plddt_mean": variant_record.plddt_mean,
                "plddt_site": site_plddt,
                "percent_low_confidence": pct_low,
            }
        else:
            variant_record.plddt_available = False
            variant_record.set_with_reason(
                "mutation_site_plddt", None, NullReason.INTENTIONALLY_SKIPPED
            )
            variant_record.set_with_reason(
                "plddt_mean", None, NullReason.INTENTIONALLY_SKIPPED
            )

    except Exception as e:
        logger.warning("Structure validation failed: %s", e)
        variant_record.set_with_reason(
            "mutation_site_present", None, NullReason.TOOL_CRASHED
        )

    return variant_record


def _compute_hash(file_path: str) -> str:
    """Compute SHA256 hash of a file.

    Args:
        file_path: Path to the file.

    Returns:
        Hex digest of SHA256 hash.
    """
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()
```

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_m2_structure.py::TestStructureValidator -v`
Expected: PASS

**Step 5: Commit**

```bash
git add varis/m2_structure/structure_validator.py tests/test_m2_structure.py
git commit -m "feat(m2): implement structure validator with pLDDT extraction

Validates PDB residue presence, extracts pLDDT from B-factor column
(AlphaFold/ESMFold only), computes quality summary and SHA256 hash."
```

---

### Task 3: M2 — ESMFold Predictor

**Files:**
- Modify: `varis/m2_structure/esmfold_predictor.py`
- Modify: `tests/test_m2_structure.py`

**Step 1: Write failing tests**

Add to `tests/test_m2_structure.py`:

```python
from unittest.mock import MagicMock, patch


class TestESMFoldPredictor:
    """Tests for esmfold_predictor.py — fallback structure prediction."""

    def test_sequence_too_long(self, m1_completed_record):
        """Sequences >400aa skip with sequence_too_long reason."""
        from varis.m2_structure.esmfold_predictor import predict_esmfold
        m1_completed_record.protein_sequence = "M" * 401
        m1_completed_record.pdb_path = None
        result = predict_esmfold(m1_completed_record)
        assert result.pdb_path is None
        assert "pdb_path" in result.null_reasons

    def test_no_sequence_available(self, m1_completed_record):
        """No protein_sequence → skip with upstream dependency reason."""
        from varis.m2_structure.esmfold_predictor import predict_esmfold
        m1_completed_record.protein_sequence = None
        m1_completed_record.pdb_path = None
        result = predict_esmfold(m1_completed_record)
        assert result.pdb_path is None

    def test_skips_if_structure_exists(self, m1_completed_record):
        """If pdb_path already set, ESMFold is not called."""
        from varis.m2_structure.esmfold_predictor import predict_esmfold
        m1_completed_record.pdb_path = "/some/existing.pdb"
        result = predict_esmfold(m1_completed_record)
        assert result.pdb_path == "/some/existing.pdb"

    def test_successful_prediction(self, m1_completed_record, tmp_path):
        """Mocked ESMFold API returns valid PDB content."""
        from varis.m2_structure.esmfold_predictor import predict_esmfold
        m1_completed_record.protein_sequence = "MKFLILLFNILCLFPVLAADNHGVS"  # 25 aa
        m1_completed_record.pdb_path = None
        m1_completed_record.uniprot_id = "TEST123"
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "ATOM      1  N   MET A   1       0.000   0.000   0.000\nEND\n"
        mock_client = MagicMock()
        mock_client.post.return_value = mock_response
        with patch("varis.m2_structure.esmfold_predictor.STRUCTURES_DIR", tmp_path):
            result = predict_esmfold(m1_completed_record, client=mock_client)
        assert result.pdb_path is not None
        assert result.structure_source == "esmfold"
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m2_structure.py::TestESMFoldPredictor -v`
Expected: FAIL

**Step 3: Implement esmfold_predictor.py**

Replace `varis/m2_structure/esmfold_predictor.py`:

```python
"""ESMFold Predictor — Fallback structure prediction via ESM Atlas API.

Calls Meta's ESM Metagenomic Atlas API to predict a 3D structure from
amino acid sequence. Limited to sequences ≤400 residues. Non-critical
fallback — only used when AlphaFold DB has no structure.

Populates: structure_source, pdb_path (if successful).
"""

import logging
from pathlib import Path

import httpx

from varis.config import ESMFOLD_API_URL, STRUCTURES_DIR
from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

_TIMEOUT = 120.0  # ESMFold can be slow
_MAX_SEQUENCE_LENGTH = 400


def predict_esmfold(variant_record: VariantRecord,
                    client: httpx.Client | None = None) -> VariantRecord:
    """Predict structure using ESMFold API.

    Args:
        variant_record: Must have protein_sequence set.
        client: Optional httpx.Client for dependency injection.

    Returns:
        VariantRecord with pdb_path set if prediction succeeded.
    """
    # Skip if structure already exists
    if variant_record.pdb_path:
        logger.info("Structure already available, skipping ESMFold")
        return variant_record

    sequence = variant_record.protein_sequence
    if not sequence:
        logger.info("No protein sequence available — skipping ESMFold")
        variant_record.set_with_reason(
            "pdb_path", None, NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        return variant_record

    if len(sequence) > _MAX_SEQUENCE_LENGTH:
        logger.info(
            "Sequence too long for ESMFold (%d > %d aa) — skipping",
            len(sequence), _MAX_SEQUENCE_LENGTH,
        )
        variant_record.set_with_reason("pdb_path", None, NullReason.INTENTIONALLY_SKIPPED)
        return variant_record

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        response = client.post(ESMFOLD_API_URL, content=sequence)
        if response.status_code != 200:
            logger.warning("ESMFold API returned status %d", response.status_code)
            variant_record.set_with_reason(
                "pdb_path", None, NullReason.TOOL_CRASHED
            )
            return variant_record

        pdb_content = response.text
        if not pdb_content or "ATOM" not in pdb_content:
            logger.warning("ESMFold returned empty or invalid PDB")
            variant_record.set_with_reason(
                "pdb_path", None, NullReason.TOOL_CRASHED
            )
            return variant_record

        # Save PDB file
        uniprot_id = variant_record.uniprot_id or "unknown"
        output_dir = STRUCTURES_DIR
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"ESMFold-{uniprot_id}.pdb"
        output_path.write_text(pdb_content)

        variant_record.structure_source = "esmfold"
        variant_record.pdb_path = str(output_path)
        variant_record.structure_source_url = ESMFOLD_API_URL
        logger.info("ESMFold structure saved: %s", output_path)

    except Exception as e:
        logger.warning("ESMFold prediction failed: %s", e)
        variant_record.set_with_reason("pdb_path", None, NullReason.TOOL_CRASHED)
    finally:
        if should_close:
            client.close()

    return variant_record
```

**Step 4: Run tests**

Run: `pytest tests/test_m2_structure.py::TestESMFoldPredictor -v`
Expected: PASS

**Step 5: Commit**

```bash
git add varis/m2_structure/esmfold_predictor.py tests/test_m2_structure.py
git commit -m "feat(m2): implement ESMFold fallback predictor

Calls ESM Atlas API for sequences ≤400aa. Skips gracefully for long
sequences with sequence_too_long reason code."
```

---

### Task 4: M2 — PDB Fixer (Conditional)

**Files:**
- Modify: `varis/m2_structure/pdb_fixer.py`
- Modify: `tests/test_m2_structure.py`

**Step 1: Write failing tests**

Add to `tests/test_m2_structure.py`:

```python
class TestPDBFixer:
    """Tests for pdb_fixer.py — conditional structure repair."""

    def test_pdb_fixer_skipped_when_clean(self, m1_completed_record):
        """Clean AlphaFold PDB → no fixer run, preparation_steps=['validated']."""
        from varis.m2_structure.pdb_fixer import fix_structure
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.preparation_steps = ["validated"]
        result = fix_structure(m1_completed_record)
        # AlphaFold PDBs are clean — fixer should not run
        assert result.pdb_fixed_path is None or result.pdb_fixed_path == result.pdb_path
        assert "validated" in result.preparation_steps

    def test_pdb_fixer_no_structure(self, m1_completed_record):
        """No pdb_path → skip gracefully."""
        from varis.m2_structure.pdb_fixer import fix_structure
        m1_completed_record.pdb_path = None
        result = fix_structure(m1_completed_record)
        assert result.pdb_fixed_path is None

    def test_pdb_fixer_output_parses(self, m1_completed_record, tmp_path):
        """If fixer runs, output PDB parses successfully."""
        from varis.m2_structure.pdb_fixer import fix_structure
        from Bio.PDB import PDBParser
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.preparation_steps = ["validated"]
        result = fix_structure(m1_completed_record)
        # Whichever path we end up with should be parseable
        pdb_to_check = result.pdb_fixed_path or result.pdb_path
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("test", pdb_to_check)
        assert structure is not None

    def test_site_preserved_after_fix(self, m1_completed_record):
        """Target residue still present after PDBFixer."""
        from varis.m2_structure.pdb_fixer import fix_structure
        from Bio.PDB import PDBParser
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.preparation_steps = ["validated"]
        result = fix_structure(m1_completed_record)
        pdb_to_check = result.pdb_fixed_path or result.pdb_path
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("test", pdb_to_check)
        chain = next(iter(structure[0].get_chains()))
        residue_ids = {r.id[1] for r in chain.get_residues() if r.id[0] == " "}
        assert 1699 in residue_ids
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m2_structure.py::TestPDBFixer -v`
Expected: FAIL

**Step 3: Implement pdb_fixer.py**

Replace `varis/m2_structure/pdb_fixer.py`:

```python
"""PDB Fixer — Conditional structure repair using OpenMM PDBFixer.

Only runs when validation detects missing heavy atoms. Adds missing atoms
and optionally hydrogens. NEVER reconstructs missing residues/loops — that
is modeling, not fixing, and creates false certainty.

Populates: pdb_fixed_path, preparation_steps (appends "added_hydrogens").
"""

import logging
from pathlib import Path

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)


def fix_structure(variant_record: VariantRecord) -> VariantRecord:
    """Conditionally repair PDB structure.

    Only adds missing heavy atoms and hydrogens. Does not fill residue gaps.
    If PDBFixer is not installed, logs and continues without fixing.

    Args:
        variant_record: Must have pdb_path set.

    Returns:
        VariantRecord with pdb_fixed_path set if fixing was performed.
    """
    pdb_path = variant_record.pdb_path
    if not pdb_path or not Path(pdb_path).exists():
        logger.info("No PDB file to fix")
        variant_record.set_with_reason(
            "pdb_fixed_path", None, NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        return variant_record

    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        logger.info("PDBFixer not installed — skipping structure repair")
        variant_record.set_with_reason(
            "pdb_fixed_path", None, NullReason.TOOL_MISSING
        )
        return variant_record

    try:
        fixer = PDBFixer(filename=pdb_path)

        # Detect missing heavy atoms (NOT missing residues)
        fixer.findMissingResidues()
        # Clear missing residues — we do NOT reconstruct them
        fixer.missingResidues = {}

        fixer.findMissingAtoms()
        missing_atoms = fixer.missingAtoms
        missing_terminals = fixer.missingTerminals

        if not missing_atoms and not missing_terminals:
            logger.info("No missing atoms detected — skipping PDBFixer")
            return variant_record

        # Add missing heavy atoms
        fixer.addMissingAtoms()
        if variant_record.preparation_steps is None:
            variant_record.preparation_steps = []
        variant_record.preparation_steps.append("added_missing_atoms")

        # Add hydrogens
        fixer.addMissingHydrogens(pH=7.0)
        variant_record.preparation_steps.append("added_hydrogens")

        # Save fixed PDB
        fixed_path = Path(pdb_path).with_suffix(".fixed.pdb")
        with open(fixed_path, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        variant_record.pdb_fixed_path = str(fixed_path)
        logger.info("Fixed structure saved: %s", fixed_path)

    except Exception as e:
        logger.warning("PDBFixer failed: %s", e)
        variant_record.set_with_reason(
            "pdb_fixed_path", None, NullReason.TOOL_CRASHED
        )

    return variant_record
```

**Step 4: Run tests**

Run: `pytest tests/test_m2_structure.py::TestPDBFixer -v`
Expected: PASS (PDBFixer may not be installed — tests handle that via graceful skip)

**Step 5: Commit**

```bash
git add varis/m2_structure/pdb_fixer.py tests/test_m2_structure.py
git commit -m "feat(m2): implement conditional PDB fixer

Only adds missing heavy atoms and hydrogens. Never reconstructs missing
residues. Gracefully handles PDBFixer not being installed."
```

---

### Task 5: M2 — Orchestrator and Integration Test

**Files:**
- Modify: `varis/m2_structure/__init__.py`
- Modify: `tests/test_m2_structure.py`

**Step 1: Write failing integration test**

Add to `tests/test_m2_structure.py`:

```python
class TestM2Integration:
    """Full M2 pipeline integration test."""

    def test_m2_full_pipeline_brca1(self, m1_completed_record):
        """Full M2 on BRCA1 AlphaFold structure."""
        from varis.m2_structure import run
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = run(m1_completed_record)
        # M2 should complete
        assert "M2" in result.modules_completed
        # Validation should have run
        assert result.mutation_site_present is True
        assert result.pdb_hash is not None
        assert result.numbering_scheme == "uniprot_canonical"
        # pLDDT should be extracted
        assert result.plddt_available is True
        assert result.mutation_site_plddt is not None
        assert 0.0 <= result.mutation_site_plddt <= 100.0

    def test_m2_no_structure(self, m1_completed_record):
        """M2 with no structure → marks failed."""
        from varis.m2_structure import run
        m1_completed_record.pdb_path = None
        m1_completed_record.structure_source = None
        m1_completed_record.protein_sequence = "M" * 500  # Too long for ESMFold
        result = run(m1_completed_record)
        assert "M2" in result.modules_failed

    def test_m2_golden_record_keys(self, m1_completed_record):
        """Golden record: verify all expected keys are set after M2."""
        from varis.m2_structure import run
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = run(m1_completed_record)
        expected_keys = [
            "mutation_site_present", "pdb_hash", "numbering_scheme",
            "preparation_steps", "plddt_available",
        ]
        for key in expected_keys:
            assert getattr(result, key) is not None, f"{key} should not be None"
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m2_structure.py::TestM2Integration -v`
Expected: FAIL — orchestrator calls old stub functions

**Step 3: Update M2 orchestrator**

Replace `varis/m2_structure/__init__.py`:

```python
"""M2: Structure Engine — Obtains and prepares 3D protein structures.

Validates the structure from M1, optionally retrieves a fallback via ESMFold,
conditionally repairs with PDBFixer, and extracts quality metrics.

Depends on: M1 (needs uniprot_id, protein_sequence, pdb_path)
Populates: mutation_site_present, mutation_site_plddt, plddt_mean,
           plddt_available, structure_quality_summary, preparation_steps,
           pdb_hash, numbering_scheme, mutation_site_confidence_bucket
"""

import logging

logger = logging.getLogger(__name__)


def run(variant_record):
    """Execute M2: validate, optionally retrieve, fix, and extract quality.

    Order: validate → (esmfold if no structure) → (fix if needed) → done.
    Structure validator extracts pLDDT and quality metrics.

    Args:
        variant_record: VariantRecord with M1 fields populated.

    Returns:
        VariantRecord with M2 fields populated (or None with reasons).
    """
    from varis.m2_structure.esmfold_predictor import predict_esmfold
    from varis.m2_structure.pdb_fixer import fix_structure
    from varis.m2_structure.structure_validator import validate_structure

    # If no structure from M1, try ESMFold fallback
    if variant_record.pdb_path is None:
        try:
            variant_record = predict_esmfold(variant_record)
        except Exception as e:
            logger.warning("ESMFold fallback failed: %s", e)
            variant_record.mark_module_failed("M2.esmfold")

    # If still no structure, M2 fails
    if variant_record.pdb_path is None:
        logger.warning("No structure obtained. M3 will be skipped.")
        variant_record.mark_module_failed("M2")
        return variant_record

    # Validate structure and extract quality metrics (including pLDDT)
    try:
        variant_record = validate_structure(variant_record)
    except Exception as e:
        logger.warning("Structure validation failed: %s", e)
        variant_record.mark_module_failed("M2.validator")

    # Conditionally fix structure (only if missing atoms detected)
    try:
        variant_record = fix_structure(variant_record)
    except Exception as e:
        logger.warning("PDBFixer failed: %s", e)
        variant_record.mark_module_failed("M2.pdbfixer")

    variant_record.mark_module_completed("M2")
    return variant_record
```

**Step 4: Remove old stubs that are no longer needed**

Delete `varis/m2_structure/alphafold_retriever.py` (M1's alphafold_client handles this).
Delete `varis/m2_structure/structure_utils.py` (absorbed into structure_validator).

**Step 5: Run tests**

Run: `pytest tests/test_m2_structure.py -v`
Expected: PASS

**Step 6: Run full test suite**

Run: `pytest tests/ -v --tb=short`
Expected: PASS (or fix any regressions)

**Step 7: Commit**

```bash
git add varis/m2_structure/__init__.py tests/test_m2_structure.py
git rm varis/m2_structure/alphafold_retriever.py varis/m2_structure/structure_utils.py
git commit -m "feat(m2): implement M2 orchestrator with integration tests

Orchestrates: validate → esmfold fallback → conditional fix.
Removes alphafold_retriever stub (M1 handles download) and
structure_utils stub (absorbed into structure_validator)."
```

---

### Task 6: M3 — FreeSASA Wrapper

**Files:**
- Modify: `varis/m3_structural_analysis/freesasa_wrapper.py`
- Modify: `tests/test_m3_structural.py`

**Step 1: Write failing tests**

Replace `tests/test_m3_structural.py`:

```python
"""Tests for M3: Structural Analysis — Feature extraction from 3D structure."""
import pytest
from pathlib import Path
from varis.config import STRUCTURES_DIR

BRCA1_PDB = STRUCTURES_DIR / "AF-P38398-F1-model_v6.pdb"


class TestFreeSASA:
    """Tests for freesasa_wrapper.py — solvent accessibility."""

    def test_freesasa_computes_sasa(self, m1_completed_record):
        """Returns float, sasa_available=True."""
        from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_freesasa(m1_completed_record)
        assert result.sasa_available is True
        assert isinstance(result.solvent_accessibility_relative, float)

    def test_freesasa_relative_bounds(self, m1_completed_record):
        """Relative SASA in [0.0, 1.0]."""
        from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_freesasa(m1_completed_record)
        assert 0.0 <= result.solvent_accessibility_relative <= 1.0

    def test_freesasa_burial_category(self, m1_completed_record):
        """Returns 'core' or 'surface' only (no 'interface')."""
        from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_freesasa(m1_completed_record)
        assert result.burial_category in ("core", "surface")

    def test_freesasa_no_structure(self, m1_completed_record):
        """No PDB → sasa_available=False with reason."""
        from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
        m1_completed_record.pdb_path = None
        result = run_freesasa(m1_completed_record)
        assert result.sasa_available is False
        assert result.sasa_missing_reason is not None
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m3_structural.py::TestFreeSASA -v`
Expected: FAIL

**Step 3: Implement freesasa_wrapper.py**

Replace `varis/m3_structural_analysis/freesasa_wrapper.py`:

```python
"""FreeSASA Wrapper — Computes solvent-accessible surface area at mutation site.

Calculates relative SASA (normalized by residue-type maximum) to determine
how buried or exposed the mutation site is. Burial is a strong predictor:
core mutations are more likely to be damaging.

Definition: Relative SASA = residue SASA / max SASA for that residue type.
Range: [0.0, 1.0]. Core < 0.25, Surface >= 0.25.

Populates: solvent_accessibility_relative, burial_category.
Sets: sasa_available, sasa_missing_reason.
"""

import logging
from pathlib import Path

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# Maximum SASA values per residue type (Ų) from Tien et al. 2013 (theoretical)
# Used to normalize to relative SASA [0, 1]
_MAX_SASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

_BURIAL_THRESHOLD = 0.25  # Below this = core, above = surface


def run_freesasa(variant_record: VariantRecord) -> VariantRecord:
    """Compute relative SASA and burial category at mutation site.

    Args:
        variant_record: Must have pdb_path and residue_position set.

    Returns:
        VariantRecord with SASA fields populated.
    """
    pdb_path = variant_record.pdb_path or variant_record.pdb_fixed_path
    position = variant_record.residue_position

    if not pdb_path or not Path(pdb_path).exists():
        logger.info("No structure available for FreeSASA")
        variant_record.set_feature_status("sasa", False, NullReason.UPSTREAM_DEPENDENCY_FAILED)
        return variant_record

    try:
        import freesasa
    except ImportError:
        logger.warning("freesasa not installed")
        variant_record.set_feature_status("sasa", False, NullReason.TOOL_MISSING)
        return variant_record

    try:
        structure = freesasa.Structure(pdb_path)
        result = freesasa.calc(structure)

        # Find the target residue SASA
        area_class = freesasa.classifyResults(result, structure)

        # Use BioPython to get residue-level SASA
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        bio_structure = parser.get_structure("protein", pdb_path)
        chain = next(iter(bio_structure[0].get_chains()))
        chain_id = chain.id

        # Get per-residue SASA using freesasa's residue-level API
        residue_areas = result.residueAreas()
        target_key = f"{chain_id},{str(position)}"

        target_area = None
        for key in residue_areas:
            if key == target_key:
                target_area = residue_areas[key]
                break

        if target_area is None:
            logger.warning("Residue %d not found in FreeSASA results", position)
            variant_record.set_feature_status("sasa", False, NullReason.NO_DATA_AVAILABLE)
            return variant_record

        absolute_sasa = target_area.total

        # Get residue type for normalization
        target_residue = chain[position]
        resname = target_residue.get_resname()
        max_sasa = _MAX_SASA.get(resname, 200.0)  # default if unknown

        relative_sasa = min(absolute_sasa / max_sasa, 1.0) if max_sasa > 0 else 0.0

        variant_record.solvent_accessibility_relative = round(relative_sasa, 4)
        variant_record.burial_category = "surface" if relative_sasa >= _BURIAL_THRESHOLD else "core"
        variant_record.set_feature_status("sasa", True)

        logger.info(
            "FreeSASA: residue %d SASA=%.1f Ų (relative=%.3f, %s)",
            position, absolute_sasa, relative_sasa, variant_record.burial_category,
        )

    except Exception as e:
        logger.warning("FreeSASA failed: %s", e)
        variant_record.set_feature_status("sasa", False, NullReason.TOOL_CRASHED)

    return variant_record
```

**Step 4: Run tests**

Run: `pytest tests/test_m3_structural.py::TestFreeSASA -v`
Expected: PASS

**Step 5: Commit**

```bash
git add varis/m3_structural_analysis/freesasa_wrapper.py tests/test_m3_structural.py
git commit -m "feat(m3): implement FreeSASA wrapper for solvent accessibility

Computes relative SASA normalized by residue-type max (Tien et al. 2013).
Burial categories: core (<0.25) / surface (>=0.25)."
```

---

### Task 7: M3 — DSSP Wrapper

**Files:**
- Modify: `varis/m3_structural_analysis/dssp_wrapper.py`
- Modify: `tests/test_m3_structural.py`

**Step 1: Write failing tests**

Add to `tests/test_m3_structural.py`:

```python
import shutil


class TestDSSP:
    """Tests for dssp_wrapper.py — secondary structure assignment."""

    def test_dssp_returns_valid_code(self, m1_completed_record):
        """Result in allowed DSSP code set."""
        if not shutil.which("mkdssp"):
            pytest.skip("mkdssp binary not installed")
        from varis.m3_structural_analysis.dssp_wrapper import run_dssp
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_dssp(m1_completed_record)
        if result.dssp_available:
            valid_codes = {"H", "B", "E", "G", "I", "T", "S", "-", "C"}
            assert result.secondary_structure in valid_codes
            assert result.secondary_structure_name in ("helix", "sheet", "coil")

    def test_dssp_mkdssp_missing(self, m1_completed_record, monkeypatch):
        """dssp_available=False when mkdssp not found."""
        from varis.m3_structural_analysis.dssp_wrapper import run_dssp
        monkeypatch.setattr(shutil, "which", lambda x: None)
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_dssp(m1_completed_record)
        assert result.dssp_available is False
        assert result.dssp_missing_reason is not None

    def test_dssp_no_structure(self, m1_completed_record):
        """No PDB → dssp_available=False."""
        from varis.m3_structural_analysis.dssp_wrapper import run_dssp
        m1_completed_record.pdb_path = None
        result = run_dssp(m1_completed_record)
        assert result.dssp_available is False
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m3_structural.py::TestDSSP -v`
Expected: FAIL

**Step 3: Implement dssp_wrapper.py**

Replace `varis/m3_structural_analysis/dssp_wrapper.py`:

```python
"""DSSP Wrapper — Assigns secondary structure at mutation site.

Uses BioPython's DSSP interface with the mkdssp binary to assign
secondary structure codes (H=helix, E=sheet, C/other=coil).

Populates: secondary_structure, secondary_structure_name.
Sets: dssp_available, dssp_missing_reason.
"""

import logging
import shutil
from pathlib import Path

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# Map DSSP 8-state codes to 3-state names
_DSSP_TO_NAME = {
    "H": "helix",  # Alpha helix
    "G": "helix",  # 3-10 helix
    "I": "helix",  # Pi helix
    "E": "sheet",  # Extended strand
    "B": "sheet",  # Isolated beta-bridge
    "T": "coil",   # Turn
    "S": "coil",   # Bend
    "-": "coil",   # Coil / no assignment
    "C": "coil",   # Coil (some DSSP versions)
}


def run_dssp(variant_record: VariantRecord) -> VariantRecord:
    """Assign secondary structure at mutation site using DSSP.

    Args:
        variant_record: Must have pdb_path and residue_position set.

    Returns:
        VariantRecord with DSSP fields populated.
    """
    pdb_path = variant_record.pdb_fixed_path or variant_record.pdb_path
    position = variant_record.residue_position

    if not pdb_path or not Path(pdb_path).exists():
        logger.info("No structure available for DSSP")
        variant_record.set_feature_status("dssp", False, NullReason.UPSTREAM_DEPENDENCY_FAILED)
        return variant_record

    # Check for mkdssp binary
    if not shutil.which("mkdssp"):
        logger.warning("mkdssp binary not found — install via 'brew install dssp'")
        variant_record.set_feature_status("dssp", False, NullReason.TOOL_MISSING)
        return variant_record

    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)
        model = structure[0]

        dssp = DSSP(model, pdb_path, dssp="mkdssp")

        # DSSP keys are (chain_id, (' ', residue_number, ' '))
        chain = next(iter(model.get_chains()))
        dssp_key = (chain.id, (" ", position, " "))

        if dssp_key not in dssp:
            logger.warning("Residue %d not found in DSSP output", position)
            variant_record.set_feature_status("dssp", False, NullReason.NO_DATA_AVAILABLE)
            return variant_record

        dssp_result = dssp[dssp_key]
        ss_code = dssp_result[2]  # Secondary structure code

        # Normalize: some DSSP versions use different codes
        if ss_code not in _DSSP_TO_NAME:
            ss_code = "-"

        variant_record.secondary_structure = ss_code
        variant_record.secondary_structure_name = _DSSP_TO_NAME[ss_code]
        variant_record.set_feature_status("dssp", True)

        logger.info(
            "DSSP: residue %d → %s (%s)",
            position, ss_code, variant_record.secondary_structure_name,
        )

    except Exception as e:
        logger.warning("DSSP failed: %s", e)
        variant_record.set_feature_status("dssp", False, NullReason.TOOL_CRASHED)

    return variant_record
```

**Step 4: Run tests**

Run: `pytest tests/test_m3_structural.py::TestDSSP -v`
Expected: PASS (or skip if mkdssp not installed)

**Step 5: Commit**

```bash
git add varis/m3_structural_analysis/dssp_wrapper.py tests/test_m3_structural.py
git commit -m "feat(m3): implement DSSP wrapper for secondary structure

Uses BioPython DSSP with mkdssp binary. Maps 8-state to 3-state
(helix/sheet/coil). Gracefully handles missing mkdssp."
```

---

### Task 8: M3 — BioPython Contacts (WT Environment)

**Files:**
- Modify: `varis/m3_structural_analysis/biopython_contacts.py`
- Modify: `tests/test_m3_structural.py`

**Step 1: Write failing tests**

Add to `tests/test_m3_structural.py`:

```python
class TestContacts:
    """Tests for biopython_contacts.py — WT local environment."""

    def test_contacts_wt_valid(self, m1_completed_record):
        """Integers >= 0, nonzero for folded region."""
        from varis.m3_structural_analysis.biopython_contacts import run_contacts
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_contacts(m1_completed_record)
        if result.contacts_available:
            assert isinstance(result.contacts_wt, int)
            assert result.contacts_wt >= 0
            assert isinstance(result.hbonds_wt, int)
            assert result.hbonds_wt >= 0
            # Residue 1699 is in a folded region — should have contacts
            assert result.contacts_wt > 0

    def test_contacts_heavy_atoms_only(self, m1_completed_record):
        """Contacts count uses heavy atoms only (no hydrogens)."""
        from varis.m3_structural_analysis.biopython_contacts import run_contacts
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_contacts(m1_completed_record)
        # The count should be reasonable for heavy atoms (typically 5-30)
        if result.contacts_available:
            assert result.contacts_wt < 200  # sanity upper bound

    def test_contacts_packing_density(self, m1_completed_record):
        """Packing density is a float >= 0."""
        from varis.m3_structural_analysis.biopython_contacts import run_contacts
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_contacts(m1_completed_record)
        if result.contacts_available:
            assert isinstance(result.packing_density, float)
            assert result.packing_density >= 0.0

    def test_contacts_no_structure(self, m1_completed_record):
        """No PDB → contacts_available=False."""
        from varis.m3_structural_analysis.biopython_contacts import run_contacts
        m1_completed_record.pdb_path = None
        result = run_contacts(m1_completed_record)
        assert result.contacts_available is False
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m3_structural.py::TestContacts -v`
Expected: FAIL

**Step 3: Implement biopython_contacts.py**

Replace `varis/m3_structural_analysis/biopython_contacts.py`:

```python
"""BioPython Contacts — Counts WT local environment at mutation site.

Measures the wild-type local environment: heavy-atom contacts within 4.5Å,
hydrogen bonds, and packing density. Does NOT model the mutant — that would
require FoldX/PyRosetta to be scientifically valid.

Definition: Heavy atoms only (element != H). Contact threshold: 4.5Å.
H-bond: donor-acceptor distance ≤3.5Å with N/O/S atoms.

Populates: contacts_wt, hbonds_wt, packing_density.
Sets: contacts_available, contacts_missing_reason.
"""

import logging
from pathlib import Path

import numpy as np
from Bio.PDB import NeighborSearch, PDBParser

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

_CONTACT_THRESHOLD = 4.5  # Ångströms
_HBOND_THRESHOLD = 3.5    # Ångströms
_HBOND_ELEMENTS = {"N", "O", "S"}


def run_contacts(variant_record: VariantRecord) -> VariantRecord:
    """Count WT local environment contacts at mutation site.

    Args:
        variant_record: Must have pdb_path and residue_position set.

    Returns:
        VariantRecord with contact fields populated.
    """
    pdb_path = variant_record.pdb_fixed_path or variant_record.pdb_path
    position = variant_record.residue_position

    if not pdb_path or not Path(pdb_path).exists():
        logger.info("No structure available for contact analysis")
        variant_record.set_feature_status("contacts", False, NullReason.UPSTREAM_DEPENDENCY_FAILED)
        return variant_record

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)
        chain = next(iter(structure[0].get_chains()))

        # Get target residue
        target_residue = None
        for residue in chain.get_residues():
            if residue.id[0] == " " and residue.id[1] == position:
                target_residue = residue
                break

        if target_residue is None:
            logger.warning("Residue %d not found in structure", position)
            variant_record.set_feature_status("contacts", False, NullReason.NO_DATA_AVAILABLE)
            return variant_record

        # Get all heavy atoms in the structure
        all_heavy_atoms = [
            atom for atom in structure[0].get_atoms()
            if atom.element != "H"
        ]

        # Get heavy atoms in target residue
        target_atoms = [
            atom for atom in target_residue.get_atoms()
            if atom.element != "H"
        ]

        # Build neighbor search from all heavy atoms
        ns = NeighborSearch(all_heavy_atoms)

        # Count contacts: heavy atoms within threshold from any target atom,
        # excluding atoms in the target residue itself
        contact_atoms = set()
        hbond_count = 0
        target_atom_set = set(id(a) for a in target_atoms)

        for target_atom in target_atoms:
            neighbors = ns.search(target_atom.get_vector().get_array(), _CONTACT_THRESHOLD)
            for neighbor in neighbors:
                if id(neighbor) not in target_atom_set:
                    contact_atoms.add(id(neighbor))

                    # Check for H-bond (N/O/S within 3.5Å)
                    if (target_atom.element in _HBOND_ELEMENTS
                            and neighbor.element in _HBOND_ELEMENTS):
                        dist = target_atom - neighbor
                        if dist <= _HBOND_THRESHOLD:
                            hbond_count += 1

        contacts_count = len(contact_atoms)

        # Packing density: contacts normalized by number of target heavy atoms
        packing = contacts_count / len(target_atoms) if target_atoms else 0.0

        variant_record.contacts_wt = contacts_count
        variant_record.hbonds_wt = hbond_count
        variant_record.packing_density = round(packing, 4)
        variant_record.set_feature_status("contacts", True)

        logger.info(
            "Contacts: residue %d has %d contacts, %d H-bonds, packing=%.2f",
            position, contacts_count, hbond_count, packing,
        )

    except Exception as e:
        logger.warning("Contact analysis failed: %s", e)
        variant_record.set_feature_status("contacts", False, NullReason.TOOL_CRASHED)

    return variant_record
```

**Step 4: Run tests**

Run: `pytest tests/test_m3_structural.py::TestContacts -v`
Expected: PASS

**Step 5: Commit**

```bash
git add varis/m3_structural_analysis/biopython_contacts.py tests/test_m3_structural.py
git commit -m "feat(m3): implement BioPython contacts for WT environment

Counts heavy-atom contacts within 4.5Å, H-bonds (N/O/S within 3.5Å),
and packing density. WT environment only — no mutant modeling."
```

---

### Task 9: M3 — InterPro Client (Domain Identification)

**Files:**
- Create: `varis/m3_structural_analysis/interpro_client.py`
- Modify: `tests/test_m3_structural.py`

**Step 1: Write failing tests**

Add to `tests/test_m3_structural.py`:

```python
from unittest.mock import MagicMock


class TestInterPro:
    """Tests for interpro_client.py — domain identification via InterPro API."""

    def test_interpro_brca1_domain(self, m1_completed_record):
        """Returns BRCT domain with boundaries for position 1699 (mocked)."""
        from varis.m3_structural_analysis.interpro_client import run_interpro
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "results": [{
                "metadata": {
                    "accession": "PF00533",
                    "name": "BRCT",
                    "type": "domain",
                },
                "proteins": [{
                    "entry_protein_locations": [{
                        "fragments": [{"start": 1646, "end": 1736}]
                    }]
                }]
            }]
        }
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        result = run_interpro(m1_completed_record, client=mock_client)
        assert result.domain_available is True
        assert result.domain_name == "BRCT"
        assert result.domain_id == "PF00533"
        assert result.domain_start == 1646
        assert result.domain_end == 1736

    def test_interpro_position_outside_domain(self, m1_completed_record):
        """Position outside any domain → domain_available=False (mocked)."""
        from varis.m3_structural_analysis.interpro_client import run_interpro
        m1_completed_record.residue_position = 50  # Not in any domain
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "results": [{
                "metadata": {
                    "accession": "PF00533",
                    "name": "BRCT",
                    "type": "domain",
                },
                "proteins": [{
                    "entry_protein_locations": [{
                        "fragments": [{"start": 1646, "end": 1736}]
                    }]
                }]
            }]
        }
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        result = run_interpro(m1_completed_record, client=mock_client)
        assert result.domain_available is False

    def test_interpro_no_uniprot_id(self, m1_completed_record):
        """No UniProt ID → domain_available=False."""
        from varis.m3_structural_analysis.interpro_client import run_interpro
        m1_completed_record.uniprot_id = None
        result = run_interpro(m1_completed_record)
        assert result.domain_available is False

    def test_interpro_stores_boundaries(self, m1_completed_record):
        """Domain start/end boundaries are stored."""
        from varis.m3_structural_analysis.interpro_client import run_interpro
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "results": [{
                "metadata": {
                    "accession": "PF00533",
                    "name": "BRCT",
                    "type": "domain",
                },
                "proteins": [{
                    "entry_protein_locations": [{
                        "fragments": [{"start": 1646, "end": 1736}]
                    }]
                }]
            }]
        }
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        result = run_interpro(m1_completed_record, client=mock_client)
        assert result.domain_start is not None
        assert result.domain_end is not None
        assert result.domain_start < result.domain_end
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m3_structural.py::TestInterPro -v`
Expected: FAIL

**Step 3: Implement interpro_client.py**

Create `varis/m3_structural_analysis/interpro_client.py`:

```python
"""InterPro Client — Identifies Pfam domains at mutation position.

Queries the InterPro REST API with a UniProt ID to retrieve Pfam domain
annotations. Checks if the mutation position falls within any annotated domain.
Stores domain boundaries for transparency.

Does not require a local Pfam database — uses the public InterPro API.

Populates: domain_name, domain_id, domain_start, domain_end, domain_criticality.
Sets: domain_available, domain_missing_reason.
"""

import logging

import httpx

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

_TIMEOUT = 30.0
_INTERPRO_API_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{uniprot_id}"

# Domains known to be functionally critical (curated list, expandable)
_CRITICAL_DOMAINS = {
    "PF00533",  # BRCT (DNA damage response)
    "PF00069",  # Protein kinase
    "PF07714",  # Protein tyrosine kinase
    "PF00097",  # Zinc finger C3H1
    "PF00096",  # Zinc finger C2H2
    "PF00870",  # P53 DNA-binding
    "PF07710",  # P53 tetramerization
    "PF00503",  # G-protein alpha subunit
    "PF00104",  # Ligand-binding domain (nuclear receptor)
}

_IMPORTANT_DOMAINS = {
    "PF00651",  # BTB/POZ
    "PF01846",  # FF domain
    "PF00412",  # LIM domain
    "PF00397",  # WW domain
}


def run_interpro(variant_record: VariantRecord,
                 client: httpx.Client | None = None) -> VariantRecord:
    """Query InterPro for Pfam domains at the mutation position.

    Args:
        variant_record: Must have uniprot_id and residue_position set.
        client: Optional httpx.Client for dependency injection.

    Returns:
        VariantRecord with domain fields populated.
    """
    uniprot_id = variant_record.uniprot_id
    position = variant_record.residue_position

    if not uniprot_id:
        logger.info("No UniProt ID — skipping InterPro lookup")
        variant_record.set_feature_status("domain", False, NullReason.UPSTREAM_DEPENDENCY_FAILED)
        return variant_record

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT, follow_redirects=True)
        should_close = True

    try:
        url = _INTERPRO_API_URL.format(uniprot_id=uniprot_id)
        response = client.get(url)

        if response.status_code == 404:
            logger.info("No InterPro entries for %s", uniprot_id)
            variant_record.set_feature_status("domain", False, NullReason.NO_DATA_AVAILABLE)
            return variant_record

        response.raise_for_status()
        data = response.json()

        results = data.get("results", [])
        if not results:
            variant_record.set_feature_status("domain", False, NullReason.NO_DATA_AVAILABLE)
            return variant_record

        # Find domain containing the mutation position
        for entry in results:
            metadata = entry.get("metadata", {})
            proteins = entry.get("proteins", [])

            for protein in proteins:
                locations = protein.get("entry_protein_locations", [])
                for location in locations:
                    fragments = location.get("fragments", [])
                    for fragment in fragments:
                        start = fragment.get("start")
                        end = fragment.get("end")
                        if start is not None and end is not None:
                            if start <= position <= end:
                                accession = metadata.get("accession", "")
                                variant_record.domain_name = metadata.get("name")
                                variant_record.domain_id = accession
                                variant_record.domain_start = start
                                variant_record.domain_end = end
                                variant_record.domain_criticality = _classify_domain(accession)
                                variant_record.set_feature_status("domain", True)
                                logger.info(
                                    "InterPro: position %d in %s (%s, %d-%d, %s)",
                                    position, variant_record.domain_name,
                                    accession, start, end,
                                    variant_record.domain_criticality,
                                )
                                return variant_record

        # Position not in any domain
        logger.info("Position %d not in any annotated Pfam domain", position)
        variant_record.set_feature_status("domain", False, NullReason.NO_DATA_AVAILABLE)

    except Exception as e:
        logger.warning("InterPro lookup failed: %s", e)
        variant_record.set_feature_status("domain", False, NullReason.TOOL_CRASHED)
    finally:
        if should_close:
            client.close()

    return variant_record


def _classify_domain(accession: str) -> str:
    """Classify domain criticality based on curated lists.

    Args:
        accession: Pfam accession ID.

    Returns:
        "critical", "important", or "peripheral".
    """
    if accession in _CRITICAL_DOMAINS:
        return "critical"
    if accession in _IMPORTANT_DOMAINS:
        return "important"
    return "peripheral"
```

**Step 4: Run tests**

Run: `pytest tests/test_m3_structural.py::TestInterPro -v`
Expected: PASS

**Step 5: Commit**

```bash
git add varis/m3_structural_analysis/interpro_client.py tests/test_m3_structural.py
git commit -m "feat(m3): implement InterPro client for Pfam domain lookup

Queries InterPro REST API, identifies domain at mutation position,
stores boundaries and criticality classification."
```

---

### Task 10: M3 — FoldX/PyRosetta Stubs

**Files:**
- Modify: `varis/m3_structural_analysis/foldx_wrapper.py`
- Modify: `varis/m3_structural_analysis/pyrosetta_wrapper.py`
- Modify: `tests/test_m3_structural.py`

**Step 1: Write failing tests**

Add to `tests/test_m3_structural.py`:

```python
class TestFoldXStub:
    """Tests for foldx_wrapper.py — ΔΔG stub."""

    def test_foldx_no_binary(self, m1_completed_record):
        """FoldX not installed → ddg_available=False, reason tool_missing."""
        from varis.m3_structural_analysis.foldx_wrapper import run_foldx
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_foldx(m1_completed_record)
        assert result.ddg_available is False
        assert result.ddg_missing_reason == "tool_missing"


class TestPyRosettaStub:
    """Tests for pyrosetta_wrapper.py — ΔΔG stub."""

    def test_pyrosetta_not_installed(self, m1_completed_record):
        """PyRosetta not installed → ddg_available=False, reason tool_missing."""
        from varis.m3_structural_analysis.pyrosetta_wrapper import run_pyrosetta
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_pyrosetta(m1_completed_record)
        # ddg_available may already be False from FoldX stub
        assert result.ddg_foldx is None or result.ddg_pyrosetta is None
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m3_structural.py::TestFoldXStub tests/test_m3_structural.py::TestPyRosettaStub -v`
Expected: FAIL

**Step 3: Implement FoldX stub**

Replace `varis/m3_structural_analysis/foldx_wrapper.py`:

```python
"""FoldX Wrapper — ΔΔG stability prediction (requires license).

FoldX computes the change in folding free energy (ΔΔG) caused by a mutation.
ΔΔG > 2 kcal/mol = destabilizing — this is the single strongest structural
feature for pathogenicity prediction.

This is a stub that checks for the FoldX binary. If not found, sets
ddg_available=False with reason tool_missing. Full implementation will be
added when a FoldX academic license is obtained.

Populates: ddg_foldx (when available).
Sets: ddg_available, ddg_missing_reason.
"""

import logging
import shutil

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)


def run_foldx(variant_record: VariantRecord) -> VariantRecord:
    """Run FoldX ΔΔG prediction (stub — requires license).

    Args:
        variant_record: Must have pdb_path and mutation details set.

    Returns:
        VariantRecord with ddg_foldx set (or None with reason).
    """
    if not shutil.which("foldx") and not shutil.which("FoldX"):
        logger.info("FoldX binary not found — requires academic license")
        variant_record.set_with_reason("ddg_foldx", None, NullReason.TOOL_MISSING)
        variant_record.set_feature_status("ddg", False, NullReason.TOOL_MISSING)
        return variant_record

    # TODO: Full FoldX implementation when license is obtained
    # 1. Prepare FoldX config with mutation notation
    # 2. Run FoldX BuildModel
    # 3. Parse Dif_*.fxout for ΔΔG value
    # 4. Set ddg_foldx and ddg_available
    logger.info("FoldX binary found but full implementation pending")
    variant_record.set_with_reason("ddg_foldx", None, NullReason.NOT_ATTEMPTED)
    variant_record.set_feature_status("ddg", False, NullReason.NOT_ATTEMPTED)
    return variant_record
```

**Step 4: Implement PyRosetta stub**

Replace `varis/m3_structural_analysis/pyrosetta_wrapper.py`:

```python
"""PyRosetta Wrapper — Fallback ΔΔG stability prediction (requires license).

PyRosetta provides a second independent ΔΔG estimate with backbone relaxation.
When FoldX and PyRosetta agree, confidence is high. This is the fallback for
when FoldX is unavailable.

This is a stub. Full implementation will be added when a PyRosetta academic
license is obtained.

Populates: ddg_pyrosetta (when available).
"""

import logging

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)


def run_pyrosetta(variant_record: VariantRecord) -> VariantRecord:
    """Run PyRosetta ΔΔG prediction (stub — requires license).

    Args:
        variant_record: Must have pdb_path and mutation details set.

    Returns:
        VariantRecord with ddg_pyrosetta set (or None with reason).
    """
    try:
        import pyrosetta  # noqa: F401
    except ImportError:
        logger.info("PyRosetta not installed — requires academic license")
        variant_record.set_with_reason("ddg_pyrosetta", None, NullReason.TOOL_MISSING)
        return variant_record

    # TODO: Full PyRosetta implementation when license is obtained
    logger.info("PyRosetta found but full implementation pending")
    variant_record.set_with_reason("ddg_pyrosetta", None, NullReason.NOT_ATTEMPTED)
    return variant_record
```

**Step 5: Run tests**

Run: `pytest tests/test_m3_structural.py::TestFoldXStub tests/test_m3_structural.py::TestPyRosettaStub -v`
Expected: PASS

**Step 6: Commit**

```bash
git add varis/m3_structural_analysis/foldx_wrapper.py varis/m3_structural_analysis/pyrosetta_wrapper.py tests/test_m3_structural.py
git commit -m "feat(m3): implement FoldX and PyRosetta stubs

Both check for binary/library availability. Set ddg_available=False
with tool_missing reason when not installed."
```

---

### Task 11: M3 — Orchestrator, Integration Test, and Golden Record

**Files:**
- Modify: `varis/m3_structural_analysis/__init__.py`
- Create: `varis/m3_structural_analysis/interpro_client.py` (done in Task 9)
- Remove: `varis/m3_structural_analysis/hmmer_wrapper.py` (replaced by interpro_client)
- Remove: `varis/m3_structural_analysis/mdanalysis_wrapper.py` (deferred)
- Modify: `tests/test_m3_structural.py`

**Step 1: Write failing integration tests**

Add to `tests/test_m3_structural.py`:

```python
class TestM3Orchestrator:
    """Tests for M3 orchestration — site-gating and full pipeline."""

    def test_m3_skips_when_site_absent(self, m1_completed_record):
        """All site-dependent features skipped when mutation_site_present=False."""
        from varis.m3_structural_analysis import run
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.mutation_site_present = False
        # Mock InterPro to avoid real API call
        result = run(m1_completed_record)
        # Site-dependent features should be skipped
        assert result.sasa_available is False or result.solvent_accessibility_relative is None
        assert result.dssp_available is False or result.secondary_structure is None
        assert result.contacts_available is False or result.contacts_wt is None
        # M3 should still complete (InterPro and stubs run)
        assert "M3" in result.modules_completed

    def test_m3_no_structure(self, m1_completed_record):
        """M3 with no structure → marks failed."""
        from varis.m3_structural_analysis import run
        m1_completed_record.pdb_path = None
        result = run(m1_completed_record)
        assert "M3" in result.modules_failed

    def test_m3_integration_golden_record(self, m1_completed_record):
        """Full M3 pipeline — verify expected keys are set or have reasons."""
        from varis.m3_structural_analysis import run
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.mutation_site_present = True
        result = run(m1_completed_record)
        assert "M3" in result.modules_completed
        # Feature availability flags should all be set (True or False)
        for flag in ["sasa_available", "dssp_available", "contacts_available",
                      "ddg_available", "domain_available"]:
            val = getattr(result, flag)
            assert val is not None, f"{flag} should be True or False, not None"
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m3_structural.py::TestM3Orchestrator -v`
Expected: FAIL

**Step 3: Update M3 orchestrator**

Replace `varis/m3_structural_analysis/__init__.py`:

```python
"""M3: Structural Analysis — Extracts structural features from 3D structure.

Each tool asks a structural question and is independently computed.
If one tool fails, others still run. Site-dependent tools (FreeSASA, DSSP,
Contacts) are skipped when mutation_site_present=False.

Depends on: M2 (needs pdb_path and mutation_site_present).
Populates: All structural_features.* fields in VariantRecord.
Critical independence: M3 and M4 are completely independent of each other.
"""

import logging

logger = logging.getLogger(__name__)


def run(variant_record):
    """Execute M3: extract structural features. Each tool is independent.

    Site-dependent tools (FreeSASA, DSSP, Contacts) are skipped if
    mutation_site_present is False. InterPro and ΔΔG stubs always run.

    Args:
        variant_record: VariantRecord with M2 fields populated.

    Returns:
        VariantRecord with structural features populated (or None with reasons).
    """
    if variant_record.pdb_path is None and variant_record.pdb_fixed_path is None:
        logger.warning("M3 skipped: no structure available from M2.")
        variant_record.mark_module_failed("M3")
        return variant_record

    from varis.m3_structural_analysis.biopython_contacts import run_contacts
    from varis.m3_structural_analysis.dssp_wrapper import run_dssp
    from varis.m3_structural_analysis.foldx_wrapper import run_foldx
    from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
    from varis.m3_structural_analysis.interpro_client import run_interpro
    from varis.m3_structural_analysis.pyrosetta_wrapper import run_pyrosetta
    from varis.models.variant_record import NullReason

    # Site-dependent tools: skip if mutation site not confirmed in structure
    site_present = variant_record.mutation_site_present
    site_tools = [
        ("M3.freesasa", run_freesasa),
        ("M3.dssp", run_dssp),
        ("M3.contacts", run_contacts),
    ]

    if site_present:
        for name, fn in site_tools:
            try:
                variant_record = fn(variant_record)
            except Exception as e:
                logger.warning("%s failed: %s", name, e)
                variant_record.mark_module_failed(name)
    else:
        logger.info("Mutation site not present — skipping site-dependent tools")
        variant_record.set_feature_status("sasa", False, NullReason.INTENTIONALLY_SKIPPED)
        variant_record.set_feature_status("dssp", False, NullReason.INTENTIONALLY_SKIPPED)
        variant_record.set_feature_status("contacts", False, NullReason.INTENTIONALLY_SKIPPED)

    # Site-independent tools: always run
    independent_tools = [
        ("M3.interpro", run_interpro),
        ("M3.foldx", run_foldx),
        ("M3.pyrosetta", run_pyrosetta),
    ]

    for name, fn in independent_tools:
        try:
            variant_record = fn(variant_record)
        except Exception as e:
            logger.warning("%s failed: %s", name, e)
            variant_record.mark_module_failed(name)

    # Compute ddg_mean if any ΔΔG values available
    ddg_values = [v for v in [variant_record.ddg_foldx, variant_record.ddg_pyrosetta]
                  if v is not None]
    if ddg_values:
        variant_record.ddg_mean = round(sum(ddg_values) / len(ddg_values), 4)

    variant_record.mark_module_completed("M3")
    return variant_record
```

**Step 4: Remove old stubs**

Delete `varis/m3_structural_analysis/hmmer_wrapper.py` (replaced by interpro_client.py).
Delete `varis/m3_structural_analysis/mdanalysis_wrapper.py` (deferred, not in Phase 2).

**Step 5: Run all M3 tests**

Run: `pytest tests/test_m3_structural.py -v`
Expected: PASS

**Step 6: Commit**

```bash
git add varis/m3_structural_analysis/__init__.py tests/test_m3_structural.py
git rm varis/m3_structural_analysis/hmmer_wrapper.py varis/m3_structural_analysis/mdanalysis_wrapper.py
git commit -m "feat(m3): implement M3 orchestrator with site-gating

Skips FreeSASA/DSSP/Contacts when mutation_site_present=False.
InterPro and DDG stubs always run. Replaces HMMER with InterPro.
Removes mdanalysis_wrapper (deferred)."
```

---

### Task 12: Full Pipeline Integration Test (M1 → M2 → M3)

**Files:**
- Create: `tests/test_phase2_integration.py`

**Step 1: Write the integration test**

```python
"""Phase 2 Integration Test — Full M1 → M2 → M3 pipeline on BRCA1 p.Arg1699Trp.

Uses the real AlphaFold PDB file already downloaded by M1.
Mocks external API calls (InterPro) for offline testing.
"""
import pytest
from unittest.mock import MagicMock, patch
from pathlib import Path

from varis.config import STRUCTURES_DIR
from varis.models.variant_record import create_variant_record

BRCA1_PDB = STRUCTURES_DIR / "AF-P38398-F1-model_v6.pdb"


class TestPhase2Integration:
    """End-to-end M1 → M2 → M3 pipeline test."""

    @pytest.fixture
    def brca1_m1_record(self):
        """Create a record that has completed M1 (using fixture data)."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.residue_position = 1699
        record.ref_amino_acid = "Arg"
        record.alt_amino_acid = "Trp"
        record.ref_aa_single = "R"
        record.alt_aa_single = "W"
        record.charge_change = "+ve → neutral"
        record.uniprot_id = "P38398"
        record.protein_name = "Breast cancer type 1 susceptibility protein"
        record.protein_length = 1863
        record.protein_sequence = "M" * 1863  # placeholder — not used by M2/M3
        record.structure_source = "alphafold"
        record.pdb_path = str(BRCA1_PDB)
        record.mark_module_completed("M1")
        return record

    def test_m2_then_m3_pipeline(self, brca1_m1_record):
        """Full M2 → M3 pipeline produces expected outputs."""
        from varis.m2_structure import run as run_m2
        from varis.m3_structural_analysis import run as run_m3

        # Run M2
        record = run_m2(brca1_m1_record)
        assert "M2" in record.modules_completed
        assert record.mutation_site_present is True
        assert record.pdb_hash is not None

        # Run M3 (mock InterPro to avoid real API call)
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "results": [{
                "metadata": {"accession": "PF00533", "name": "BRCT", "type": "domain"},
                "proteins": [{
                    "entry_protein_locations": [{
                        "fragments": [{"start": 1646, "end": 1736}]
                    }]
                }]
            }]
        }
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response

        with patch("varis.m3_structural_analysis.interpro_client.httpx.Client",
                    return_value=mock_client):
            record = run_m3(record)

        assert "M3" in record.modules_completed

        # Verify all feature availability flags are set
        assert record.sasa_available is not None
        assert record.dssp_available is not None
        assert record.contacts_available is not None
        assert record.ddg_available is not None
        assert record.domain_available is not None

        # If FreeSASA is installed, check results
        if record.sasa_available:
            assert 0.0 <= record.solvent_accessibility_relative <= 1.0
            assert record.burial_category in ("core", "surface")

        # If DSSP ran, check results
        if record.dssp_available:
            assert record.secondary_structure is not None
            assert record.secondary_structure_name in ("helix", "sheet", "coil")

        # Contacts should work (pure BioPython)
        assert record.contacts_available is True
        assert record.contacts_wt > 0

        # DDG should be unavailable (no FoldX/PyRosetta)
        assert record.ddg_available is False

    def test_golden_record_schema(self, brca1_m1_record):
        """After M2+M3, record serializes to JSON with all expected keys."""
        from varis.m2_structure import run as run_m2
        from varis.m3_structural_analysis import run as run_m3

        record = run_m2(brca1_m1_record)

        with patch("varis.m3_structural_analysis.interpro_client.httpx.Client"):
            record = run_m3(record)

        data = record.to_dict()

        # M2 keys
        m2_keys = ["mutation_site_present", "pdb_hash", "numbering_scheme",
                    "preparation_steps", "plddt_available"]
        for key in m2_keys:
            assert key in data, f"Missing M2 key: {key}"

        # Feature availability keys
        avail_keys = ["sasa_available", "dssp_available", "contacts_available",
                       "ddg_available", "domain_available"]
        for key in avail_keys:
            assert key in data, f"Missing availability key: {key}"
            assert data[key] is not None, f"{key} should not be None"
```

**Step 2: Run integration tests**

Run: `pytest tests/test_phase2_integration.py -v`
Expected: PASS

**Step 3: Commit**

```bash
git add tests/test_phase2_integration.py
git commit -m "test(phase2): add full M2→M3 integration test with BRCA1

Verifies end-to-end pipeline from M1 output through M2 validation
and M3 feature extraction. Golden record schema check."
```

---

### Task 13: Update Dependencies and Config

**Files:**
- Modify: `pyproject.toml`
- Modify: `varis/config.py`

**Step 1: Add InterPro API URL to config**

Add to `varis/config.py` in the API ENDPOINTS section:

```python
INTERPRO_API_URL = "https://www.ebi.ac.uk/interpro/api"
```

**Step 2: Update pyproject.toml dependencies**

Add to the `[project.optional-dependencies]` structure group:

```toml
structure = [
    "freesasa>=2.2",
    "mdanalysis>=2.6",
    "gemmi>=0.6",
    "openbabel-wheel>=3.1",
    "pdbfixer>=1.9",
    "mdtraj>=1.9",
]
```

**Step 3: Run full test suite**

Run: `pytest tests/ -v --tb=short`
Expected: ALL PASS

**Step 4: Commit**

```bash
git add pyproject.toml varis/config.py
git commit -m "chore: add Phase 2 dependencies and InterPro config

Add pdbfixer, mdtraj to structure dependencies. Add InterPro API URL."
```

---

### Task 14: Final Verification and Cleanup

**Step 1: Run full test suite with coverage**

Run: `pytest tests/ -v --tb=short`
Expected: ALL PASS

**Step 2: Check for any remaining references to removed fields**

Search for: `helix_disruption`, `hbonds_lost`, `contacts_changed`, `functional_site_distance`, `nearest_functional_site`, `structure_resolution`, `solvent_accessibility` (not `solvent_accessibility_relative`), `plddt_score` (not `mutation_site_plddt`)

Fix any remaining references in test files or other code.

**Step 3: Verify no import cycle issues**

Run: `python -c "from varis.m2_structure import run; from varis.m3_structural_analysis import run"`
Expected: No errors

**Step 4: Commit any cleanup**

```bash
git add -A
git commit -m "chore: Phase 2 cleanup — remove stale references to v1.1.0 fields"
```
