"""M3: Structural Analysis — Extracts structural features from 3D structure.

This is the scientific core of Varis. Every tool asks a structural question.
Build in priority order: FreeSASA (easy) → DSSP → BioPython → FoldX → PyRosetta (hard).

Depends on: M2 (needs pdb_path). If M2 failed, M3 is skipped entirely.
Populates: All structural_features.* fields in VariantRecord.
Critical independence: M3 and M4 are completely independent of each other.
"""
import logging
logger = logging.getLogger(__name__)

def run(variant_record):
    """Execute M3: extract structural features. Each tool is independent."""
    if variant_record.pdb_path is None and variant_record.pdb_fixed_path is None:
        logger.warning("M3 skipped: no structure available from M2.")
        variant_record.mark_module_failed("M3")
        return variant_record

    from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
    from varis.m3_structural_analysis.dssp_wrapper import run_dssp
    from varis.m3_structural_analysis.biopython_contacts import run_contacts
    from varis.m3_structural_analysis.foldx_wrapper import run_foldx
    from varis.m3_structural_analysis.pyrosetta_wrapper import run_pyrosetta
    from varis.m3_structural_analysis.interpro_client import run_interpro
    from varis.m3_structural_analysis.mdanalysis_wrapper import run_mdanalysis

    tools = [
        ("M3.freesasa", run_freesasa),
        ("M3.dssp", run_dssp),
        ("M3.contacts", run_contacts),
        ("M3.foldx", run_foldx),
        ("M3.pyrosetta", run_pyrosetta),
        ("M3.interpro", run_interpro),
        ("M3.mdanalysis", run_mdanalysis),
    ]
    for name, fn in tools:
        try: variant_record = fn(variant_record)
        except Exception as e:
            logger.warning(f"{name} failed: {e}"); variant_record.mark_module_failed(name)

    variant_record.mark_module_completed("M3")
    return variant_record
