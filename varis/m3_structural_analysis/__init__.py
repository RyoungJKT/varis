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
    mutation_site_present is False. InterPro and DDG stubs always run.

    Args:
        variant_record: VariantRecord with M2 fields populated.

    Returns:
        VariantRecord with structural features populated (or None with reasons).
    """
    if variant_record.pdb_path is None and variant_record.pdb_fixed_path is None:
        logger.warning("M3 skipped: no structure available from M2.")
        variant_record.mark_module_failed("M3")
        return variant_record

    from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
    from varis.m3_structural_analysis.dssp_wrapper import run_dssp
    from varis.m3_structural_analysis.biopython_contacts import run_contacts
    from varis.m3_structural_analysis.interpro_client import run_interpro
    from varis.m3_structural_analysis.evoef2_wrapper import run_evoef2
    from varis.m3_structural_analysis.foldx_wrapper import run_foldx
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
        ("M3.evoef2", run_evoef2),
        ("M3.foldx", run_foldx),
        ("M3.pyrosetta", run_pyrosetta),
    ]

    for name, fn in independent_tools:
        try:
            variant_record = fn(variant_record)
        except Exception as e:
            logger.warning("%s failed: %s", name, e)
            variant_record.mark_module_failed(name)

    # Compute ddg_mean if any DDG values available
    ddg_values = [v for v in [
        variant_record.ddg_evoef2,
        variant_record.ddg_foldx,
        variant_record.ddg_pyrosetta,
    ] if v is not None]
    if ddg_values:
        variant_record.ddg_mean = round(sum(ddg_values) / len(ddg_values), 4)

    # Set ddg_available based on whether any DDG tool succeeded
    if variant_record.ddg_mean is not None:
        variant_record.set_feature_status("ddg", True)
    else:
        variant_record.set_feature_status("ddg", False, NullReason.TOOL_MISSING)

    variant_record.mark_module_completed("M3")
    return variant_record
