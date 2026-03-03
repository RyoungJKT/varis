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
