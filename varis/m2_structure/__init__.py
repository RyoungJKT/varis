"""M2: Structure Engine — Obtains and prepares 3D protein structures.

Depends on: M1 (needs uniprot_id or protein_sequence)
Populates: pdb_path, pdb_fixed_path, plddt_score, plddt_mean, structure_source
"""
import logging
logger = logging.getLogger(__name__)

def run(variant_record):
    """Execute M2: obtain and prepare 3D structure.
    Tries AlphaFold DB first, ESMFold fallback, then PDBFixer repair."""
    from varis.m2_structure.alphafold_retriever import retrieve_alphafold
    from varis.m2_structure.esmfold_predictor import predict_esmfold
    from varis.m2_structure.pdb_fixer import fix_structure
    from varis.m2_structure.structure_utils import extract_plddt

    try:
        variant_record = retrieve_alphafold(variant_record)
    except Exception as e:
        logger.warning(f"AlphaFold retrieval failed: {e}")
        variant_record.mark_module_failed("M2.alphafold")

    if variant_record.pdb_path is None:
        try:
            variant_record = predict_esmfold(variant_record)
        except Exception as e:
            logger.warning(f"ESMFold fallback failed: {e}")
            variant_record.mark_module_failed("M2.esmfold")

    if variant_record.pdb_path is not None:
        try: variant_record = fix_structure(variant_record)
        except Exception as e:
            logger.warning(f"PDBFixer failed: {e}"); variant_record.mark_module_failed("M2.pdbfixer")
        try: variant_record = extract_plddt(variant_record)
        except Exception as e:
            logger.warning(f"pLDDT extraction failed: {e}"); variant_record.mark_module_failed("M2.plddt")
        variant_record.mark_module_completed("M2")
    else:
        logger.warning("No structure obtained. M3 structural analysis will be skipped.")
        variant_record.mark_module_failed("M2")
    return variant_record
