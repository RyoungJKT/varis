"""M4: Conservation Engine — Calculates evolutionary conservation.

CRITICAL INDEPENDENCE: M4 depends on M1 (needs protein sequence), NOT on M2 or M3.
If all structural analysis fails, conservation still works.
If conservation fails, structural analysis still works.

Depends on: M1 (needs protein_sequence)
Populates: conservation_score, conservation_method, num_orthologs, position_entropy, conserved_across_mammals
"""
import logging
logger = logging.getLogger(__name__)

def run(variant_record):
    """Execute M4: evolutionary conservation analysis. Independent of M2/M3."""
    from varis.m4_conservation.blast_client import run_blast
    from varis.m4_conservation.clustal_client import run_alignment
    from varis.m4_conservation.conservation_scorer import score_conservation
    from varis.m4_conservation.consurf_fallback import fetch_consurf

    # Step 1: Find homologs via BLAST
    homologs = None
    try:
        variant_record, homologs = run_blast(variant_record)
    except Exception as e:
        logger.warning(f"BLAST failed: {e}")
        variant_record.mark_module_failed("M4.blast")

    # Step 2: Multiple sequence alignment
    alignment = None
    if homologs is not None:
        try:
            variant_record, alignment = run_alignment(variant_record, homologs)
        except Exception as e:
            logger.warning(f"Alignment failed: {e}")
            variant_record.mark_module_failed("M4.alignment")

    # Step 3: Score conservation from alignment
    if alignment is not None:
        try:
            variant_record = score_conservation(variant_record, alignment)
        except Exception as e:
            logger.warning(f"Conservation scoring failed: {e}")
            variant_record.mark_module_failed("M4.scoring")

    # Fallback: ConSurf pre-computed scores if BLAST/alignment failed
    if variant_record.conservation_score is None:
        try:
            variant_record = fetch_consurf(variant_record)
        except Exception as e:
            logger.warning(f"ConSurf fallback failed: {e}")
            variant_record.mark_module_failed("M4.consurf")

    if variant_record.conservation_score is not None:
        variant_record.mark_module_completed("M4")
    else:
        variant_record.mark_module_failed("M4")
    return variant_record
