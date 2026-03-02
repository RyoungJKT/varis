"""Varis Pipeline — Orchestrates the full investigation workflow.

Runs M1 through M5 in sequence (M6 and M7 are separate services).
Each module is wrapped in try/except — no single failure stops the pipeline.
The VariantRecord accumulates results from each module.
"""

import logging
import time
from varis.models.variant_record import VariantRecord, create_variant_record
from varis.config import PIPELINE_VERSION
from varis import m1_ingestion
from varis import m2_structure
from varis import m3_structural_analysis
from varis import m4_conservation
from varis import m5_scoring

logger = logging.getLogger(__name__)


def investigate(gene: str, hgvs_protein: str) -> VariantRecord:
    """Run the full Varis investigation pipeline on a single variant.

    This is the main entry point. Creates a VariantRecord and passes it
    through M1 → M2 → M3 → M4 → M5, collecting structural evidence
    at each stage. If any module fails, the pipeline continues with
    whatever data it has.

    Args:
        gene: Gene symbol, e.g., "BRCA1"
        hgvs_protein: Protein-level HGVS, e.g., "p.Arg1699Trp"

    Returns:
        Completed VariantRecord with all available evidence.
    """
    start = time.time()
    record = create_variant_record(gene, hgvs_protein)
    record.pipeline_version = PIPELINE_VERSION

    logger.info(f"Starting investigation: {gene} {hgvs_protein}")

    # M1: Data Ingestion (always runs)
    record = _run_module("M1", m1_ingestion.run, record)

    # M2: Structure Engine (needs M1)
    record = _run_module("M2", m2_structure.run, record)

    # M3: Structural Analysis (needs M2, independent of M4)
    record = _run_module("M3", m3_structural_analysis.run, record)

    # M4: Conservation (needs M1, independent of M2/M3)
    record = _run_module("M4", m4_conservation.run, record)

    # M5: ML Scoring (takes whatever features are available)
    record = _run_module("M5", m5_scoring.run, record)

    record.processing_time_seconds = round(time.time() - start, 2)
    logger.info(
        f"Investigation complete: {gene} {hgvs_protein} | "
        f"Features: {record.count_available_features()}/15 | "
        f"Score: {record.score_ensemble} | "
        f"Time: {record.processing_time_seconds}s"
    )
    return record


def _run_module(name: str, module_fn, record: VariantRecord) -> VariantRecord:
    """Run a module with top-level exception handling.

    If the module itself crashes (beyond its internal error handling),
    catch it here so the pipeline continues.
    """
    try:
        record = module_fn(record)
    except Exception as e:
        logger.error(f"Module {name} crashed unexpectedly: {e}")
        record.mark_module_failed(name)
    return record
