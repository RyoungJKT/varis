"""M1: Data Ingestion — Parses variants and retrieves data from public databases.

This is the entry point for every investigation. It takes a gene symbol and
HGVS protein notation, and populates the VariantRecord with identifiers,
protein data, population frequencies, and external scores.

Depends on: Nothing (this is the entry point)
Populates: identifiers, variant details, protein, population, external scores
"""

from varis.m1_ingestion.hgvs_parser import parse_hgvs
from varis.m1_ingestion.variant_normalizer import normalize_variant
from varis.m1_ingestion.clinvar_client import fetch_clinvar
from varis.m1_ingestion.gnomad_client import fetch_gnomad
from varis.m1_ingestion.uniprot_client import fetch_uniprot
from varis.m1_ingestion.alphafold_client import fetch_alphafold_structure
from varis.m1_ingestion.alphamissense_client import fetch_alphamissense

import logging

logger = logging.getLogger(__name__)


def run(variant_record: dict) -> dict:
    """Execute all M1 sub-modules in sequence.

    Each sub-module is wrapped in try/except. If any fails,
    the others still run. Fields from failed sub-modules remain None
    with reason codes recorded in null_reasons.

    Order matters: hgvs_parser → uniprot → normalizer → clinvar → gnomad → alphafold → alphamissense.
    UniProt runs before normalizer because normalizer validates ref AA against UniProt sequence.
    ClinVar runs before gnomAD because gnomAD needs ClinVar's genomic coordinates.

    Args:
        variant_record: VariantRecord with gene_symbol and hgvs_protein set.

    Returns:
        VariantRecord with M1 fields populated (or None with reasons).
    """
    steps = [
        ("M1.hgvs_parser", parse_hgvs),
        ("M1.uniprot", fetch_uniprot),
        ("M1.normalizer", normalize_variant),
        ("M1.clinvar", fetch_clinvar),
        ("M1.gnomad", fetch_gnomad),
        ("M1.alphafold", fetch_alphafold_structure),
        ("M1.alphamissense", fetch_alphamissense),
    ]

    for step_name, step_fn in steps:
        try:
            variant_record = step_fn(variant_record)
        except Exception as e:
            logger.warning(f"{step_name} failed: {e}")
            variant_record.mark_module_failed(step_name)

    variant_record.mark_module_completed("M1")
    return variant_record
