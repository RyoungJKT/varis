"""Variant Normalizer — Coordinate mapping and notation disambiguation.

This is the hardest part of M1 and one of the most common sources of silent errors
in bioinformatics pipelines. The same mutation can be expressed differently across
databases, and residue numbering varies between transcripts, isoforms, and structures.

Problems this module solves:
  - Same variant, different notations: p.Arg1699Trp vs p.R1699W vs c.5095C>T
  - Transcript differences: NM_007294.4 vs NM_007297.4 give different positions
  - Isoform numbering: UniProt canonical position != AlphaFold structure position
  - ClinVar uses genomic coordinates; UniProt uses protein coordinates
  - AlphaFold structures may number from 1 based on a different isoform

Strategy:
  1. Parse raw input (any notation format)
  2. Map to MANE Select canonical transcript (NCBI standard)
  3. Map protein position to UniProt canonical isoform
  4. Map protein position to structure residue number (via SIFTS or BLAST alignment)
  5. Validate: check that ref amino acid at mapped position matches expected residue
  6. Record confidence level and any warnings

Known-bad test cases (build tests for these):
  - BRCA1: multiple transcripts with different exon structures
  - TP53: common isoform confusion between p53alpha and p53beta
  - CFTR: c.1521_1523delCTT (F508del) — deletion, not missense, should be rejected
  - Variants where ClinVar and UniProt disagree on residue numbering

Priority: Build this BEFORE any structural analysis. A wrong position = silently
wrong ΔΔG, wrong SASA, wrong conservation score. Everything downstream is garbage.
"""

import logging
from varis.models.variant_record import VariantRecord, NullReason

logger = logging.getLogger(__name__)


def normalize_variant(variant_record: VariantRecord) -> VariantRecord:
    """Normalize variant notation and map coordinates across databases.

    This runs AFTER hgvs_parser (which extracts basic fields) and BEFORE
    any database lookups. It resolves the canonical transcript, maps
    positions across coordinate systems, and validates consistency.

    Populates:
      - canonical_transcript, canonical_protein
      - uniprot_residue_position, structure_residue_position
      - coordinate_mapping_method, coordinate_mapping_confidence
      - normalization_warnings
      - input_notation_normalized

    Args:
        variant_record: Must have gene_symbol and hgvs_protein set.

    Returns:
        VariantRecord with normalization fields populated.
    """
    pass


def _resolve_canonical_transcript(gene_symbol: str) -> tuple[str, str] | None:
    """Look up the MANE Select transcript for a gene.

    MANE (Matched Annotation from NCBI and EBI) Select provides a single
    canonical transcript agreed on by both RefSeq and Ensembl for each gene.
    This eliminates transcript ambiguity.

    Returns:
        Tuple of (transcript_id, protein_id) or None if not found.
    """
    pass


def _map_to_uniprot_position(residue_position: int, transcript_id: str,
                              uniprot_id: str) -> tuple[int, str, str]:
    """Map a transcript-based residue position to UniProt canonical coordinates.

    Strategy:
    1. Try direct mapping (positions often match for canonical isoform)
    2. If mismatch, use SIFTS database (PDBe structure integration)
    3. If SIFTS unavailable, use pairwise BLAST alignment

    Returns:
        Tuple of (uniprot_position, method, confidence).
        method: "direct" / "sifts" / "blast_alignment"
        confidence: "exact" / "high" / "low" / "failed"
    """
    pass


def _map_to_structure_position(uniprot_position: int, uniprot_id: str,
                                pdb_path: str = None) -> tuple[int, str, str]:
    """Map UniProt position to the 3D structure residue number.

    AlphaFold structures are usually numbered by UniProt position (direct),
    but experimental PDB structures may have different numbering, missing
    residues, or insertion codes.

    Returns:
        Tuple of (structure_position, method, confidence).
    """
    pass


def _validate_residue_match(protein_sequence: str, position: int,
                             expected_aa: str) -> bool:
    """Verify that the amino acid at the mapped position matches expectations.

    This is the critical safety check. If the reference amino acid at the
    mapped position doesn't match what we expect, the mapping is wrong
    and all downstream analysis would be on the wrong residue.

    Returns:
        True if match confirmed, False if mismatch detected.
    """
    pass
