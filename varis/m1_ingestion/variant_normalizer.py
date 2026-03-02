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
wrong DDG, wrong SASA, wrong conservation score. Everything downstream is garbage.

Phase 1 Implementation:
  Simple direct mapping with ref AA validation against UniProt sequence.
  Transcript resolution and cross-database alignment stubs return defaults.
"""

import logging
from typing import Optional

from varis.models.variant_record import VariantRecord, NullReason

logger = logging.getLogger(__name__)


def normalize_variant(variant_record: VariantRecord) -> VariantRecord:
    """Normalize variant notation and map coordinates across databases.

    This runs AFTER hgvs_parser (which extracts basic fields) and BEFORE
    any database lookups. It resolves the canonical transcript, maps
    positions across coordinate systems, and validates consistency.

    Phase 1 performs direct mapping: assumes the parsed residue_position from
    HGVS notation corresponds directly to the UniProt canonical isoform position,
    then validates by checking the reference amino acid against the UniProt
    protein sequence.

    Populates:
      - coordinate_mapping_confidence: "high" on match, "failed" on mismatch
      - coordinate_mapping_method: "direct" when mapping succeeds
      - uniprot_residue_position: the validated position in UniProt coordinates
      - normalization_warnings: list of any issues encountered

    Args:
        variant_record: Must have residue_position and ref_aa_single set by
            the HGVS parser, and protein_sequence set by the UniProt fetcher.

    Returns:
        VariantRecord with normalization fields populated. Never raises;
        on failure, sets coordinate_mapping_confidence to None or "failed"
        with appropriate reason codes.
    """
    position = variant_record.residue_position
    ref_aa = variant_record.ref_aa_single
    sequence = variant_record.protein_sequence

    # --- Guard: position or ref_aa missing (upstream parser failed) -----------
    if position is None or ref_aa is None:
        logger.warning(
            "Cannot normalize variant %s: residue_position=%s, ref_aa_single=%s "
            "(upstream HGVS parsing likely failed)",
            variant_record.variant_id,
            position,
            ref_aa,
        )
        variant_record.set_with_reason(
            "coordinate_mapping_confidence",
            None,
            NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    # --- Guard: protein sequence missing (upstream UniProt fetch failed) ------
    if sequence is None:
        logger.warning(
            "Cannot normalize variant %s: protein_sequence is None "
            "(upstream UniProt fetch likely failed)",
            variant_record.variant_id,
        )
        variant_record.set_with_reason(
            "coordinate_mapping_confidence",
            None,
            NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    # --- Guard: position out of range -----------------------------------------
    if position < 1 or position > len(sequence):
        warning_msg = (
            f"Residue position {position} is out of range for protein sequence "
            f"of length {len(sequence)} (variant {variant_record.variant_id})"
        )
        logger.warning(warning_msg)
        variant_record.coordinate_mapping_confidence = "failed"
        if variant_record.normalization_warnings is None:
            variant_record.normalization_warnings = []
        variant_record.normalization_warnings.append(warning_msg)
        return variant_record

    # --- Core validation: ref AA vs sequence ----------------------------------
    if _validate_residue_match(sequence, position, ref_aa):
        # Match confirmed — direct mapping is valid
        variant_record.coordinate_mapping_confidence = "high"
        variant_record.coordinate_mapping_method = "direct"
        variant_record.uniprot_residue_position = position
        logger.info(
            "Variant %s: direct mapping confirmed — %s at position %d matches "
            "UniProt sequence",
            variant_record.variant_id,
            ref_aa,
            position,
        )
    else:
        # Mismatch — the position does not correspond to the expected residue
        actual_aa = sequence[position - 1]
        warning_msg = (
            f"Reference amino acid mismatch at position {position}: "
            f"expected '{ref_aa}' from HGVS notation but found '{actual_aa}' "
            f"in UniProt sequence (variant {variant_record.variant_id}). "
            f"This may indicate an isoform mismatch or incorrect transcript."
        )
        logger.warning(warning_msg)
        variant_record.coordinate_mapping_confidence = "failed"
        if variant_record.normalization_warnings is None:
            variant_record.normalization_warnings = []
        variant_record.normalization_warnings.append(warning_msg)

    return variant_record


def _resolve_canonical_transcript(gene_symbol: str) -> Optional[tuple[str, str]]:
    """Look up the MANE Select transcript for a gene.

    MANE (Matched Annotation from NCBI and EBI) Select provides a single
    canonical transcript agreed on by both RefSeq and Ensembl for each gene.
    This eliminates transcript ambiguity.

    Phase 1 stub: Returns None. Future phases will query the MANE Select
    database (available from NCBI FTP) to resolve the canonical transcript
    and corresponding protein accession for a given gene symbol.

    Args:
        gene_symbol: HGNC gene symbol, e.g. "BRCA1", "TP53".

    Returns:
        Tuple of (transcript_id, protein_id) if found, e.g.
        ("NM_007294.4", "NP_009225.1"), or None if the gene is not in
        the MANE Select set.
    """
    return None


def _map_to_uniprot_position(
    residue_position: int,
    transcript_id: str,
    uniprot_id: str,
) -> tuple[int, str, str]:
    """Map a transcript-based residue position to UniProt canonical coordinates.

    Strategy (future phases):
      1. Try direct mapping (positions often match for canonical isoform)
      2. If mismatch, use SIFTS database (PDBe structure integration)
      3. If SIFTS unavailable, use pairwise BLAST alignment

    Phase 1 stub: Returns the input position unchanged with "direct" method
    and "high" confidence. This assumes the HGVS-parsed position corresponds
    directly to the UniProt canonical isoform, which is validated downstream
    by checking the reference amino acid against the sequence.

    Args:
        residue_position: Residue position from the parsed HGVS notation.
        transcript_id: RefSeq transcript accession, e.g. "NM_007294.4".
        uniprot_id: UniProt accession, e.g. "P38398".

    Returns:
        Tuple of (uniprot_position, method, confidence) where:
          - uniprot_position: mapped position in UniProt coordinates
          - method: "direct", "sifts", or "blast_alignment"
          - confidence: "exact", "high", "low", or "failed"
    """
    return (residue_position, "direct", "high")


def _map_to_structure_position(
    uniprot_position: int,
    uniprot_id: str,
    pdb_path: Optional[str] = None,
) -> tuple[int, str, str]:
    """Map UniProt position to the 3D structure residue number.

    AlphaFold structures are usually numbered by UniProt position (direct),
    but experimental PDB structures may have different numbering, missing
    residues, or insertion codes.

    Phase 1 stub: Returns the input position unchanged with "direct" method
    and "high" confidence. This is correct for AlphaFold structures, which
    use UniProt numbering. Future phases will handle experimental PDB structures
    via SIFTS residue-level mapping.

    Args:
        uniprot_position: Residue position in UniProt canonical isoform
            coordinates.
        uniprot_id: UniProt accession, e.g. "P38398".
        pdb_path: Optional path to a PDB/mmCIF file. If provided, future
            phases will extract residue numbering from the structure itself.

    Returns:
        Tuple of (structure_position, method, confidence) where:
          - structure_position: residue number in the 3D structure
          - method: "direct", "sifts", or "blast_alignment"
          - confidence: "exact", "high", "low", or "failed"
    """
    return (uniprot_position, "direct", "high")


def _validate_residue_match(
    protein_sequence: str,
    position: int,
    expected_aa: str,
) -> bool:
    """Verify that the amino acid at the mapped position matches expectations.

    This is the critical safety check. If the reference amino acid at the
    mapped position doesn't match what we expect from the HGVS notation,
    the mapping is wrong and all downstream analysis would be on the wrong
    residue.

    Converts from 1-indexed biological position to 0-indexed Python string
    index (position 1 = sequence[0]).

    Args:
        protein_sequence: Full protein amino acid sequence (single-letter
            codes), e.g. "MDLSALREVE...".
        position: 1-indexed residue position to check.
        expected_aa: Single-letter amino acid code expected at this position,
            e.g. "R" for Arginine.

    Returns:
        True if the amino acid at sequence[position - 1] matches expected_aa,
        False if there is a mismatch.
    """
    return protein_sequence[position - 1] == expected_aa
