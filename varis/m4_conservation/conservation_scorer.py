"""Conservation Scorer — Calculates per-position conservation from alignment.

Uses Shannon entropy and percent identity at the mutation position.
Priority: 3 (easy — custom Python, no external tool)
Populates: conservation_score, position_entropy, conservation_method,
           msa_column_index, msa_gap_fraction_at_site, msa_num_sequences,
           num_orthologs, conserved_across_mammals, conservation_available
"""

import logging
import math
from typing import Optional

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# =============================================================================
# CONSTANTS
# =============================================================================

_STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")
_MAX_ENTROPY = math.log2(20)  # ≈ 4.322 bits

# Common mammalian NCBI taxon IDs
_MAMMAL_TAXON_IDS = {
    9606,   # Homo sapiens (human)
    9615,   # Canis lupus familiaris (dog)
    10090,  # Mus musculus (mouse)
    9913,   # Bos taurus (cow)
    9823,   # Sus scrofa (pig)
    10116,  # Rattus norvegicus (rat)
    9986,   # Oryctolagus cuniculus (rabbit)
    9685,   # Felis catus (cat)
    9796,   # Equus caballus (horse)
    13616,  # Monodelphis domestica (opossum)
    9544,   # Macaca mulatta (rhesus macaque)
    9598,   # Pan troglodytes (chimpanzee)
    9601,   # Pongo abelii (orangutan)
    9593,   # Gorilla gorilla (gorilla)
    30611,  # Otolemur garnettii (bushbaby)
    9669,   # Mustela putorius furo (ferret)
    9785,   # Loxodonta africana (elephant)
    9739,   # Tursiops truncatus (dolphin)
    9813,   # Procavia capensis (rock hyrax)
    132908, # Pteropus vampyrus (flying fox)
    59463,  # Myotis lucifugus (little brown bat)
    43179,  # Ictidomys tridecemlineatus (ground squirrel)
    10141,  # Cavia porcellus (guinea pig)
    9361,   # Dasypus novemcinctus (armadillo)
    9371,   # Echinops telfairi (tenrec)
    9478,   # Tarsius syrichta (tarsier)
    9555,   # Papio anubis (olive baboon)
    9541,   # Macaca fascicularis (crab-eating macaque)
    37347,  # Tupaia belangeri (tree shrew)
}

_MIN_MAMMALS = 5
_MAMMAL_CONSERVATION_THRESHOLD = 0.9


# =============================================================================
# PUBLIC API
# =============================================================================

def score_conservation(
    variant_record: VariantRecord,
    alignment: dict,
) -> VariantRecord:
    """Calculate conservation score at the mutation position from MSA.

    Reads the query sequence from the alignment, maps the variant's residue
    position to the alignment column (accounting for gaps), computes Shannon
    entropy at that column, and derives a conservation score normalized to
    [0, 1] where 1.0 = fully conserved and 0.0 = maximally variable.

    Also checks mammalian conservation if taxonomy data is available.

    Args:
        variant_record: The shared VariantRecord (needs residue_position,
            ref_aa_single populated by M1).
        alignment: Dict with keys "sequences" (id->aligned_seq),
            "query_id" (str), and optionally "taxonomy" (id->taxon_id).

    Returns:
        The variant_record with conservation fields populated, or set to
        None with reason codes on failure.
    """
    try:
        sequences = alignment.get("sequences", {})
        query_id = alignment.get("query_id")
        taxonomy = alignment.get("taxonomy", {})

        if not sequences or query_id is None or query_id not in sequences:
            logger.warning(
                "Alignment missing sequences or query_id for %s",
                variant_record.variant_id,
            )
            _set_all_null(variant_record, NullReason.NO_DATA_AVAILABLE)
            return variant_record

        query_row = sequences[query_id]
        position = variant_record.residue_position
        ref_aa = variant_record.ref_aa_single

        if position is None:
            logger.warning(
                "residue_position is None for %s — cannot score conservation",
                variant_record.variant_id,
            )
            _set_all_null(variant_record, NullReason.UPSTREAM_DEPENDENCY_FAILED)
            return variant_record

        # Map position to alignment column (gap-aware)
        col_idx = _map_position_to_column(query_row, position, expected_aa=ref_aa)
        if col_idx is None:
            logger.warning(
                "Position %d (ref=%s) could not be mapped to alignment column "
                "for %s — possible ref AA mismatch",
                position, ref_aa, variant_record.variant_id,
            )
            _set_all_null(variant_record, NullReason.COORDINATE_MAPPING_FAILED)
            return variant_record

        # Extract the column across all sequences
        column = _extract_column(sequences, col_idx)
        num_sequences = len(sequences)
        num_orthologs = num_sequences - 1  # exclude query

        # Compute entropy and conservation score
        entropy = _shannon_entropy(column)
        conservation_score = 1.0 - (entropy / _MAX_ENTROPY) if _MAX_ENTROPY > 0 else 1.0
        # Clamp to [0, 1] for numerical safety
        conservation_score = max(0.0, min(1.0, conservation_score))

        # Gap fraction at this site
        gap_count = sum(1 for aa in column if aa not in _STANDARD_AAS)
        gap_fraction = gap_count / len(column) if column else 0.0

        # Mammal conservation
        conserved_mammals = _mammal_conservation(
            sequences, taxonomy, col_idx, ref_aa,
        )

        # Populate the variant record
        variant_record.conservation_score = conservation_score
        variant_record.position_entropy = entropy
        variant_record.conservation_method = "clustal_omega"
        variant_record.msa_column_index = col_idx
        variant_record.msa_gap_fraction_at_site = gap_fraction
        variant_record.msa_num_sequences = num_sequences
        variant_record.num_orthologs = num_orthologs
        variant_record.conserved_across_mammals = conserved_mammals
        variant_record.set_feature_status("conservation", True)

        logger.info(
            "Conservation scored for %s: score=%.3f, entropy=%.3f, "
            "orthologs=%d, mammal_conserved=%s",
            variant_record.variant_id,
            conservation_score,
            entropy,
            num_orthologs,
            conserved_mammals,
        )

    except Exception as e:
        logger.warning(
            "Conservation scoring failed for %s: %s",
            variant_record.variant_id if variant_record else "unknown",
            e,
        )
        _set_all_null(variant_record, NullReason.TOOL_CRASHED)

    return variant_record


# =============================================================================
# PRIVATE HELPERS
# =============================================================================

def _map_position_to_column(
    query_row: str,
    position: int,
    expected_aa: Optional[str] = None,
) -> Optional[int]:
    """Map a 1-indexed residue position to a 0-indexed alignment column.

    Walks the aligned query sequence, counting only non-gap characters.
    When the count reaches ``position``, the current column index is returned.
    Gaps ('-', '.') in the query are skipped when counting residues.

    Args:
        query_row: The aligned query sequence (may contain gap characters).
        position: 1-indexed residue position in the unaligned sequence.
        expected_aa: If provided, verify the amino acid at that position
            matches. Returns None on mismatch.

    Returns:
        0-indexed column index in the alignment, or None if the position
        is out of range or the reference amino acid does not match.
    """
    residue_count = 0
    for col_idx, char in enumerate(query_row):
        if char not in ("-", "."):
            residue_count += 1
            if residue_count == position:
                if expected_aa is not None and char.upper() != expected_aa.upper():
                    return None
                return col_idx
    # Position beyond sequence length
    return None


def _shannon_entropy(column: list[str]) -> float:
    """Calculate Shannon entropy (in bits) for an alignment column.

    Only standard amino acids (the 20 canonical residues) are counted.
    Gaps ('-', '.') and non-standard characters are excluded from
    frequency calculation.

    Args:
        column: List of single-character amino acids (may include gaps).

    Returns:
        Shannon entropy in bits. 0.0 for a fully conserved column,
        log2(20) ≈ 4.322 for a maximally variable column.
    """
    # Filter to standard AAs only
    filtered = [aa.upper() for aa in column if aa.upper() in _STANDARD_AAS]
    if not filtered:
        return 0.0

    total = len(filtered)
    # Count frequencies
    freq: dict[str, int] = {}
    for aa in filtered:
        freq[aa] = freq.get(aa, 0) + 1

    # H = -Σ p_i * log2(p_i)
    entropy = 0.0
    for count in freq.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)

    return entropy


def _extract_column(sequences: dict[str, str], col_idx: int) -> list[str]:
    """Extract a single column from the alignment.

    Args:
        sequences: Dict of sequence_id -> aligned_sequence.
        col_idx: 0-indexed column index.

    Returns:
        List of characters at that column, one per sequence.
    """
    column = []
    for seq in sequences.values():
        if col_idx < len(seq):
            column.append(seq[col_idx])
    return column


def _mammal_conservation(
    sequences: dict[str, str],
    taxonomy: dict[str, int],
    col_idx: int,
    ref_aa: Optional[str],
) -> Optional[bool]:
    """Check whether the reference amino acid is conserved across mammals.

    Requires at least ``_MIN_MAMMALS`` (5) mammalian sequences in the
    alignment. If met, returns True when >= 90% of mammals have the
    reference amino acid at the given column.

    Args:
        sequences: Dict of sequence_id -> aligned_sequence.
        taxonomy: Dict of sequence_id -> NCBI taxon ID.
        col_idx: 0-indexed alignment column to check.
        ref_aa: The reference amino acid (single letter) to look for.

    Returns:
        True if >= 90% of mammals share the ref AA at this column,
        False if < 90%, or None if fewer than _MIN_MAMMALS mammals
        are available.
    """
    if ref_aa is None:
        return None

    ref_aa_upper = ref_aa.upper()
    mammal_aas: list[str] = []

    for seq_id, taxon_id in taxonomy.items():
        if taxon_id in _MAMMAL_TAXON_IDS and seq_id in sequences:
            seq = sequences[seq_id]
            if col_idx < len(seq):
                mammal_aas.append(seq[col_idx].upper())

    if len(mammal_aas) < _MIN_MAMMALS:
        return None

    matching = sum(1 for aa in mammal_aas if aa == ref_aa_upper)
    fraction = matching / len(mammal_aas)
    return fraction >= _MAMMAL_CONSERVATION_THRESHOLD


def _set_all_null(variant_record: VariantRecord, reason: str) -> None:
    """Set all conservation fields to None with the given reason code.

    Args:
        variant_record: The VariantRecord to update.
        reason: A NullReason constant explaining why the fields are null.
    """
    conservation_fields = [
        "conservation_score",
        "position_entropy",
        "conservation_method",
        "msa_column_index",
        "msa_gap_fraction_at_site",
        "msa_num_sequences",
        "num_orthologs",
        "conserved_across_mammals",
    ]
    for field_name in conservation_fields:
        variant_record.set_with_reason(field_name, None, reason)
    variant_record.set_feature_status("conservation", False, reason)
