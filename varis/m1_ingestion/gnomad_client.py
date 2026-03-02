"""gnomAD Client -- Retrieves population allele frequencies from gnomAD (Broad Institute).

Queries the gnomAD GraphQL API using genomic coordinates provided by ClinVar.
Only queries when ClinVar has provided genomic coordinates AND the reference
build matches GRCh38 (the build used by gnomAD v4).

Populates: gnomad_frequency, gnomad_popmax, gnomad_homozygotes.
"""

from __future__ import annotations

import logging
from typing import Optional

import httpx

from varis.models.variant_record import VariantRecord, NullReason
from varis.config import GNOMAD_API_URL

logger = logging.getLogger(__name__)

_TIMEOUT = 15.0
_GNOMAD_BUILD = "GRCh38"

_GNOMAD_FIELDS = ("gnomad_frequency", "gnomad_popmax", "gnomad_homozygotes")

_GNOMAD_QUERY = """
query GnomadVariant($variantId: String!) {
  variant(variantId: $variantId, dataset: gnomad_r4) {
    genome {
      ac
      an
      homozygote_count
      populations {
        id
        ac
        an
      }
    }
    exome {
      ac
      an
      homozygote_count
      populations {
        id
        ac
        an
      }
    }
  }
}
"""


def fetch_gnomad(
    variant_record: VariantRecord,
    client: Optional[httpx.Client] = None,
) -> VariantRecord:
    """Query gnomAD for allele frequency across populations.

    Only queries gnomAD when ClinVar has provided genomic coordinates
    and the reference build matches GRCh38.

    Args:
        variant_record: Must have clinvar_chrom, clinvar_pos, clinvar_ref,
            and clinvar_alt set (populated by ClinVar client).
        client: Optional httpx.Client for dependency injection (e.g., testing).
            If None, a new client is created internally.

    Returns:
        VariantRecord with gnomAD frequency fields populated, or None fields
        with reasons on failure.
    """
    # Gate 1: Check that genomic coordinates are available
    coords = (
        variant_record.clinvar_chrom,
        variant_record.clinvar_pos,
        variant_record.clinvar_ref,
        variant_record.clinvar_alt,
    )
    if not all(c is not None for c in coords):
        logger.info(
            "Skipping gnomAD lookup for %s: no genomic coordinates available",
            variant_record.variant_id,
        )
        for field_name in _GNOMAD_FIELDS:
            variant_record.set_with_reason(
                field_name, None, "no_genomic_coordinates"
            )
        return variant_record

    # Gate 2: Check that reference build matches gnomAD endpoint
    if variant_record.reference_build != _GNOMAD_BUILD:
        logger.info(
            "Skipping gnomAD lookup for %s: reference build %s != %s",
            variant_record.variant_id,
            variant_record.reference_build,
            _GNOMAD_BUILD,
        )
        for field_name in _GNOMAD_FIELDS:
            variant_record.set_with_reason(
                field_name, None, NullReason.VALIDATION_FAILED
            )
        return variant_record

    try:
        result = _query_gnomad_graphql(
            chrom=variant_record.clinvar_chrom,
            pos=variant_record.clinvar_pos,
            ref=variant_record.clinvar_ref,
            alt=variant_record.clinvar_alt,
            client=client,
        )

        if result is None:
            logger.info(
                "gnomAD returned no data for %s",
                variant_record.variant_id,
            )
            for field_name in _GNOMAD_FIELDS:
                variant_record.set_with_reason(
                    field_name, None, NullReason.NO_DATA_AVAILABLE
                )
            return variant_record

        # Populate fields from gnomAD response
        variant_record.set_with_reason("gnomad_frequency", result["frequency"])
        variant_record.set_with_reason("gnomad_popmax", result["popmax"])
        variant_record.set_with_reason("gnomad_homozygotes", result["homozygotes"])

        logger.info(
            "gnomAD data retrieved for %s: frequency=%.2e, popmax=%.2e, "
            "homozygotes=%d",
            variant_record.variant_id,
            result["frequency"],
            result["popmax"],
            result["homozygotes"],
        )

    except Exception as e:
        logger.warning(
            "gnomAD lookup failed for %s: %s",
            variant_record.variant_id,
            e,
        )
        for field_name in _GNOMAD_FIELDS:
            variant_record.set_with_reason(
                field_name, None, NullReason.TOOL_CRASHED
            )

    return variant_record


def _query_gnomad_graphql(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    client: Optional[httpx.Client] = None,
) -> Optional[dict]:
    """Execute GraphQL query against gnomAD API.

    Constructs a variant ID in gnomAD format (chrom-pos-ref-alt) and queries
    both genome and exome datasets, combining allele counts to compute
    overall frequency and population-maximum frequency.

    Args:
        chrom: Chromosome (e.g., "17").
        pos: Genomic position.
        ref: Reference allele (e.g., "G").
        alt: Alternate allele (e.g., "A").
        client: Optional httpx.Client for dependency injection.

    Returns:
        Dict with keys "frequency", "popmax", "homozygotes", or None if the
        variant is not found in gnomAD.
    """
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"
    variables = {"variantId": variant_id}

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        response = client.post(
            GNOMAD_API_URL,
            json={"query": _GNOMAD_QUERY, "variables": variables},
        )
        response.raise_for_status()
        data = response.json()
    finally:
        if should_close:
            client.close()

    # Navigate the GraphQL response
    variant_data = data.get("data", {}).get("variant")
    if variant_data is None:
        return None

    genome = variant_data.get("genome")
    exome = variant_data.get("exome")

    # If neither genome nor exome data exists, variant is not in gnomAD
    if genome is None and exome is None:
        return None

    # Aggregate allele counts across genome and exome
    total_ac = 0
    total_an = 0
    total_homozygotes = 0

    # Population-level counts: accumulate ac and an per population across
    # genome and exome datasets for popmax calculation
    pop_counts: dict[str, dict[str, int]] = {}

    for dataset in (genome, exome):
        if dataset is None:
            continue

        total_ac += dataset.get("ac", 0)
        total_an += dataset.get("an", 0)
        total_homozygotes += dataset.get("homozygote_count", 0)

        for pop in dataset.get("populations", []):
            pop_id = pop.get("id", "")
            if not pop_id:
                continue
            if pop_id not in pop_counts:
                pop_counts[pop_id] = {"ac": 0, "an": 0}
            pop_counts[pop_id]["ac"] += pop.get("ac", 0)
            pop_counts[pop_id]["an"] += pop.get("an", 0)

    # Calculate overall frequency
    frequency = total_ac / total_an if total_an > 0 else 0.0

    # Calculate popmax (maximum allele frequency across populations)
    popmax = 0.0
    for counts in pop_counts.values():
        if counts["an"] > 0:
            pop_freq = counts["ac"] / counts["an"]
            if pop_freq > popmax:
                popmax = pop_freq

    return {
        "frequency": frequency,
        "popmax": popmax,
        "homozygotes": total_homozygotes,
    }
