"""ClinVar Client — Retrieves variant classification and clinical data from ClinVar (NIH).

Uses NCBI E-utilities API via httpx to search and fetch variant records.
Also extracts genomic coordinates from ClinVar for gnomAD lookups.

Populates: clinvar_id, clinvar_classification, clinvar_review_status, clinvar_conditions.
Populates genomic coords: reference_build, clinvar_chrom, clinvar_pos, clinvar_ref,
    clinvar_alt, coordinate_source.
"""

import logging
import time
import xml.etree.ElementTree as ET

import httpx

from varis.models.variant_record import VariantRecord, NullReason
from varis.config import CLINVAR_BASE_URL, NCBI_API_KEY

logger = logging.getLogger(__name__)

_TIMEOUT = 15.0

# Known ClinVar classification values for matching in XML parsing
_KNOWN_CLASSIFICATIONS = {
    "Pathogenic",
    "Likely pathogenic",
    "Uncertain significance",
    "Likely benign",
    "Benign",
    "Pathogenic/Likely pathogenic",
    "Benign/Likely benign",
    "Conflicting classifications of pathogenicity",
    "Conflicting interpretations of pathogenicity",
    "not provided",
    "drug response",
    "risk factor",
    "association",
    "protective",
    "Affects",
    "other",
}

# Fields that constitute the ClinVar identifier group
_CLINVAR_ID_FIELDS = (
    "clinvar_id",
    "clinvar_classification",
    "clinvar_review_status",
    "clinvar_conditions",
)

# Fields that constitute the genomic coordinate group
_COORD_FIELDS = (
    "reference_build",
    "clinvar_chrom",
    "clinvar_pos",
    "clinvar_ref",
    "clinvar_alt",
    "coordinate_source",
)


def fetch_clinvar(
    variant_record: VariantRecord, client: httpx.Client | None = None
) -> VariantRecord:
    """Query ClinVar for variant classification and clinical significance.

    Also extracts genomic coordinates (chrom/pos/ref/alt) from the ClinVar
    record for downstream gnomAD population frequency lookups.

    Args:
        variant_record: Must have gene_symbol and hgvs_protein set.
        client: Optional httpx.Client for dependency injection in tests.

    Returns:
        VariantRecord with ClinVar fields populated, or None fields with
        reason codes on failure.
    """
    gene = variant_record.gene_symbol
    hgvs = variant_record.hgvs_protein

    if not gene or not hgvs:
        logger.warning("ClinVar lookup requires gene_symbol and hgvs_protein")
        for f in _CLINVAR_ID_FIELDS + _COORD_FIELDS:
            variant_record.set_with_reason(
                f, None, NullReason.UPSTREAM_DEPENDENCY_FAILED
            )
        return variant_record

    try:
        variation_id = _search_clinvar(gene, hgvs, client=client)
        if not variation_id:
            logger.info(f"No ClinVar entry found for {gene} {hgvs}")
            for f in _CLINVAR_ID_FIELDS + _COORD_FIELDS:
                variant_record.set_with_reason(
                    f, None, NullReason.NO_DATA_AVAILABLE
                )
            return variant_record

        record_data = _fetch_clinvar_record(variation_id, client=client)
        if not record_data:
            logger.warning(
                f"ClinVar record fetch failed for variation {variation_id}"
            )
            for f in _CLINVAR_ID_FIELDS + _COORD_FIELDS:
                variant_record.set_with_reason(
                    f, None, NullReason.NO_DATA_AVAILABLE
                )
            return variant_record

        # Populate ClinVar classification fields
        variant_record.clinvar_id = record_data.get("variation_id")
        variant_record.clinvar_classification = record_data.get("classification")
        variant_record.clinvar_review_status = record_data.get("review_status")
        variant_record.clinvar_conditions = record_data.get("conditions")

        # Populate genomic coordinates if available
        assembly = record_data.get("assembly")
        chrom = record_data.get("chrom")
        pos = record_data.get("pos")
        ref = record_data.get("ref")
        alt = record_data.get("alt")

        if chrom and pos is not None:
            variant_record.reference_build = assembly
            variant_record.clinvar_chrom = chrom
            variant_record.clinvar_pos = int(pos) if pos is not None else None
            variant_record.clinvar_ref = ref
            variant_record.clinvar_alt = alt
            variant_record.coordinate_source = "clinvar"
        else:
            logger.info(
                f"No genomic coordinates in ClinVar for {gene} {hgvs}"
            )
            for f in _COORD_FIELDS:
                variant_record.set_with_reason(
                    f, None, NullReason.NO_DATA_AVAILABLE
                )

        logger.info(
            f"ClinVar: {gene} {hgvs} → {variant_record.clinvar_classification} "
            f"(VCV{variant_record.clinvar_id})"
        )
        return variant_record

    except Exception as e:
        logger.warning(f"ClinVar lookup failed for {gene} {hgvs}: {e}")
        for f in _CLINVAR_ID_FIELDS + _COORD_FIELDS:
            variant_record.set_with_reason(f, None, NullReason.TOOL_CRASHED)
        return variant_record


def _search_clinvar(
    gene: str, variant: str, client: httpx.Client | None = None
) -> str | None:
    """Search ClinVar for a variant and return the ClinVar variation ID.

    Args:
        gene: Gene symbol, e.g., "BRCA1".
        variant: HGVS protein notation, e.g., "p.Arg1699Trp".
        client: Optional httpx.Client for dependency injection.

    Returns:
        ClinVar variation ID string, or None if not found.
    """
    url = f"{CLINVAR_BASE_URL}/esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": f"{gene}[gene] AND {variant}[variant name]",
        "retmode": "json",
        "retmax": "5",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        if client is not None:
            response = client.get(url, params=params, timeout=_TIMEOUT)
        else:
            response = httpx.get(url, params=params, timeout=_TIMEOUT)

        response.raise_for_status()
        data = response.json()

        id_list = data.get("esearchresult", {}).get("idlist", [])
        if id_list:
            logger.debug(
                f"ClinVar search found {len(id_list)} result(s) for "
                f"{gene} {variant}, using first: {id_list[0]}"
            )
            return id_list[0]

        logger.debug(f"ClinVar search returned no results for {gene} {variant}")
        return None

    except httpx.TimeoutException:
        logger.warning(f"ClinVar search timed out for {gene} {variant}")
        return None
    except httpx.HTTPStatusError as e:
        logger.warning(
            f"ClinVar search HTTP error for {gene} {variant}: "
            f"{e.response.status_code}"
        )
        return None
    except Exception as e:
        logger.warning(f"ClinVar search failed for {gene} {variant}: {e}")
        return None


def _fetch_clinvar_record(
    variation_id: str, client: httpx.Client | None = None
) -> dict | None:
    """Fetch full ClinVar record for a given variation ID.

    Includes rate limiting (0.35s sleep) to stay within NCBI's rate limits
    (3 requests/sec without API key, 10 requests/sec with key).

    Args:
        variation_id: ClinVar variation ID.
        client: Optional httpx.Client for dependency injection.

    Returns:
        Parsed ClinVar record as dict, or None on failure.
    """
    url = f"{CLINVAR_BASE_URL}/efetch.fcgi"
    params = {
        "db": "clinvar",
        "id": variation_id,
        "rettype": "vcv",
        "is_variationid": "true",
        "retmode": "xml",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    # Rate limiting: NCBI allows 3 req/s without key, 10 req/s with key
    time.sleep(0.35)

    try:
        if client is not None:
            response = client.get(url, params=params, timeout=_TIMEOUT)
        else:
            response = httpx.get(url, params=params, timeout=_TIMEOUT)

        response.raise_for_status()
        xml_text = response.text

        return _parse_clinvar_xml(xml_text, fallback_id=variation_id)

    except httpx.TimeoutException:
        logger.warning(
            f"ClinVar fetch timed out for variation {variation_id}"
        )
        return None
    except httpx.HTTPStatusError as e:
        logger.warning(
            f"ClinVar fetch HTTP error for variation {variation_id}: "
            f"{e.response.status_code}"
        )
        return None
    except Exception as e:
        logger.warning(
            f"ClinVar fetch failed for variation {variation_id}: {e}"
        )
        return None


def _parse_clinvar_xml(xml_text: str, fallback_id: str) -> dict | None:
    """Parse ClinVar VCV XML response to extract classification and coordinates.

    Handles both the newer (2023+) XML format with GermlineClassification
    and the older format with plain Description elements.

    Args:
        xml_text: Raw XML response from ClinVar efetch.
        fallback_id: Variation ID to use if not found in XML.

    Returns:
        Dict with keys: variation_id, classification, review_status,
        conditions, assembly, chrom, pos, ref, alt. Or None on parse failure.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        logger.warning(f"ClinVar XML parse error: {e}")
        return None

    result = {
        "variation_id": fallback_id,
        "classification": None,
        "review_status": None,
        "conditions": [],
        "assembly": None,
        "chrom": None,
        "pos": None,
        "ref": None,
        "alt": None,
    }

    # --- Extract variation ID from VariationArchive Accession ---
    variation_archive = root.find(".//VariationArchive")
    if variation_archive is not None:
        accession = variation_archive.get("Accession", "")
        if accession:
            # Strip the "VCV" prefix to get the numeric ID, or use as-is
            result["variation_id"] = accession

    # --- Extract classification ---
    # Newer format (2023+): GermlineClassification > Description
    germline_desc = root.find(
        ".//GermlineClassification/Description"
    )
    if germline_desc is not None and germline_desc.text:
        result["classification"] = germline_desc.text.strip()
    else:
        # Older format: ClassifiedRecord > Classifications > ...
        # or just look for Description elements with known classification values
        for desc_elem in root.iter("Description"):
            if desc_elem.text and desc_elem.text.strip() in _KNOWN_CLASSIFICATIONS:
                result["classification"] = desc_elem.text.strip()
                break

    # --- Extract review status ---
    # Newer format
    review_elem = root.find(".//GermlineClassification/ReviewStatus")
    if review_elem is not None and review_elem.text:
        result["review_status"] = review_elem.text.strip()
    else:
        # Older format: look for ReviewStatus anywhere
        for rs_elem in root.iter("ReviewStatus"):
            if rs_elem.text:
                result["review_status"] = rs_elem.text.strip()
                break

    # --- Extract conditions (trait names) ---
    conditions = []
    # Try TraitName elements first (newer format)
    for trait_name in root.iter("TraitName"):
        if trait_name.text and trait_name.text.strip():
            name = trait_name.text.strip()
            if name not in conditions:
                conditions.append(name)

    # If no TraitName, try Trait > Name > ElementValue
    if not conditions:
        for trait in root.iter("Trait"):
            for name_elem in trait.findall("Name"):
                ev = name_elem.find("ElementValue")
                if ev is not None and ev.text and ev.text.strip():
                    name = ev.text.strip()
                    if name not in conditions:
                        conditions.append(name)

    result["conditions"] = conditions if conditions else None

    # --- Extract genomic coordinates from SequenceLocation ---
    # Prefer GRCh38 over GRCh37
    grch38_loc = None
    grch37_loc = None

    for seq_loc in root.iter("SequenceLocation"):
        assembly = seq_loc.get("Assembly", "")
        # Prefer entries that have allele info (variant-level, not gene-level)
        has_alleles = (
            seq_loc.get("referenceAlleleVCF") is not None
            or seq_loc.get("referenceAllele") is not None
        )
        if assembly == "GRCh38" and (has_alleles or grch38_loc is None):
            if has_alleles or grch38_loc is None:
                grch38_loc = seq_loc
        elif assembly == "GRCh37" and (has_alleles or grch37_loc is None):
            if has_alleles or grch37_loc is None:
                grch37_loc = seq_loc

    chosen_loc = grch38_loc if grch38_loc is not None else grch37_loc

    if chosen_loc is not None:
        result["assembly"] = chosen_loc.get("Assembly")
        result["chrom"] = chosen_loc.get("Chr")

        # Position: try "start" then "Start" (different ClinVar XML versions)
        pos = chosen_loc.get("start") or chosen_loc.get("Start")
        if pos is not None:
            try:
                result["pos"] = int(pos)
            except (ValueError, TypeError):
                result["pos"] = None

        # Reference allele: try VCF-specific then general attribute names
        result["ref"] = (
            chosen_loc.get("referenceAlleleVCF")
            or chosen_loc.get("referenceAllele")
        )

        # Alternate allele: try VCF-specific then general attribute names
        result["alt"] = (
            chosen_loc.get("alternateAlleleVCF")
            or chosen_loc.get("alternateAllele")
        )

    return result
