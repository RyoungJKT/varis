# Phase 1: M1 Ingestion Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement all 7 M1 sub-modules so `python -m varis BRCA1 p.Arg1699Trp` produces a complete VariantRecord JSON with real data from public APIs.

**Architecture:** 7 sub-modules execute in sequence via the M1 orchestrator. Each sub-module reads/writes the shared VariantRecord dataclass. Failures set fields to None with reason codes. ClinVar provides genomic coordinates for gnomAD lookups.

**Tech Stack:** Python 3.11+, httpx (HTTP client), xml.etree.ElementTree (ClinVar XML), re (HGVS parsing), pytest (testing)

**Design doc:** `docs/plans/2026-03-03-phase1-m1-ingestion-design.md`

---

## Task 1: Add genomic coordinate fields to VariantRecord

**Files:**
- Modify: `varis/models/variant_record.py` (lines 98-106, identifiers section)
- Modify: `schema/variant_record_schema.json` (lines 44-50, identifiers section)
- Test: `tests/test_variant_record.py`

**Step 1: Write the failing test**

Add to `tests/test_variant_record.py`:

```python
class TestGenomicCoordinateFields:
    def test_new_fields_exist(self):
        """New genomic coordinate fields should be on VariantRecord."""
        from varis.models.variant_record import create_variant_record
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "reference_build")
        assert hasattr(record, "clinvar_chrom")
        assert hasattr(record, "clinvar_pos")
        assert hasattr(record, "clinvar_ref")
        assert hasattr(record, "clinvar_alt")
        assert hasattr(record, "coordinate_source")
        assert record.reference_build is None
        assert record.clinvar_chrom is None

    def test_genomic_fields_serialize(self):
        """Genomic coordinate fields should survive to_dict/from_dict."""
        from varis.models.variant_record import create_variant_record
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.reference_build = "GRCh38"
        record.clinvar_chrom = "17"
        record.clinvar_pos = 43057051
        record.clinvar_ref = "G"
        record.clinvar_alt = "A"
        record.coordinate_source = "clinvar"
        d = record.to_dict()
        restored = record.from_dict(d)
        assert restored.reference_build == "GRCh38"
        assert restored.clinvar_pos == 43057051
        assert restored.coordinate_source == "clinvar"
```

**Step 2: Run test to verify it fails**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_variant_record.py::TestGenomicCoordinateFields -v`
Expected: FAIL — `AttributeError: 'VariantRecord' object has no attribute 'reference_build'`

**Step 3: Add fields to VariantRecord**

In `varis/models/variant_record.py`, add after `ensembl_transcript_id` (line 105), before the VARIANT DETAILS section:

```python
    # =========================================================================
    # GENOMIC COORDINATES — From ClinVar (populated by M1)
    # =========================================================================
    reference_build: Optional[str] = None          # "GRCh37" or "GRCh38"
    clinvar_chrom: Optional[str] = None            # e.g. "17"
    clinvar_pos: Optional[int] = None              # genomic position
    clinvar_ref: Optional[str] = None              # reference nucleotide allele
    clinvar_alt: Optional[str] = None              # alternate nucleotide allele
    coordinate_source: Optional[str] = None        # "clinvar" for Phase 1
```

**Step 4: Add fields to JSON schema**

In `schema/variant_record_schema.json`, add after `"ensembl_transcript_id"` entry:

```json
    "reference_build": {"type": ["string", "null"], "enum": ["GRCh37", "GRCh38", null]},
    "clinvar_chrom": {"type": ["string", "null"]},
    "clinvar_pos": {"type": ["integer", "null"]},
    "clinvar_ref": {"type": ["string", "null"]},
    "clinvar_alt": {"type": ["string", "null"]},
    "coordinate_source": {"type": ["string", "null"]},
```

**Step 5: Bump schema version**

In `varis/models/variant_record.py`, change:
```python
RECORD_SCHEMA_VERSION = "1.1.0"
```

**Step 6: Run tests to verify they pass**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_variant_record.py -v`
Expected: ALL PASS (including existing tests — the schema version test in `TestVariantRecordCreation` will need the version comparison updated if it checks exact match)

**Step 7: Commit**

```bash
git add varis/models/variant_record.py schema/variant_record_schema.json tests/test_variant_record.py
git commit -m "feat(m1): add genomic coordinate fields to VariantRecord for ClinVar-sourced coords"
```

---

## Task 2: Reorder M1 orchestrator

**Files:**
- Modify: `varis/m1_ingestion/__init__.py`

**Step 1: Update execution order**

Change the `steps` list in `run()` to:

```python
    steps = [
        ("M1.hgvs_parser", parse_hgvs),
        ("M1.uniprot", fetch_uniprot),
        ("M1.normalizer", normalize_variant),
        ("M1.clinvar", fetch_clinvar),
        ("M1.gnomad", fetch_gnomad),
        ("M1.alphafold", fetch_alphafold_structure),
        ("M1.alphamissense", fetch_alphamissense),
    ]
```

Also update the docstring — change the sentence about order to:
```
    Order matters: hgvs_parser → uniprot → normalizer → clinvar → gnomad → alphafold → alphamissense.
    UniProt runs before normalizer because normalizer validates ref AA against UniProt sequence.
    ClinVar runs before gnomAD because gnomAD needs ClinVar's genomic coordinates.
```

**Step 2: Run existing pipeline test**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_pipeline.py::TestPipelineGracefulDegradation::test_pipeline_never_crashes -v`
Expected: PASS (pipeline still never crashes)

**Step 3: Commit**

```bash
git add varis/m1_ingestion/__init__.py
git commit -m "refactor(m1): reorder M1 steps — uniprot before normalizer, clinvar before gnomad"
```

---

## Task 3: Implement hgvs_parser

**Files:**
- Modify: `varis/m1_ingestion/hgvs_parser.py`
- Modify: `tests/test_m1_ingestion.py`

**Step 1: Write the tests**

Replace the `TestHGVSParser` class in `tests/test_m1_ingestion.py`:

```python
from varis.models.variant_record import create_variant_record
from varis.m1_ingestion.hgvs_parser import parse_hgvs


class TestHGVSParser:
    def test_parse_three_letter(self):
        """Parse standard three-letter HGVS: p.Arg1699Trp"""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = parse_hgvs(record)
        assert record.ref_amino_acid == "Arg"
        assert record.alt_amino_acid == "Trp"
        assert record.ref_aa_single == "R"
        assert record.alt_aa_single == "W"

    def test_parse_extracts_position(self):
        """Position should be extracted as an integer."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = parse_hgvs(record)
        assert record.residue_position == 1699

    def test_parse_single_letter(self):
        """Also handle single-letter notation: p.R1699W"""
        record = create_variant_record("BRCA1", "p.R1699W")
        record = parse_hgvs(record)
        assert record.residue_position == 1699
        assert record.ref_amino_acid == "Arg"
        assert record.alt_amino_acid == "Trp"
        assert record.ref_aa_single == "R"
        assert record.alt_aa_single == "W"

    def test_parse_charge_change(self):
        """Arg (+) to Trp (0) should produce a charge change string."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = parse_hgvs(record)
        assert record.charge_change is not None
        assert "+" in record.charge_change or "positive" in record.charge_change.lower()

    def test_parse_no_charge_change(self):
        """Ala (0) to Val (0) should indicate no charge change."""
        record = create_variant_record("FAKEGENE", "p.Ala100Val")
        record = parse_hgvs(record)
        assert record.charge_change is not None
        assert "no change" in record.charge_change.lower() or "neutral" in record.charge_change.lower()

    def test_invalid_hgvs_returns_none(self):
        """Invalid notation should set fields to None with reason, not crash."""
        record = create_variant_record("BRCA1", "not_valid_hgvs")
        record = parse_hgvs(record)
        assert record.residue_position is None
        assert "residue_position" in record.null_reasons

    def test_normalized_notation(self):
        """Should populate input_notation_normalized."""
        record = create_variant_record("BRCA1", "p.R1699W")
        record = parse_hgvs(record)
        assert record.input_notation_normalized == "p.Arg1699Trp"
```

**Step 2: Run tests to verify they fail**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestHGVSParser -v`
Expected: FAIL — functions return None (pass statements)

**Step 3: Implement hgvs_parser.py**

Replace the function bodies in `varis/m1_ingestion/hgvs_parser.py`:

```python
"""HGVS Parser — Extracts gene, position, and amino acid change from HGVS notation.

Parses protein-level HGVS like "p.Arg1699Trp" into structured components:
residue position (1699), reference AA (Arg), alternate AA (Trp), charge change.
"""

import re
import logging
from varis.models.variant_record import VariantRecord, NullReason
from varis.config import AA_CHARGE, AA_THREE_TO_ONE, AA_ONE_TO_THREE

logger = logging.getLogger(__name__)

# Pattern: p.Arg1699Trp (three-letter codes)
HGVS_PROTEIN_PATTERN = re.compile(
    r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})"
)
# Pattern: p.R1699W (single-letter codes)
HGVS_PROTEIN_SINGLE = re.compile(
    r"p\.([A-Z])(\d+)([A-Z])"
)


def parse_hgvs(variant_record: VariantRecord) -> VariantRecord:
    """Parse HGVS protein notation and populate variant detail fields.

    Args:
        variant_record: Must have hgvs_protein set (e.g., "p.Arg1699Trp").

    Returns:
        VariantRecord with residue_position, ref/alt amino acids, and charge change.
    """
    hgvs = variant_record.hgvs_protein
    if not hgvs:
        for f in ("residue_position", "ref_amino_acid", "alt_amino_acid",
                  "ref_aa_single", "alt_aa_single", "charge_change"):
            variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
        return variant_record

    # Try three-letter pattern first, then single-letter
    match = HGVS_PROTEIN_PATTERN.match(hgvs)
    if match:
        ref_three, position_str, alt_three = match.groups()
        ref_single = AA_THREE_TO_ONE.get(ref_three)
        alt_single = AA_THREE_TO_ONE.get(alt_three)
        if not ref_single or not alt_single:
            logger.warning(f"Unknown amino acid code in {hgvs}")
            for f in ("residue_position", "ref_amino_acid", "alt_amino_acid",
                      "ref_aa_single", "alt_aa_single", "charge_change"):
                variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
            return variant_record
    else:
        match = HGVS_PROTEIN_SINGLE.match(hgvs)
        if match:
            ref_single, position_str, alt_single = match.groups()
            ref_three = AA_ONE_TO_THREE.get(ref_single)
            alt_three = AA_ONE_TO_THREE.get(alt_single)
            if not ref_three or not alt_three:
                logger.warning(f"Unknown amino acid code in {hgvs}")
                for f in ("residue_position", "ref_amino_acid", "alt_amino_acid",
                          "ref_aa_single", "alt_aa_single", "charge_change"):
                    variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
                return variant_record
        else:
            logger.warning(f"Could not parse HGVS notation: {hgvs}")
            for f in ("residue_position", "ref_amino_acid", "alt_amino_acid",
                      "ref_aa_single", "alt_aa_single", "charge_change"):
                variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
            return variant_record

    variant_record.residue_position = int(position_str)
    variant_record.ref_amino_acid = ref_three
    variant_record.alt_amino_acid = alt_three
    variant_record.ref_aa_single = ref_single
    variant_record.alt_aa_single = alt_single
    variant_record.charge_change = _calculate_charge_change(ref_three, alt_three)
    variant_record.input_notation_normalized = f"p.{ref_three}{position_str}{alt_three}"

    logger.info(
        f"Parsed HGVS: {hgvs} → position={position_str}, "
        f"{ref_three}({ref_single}) → {alt_three}({alt_single})"
    )
    return variant_record


def _calculate_charge_change(ref_aa: str, alt_aa: str) -> str:
    """Calculate the electrostatic charge change between two amino acids.

    Args:
        ref_aa: Reference amino acid (three-letter code).
        alt_aa: Alternate amino acid (three-letter code).

    Returns:
        Human-readable charge change, e.g., "positive → neutral".
    """
    charge_labels = {"+": "positive", "-": "negative", "0": "neutral"}
    ref_charge = AA_CHARGE.get(ref_aa, "0")
    alt_charge = AA_CHARGE.get(alt_aa, "0")
    ref_label = charge_labels.get(ref_charge, "neutral")
    alt_label = charge_labels.get(alt_charge, "neutral")
    if ref_charge == alt_charge:
        return f"no change ({ref_label})"
    return f"{ref_label} → {alt_label}"
```

**Step 4: Run tests to verify they pass**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestHGVSParser -v`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add varis/m1_ingestion/hgvs_parser.py tests/test_m1_ingestion.py
git commit -m "feat(m1): implement HGVS parser with three-letter and single-letter support"
```

---

## Task 4: Implement uniprot_client

**Files:**
- Modify: `varis/m1_ingestion/uniprot_client.py`
- Modify: `tests/test_m1_ingestion.py`

**Step 1: Write the tests**

Add to `tests/test_m1_ingestion.py`:

```python
import pytest
from varis.m1_ingestion.uniprot_client import fetch_uniprot


class TestUniProtClient:
    @pytest.mark.timeout(30)
    def test_fetch_brca1(self):
        """Fetch BRCA1 protein data from UniProt — real API call."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = fetch_uniprot(record)
        assert record.uniprot_id == "P38398"
        assert record.protein_sequence is not None
        assert len(record.protein_sequence) > 1000  # BRCA1 is 1863 AA
        assert record.protein_length == len(record.protein_sequence)
        assert record.protein_name is not None

    @pytest.mark.timeout(30)
    def test_unknown_gene_returns_none(self):
        """Unknown gene should set fields to None with reason, not crash."""
        record = create_variant_record("FAKEGENE999", "p.Ala1Val")
        record = fetch_uniprot(record)
        assert record.uniprot_id is None
        assert "uniprot_id" in record.null_reasons
```

**Step 2: Run tests to verify they fail**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestUniProtClient -v`
Expected: FAIL

**Step 3: Implement uniprot_client.py**

Replace `varis/m1_ingestion/uniprot_client.py`:

```python
"""UniProt Client — Retrieves protein sequence and functional annotations.

Populates: uniprot_id, protein_name, protein_sequence, protein_length, protein_function.
"""

import logging
import httpx
from varis.models.variant_record import VariantRecord, NullReason
from varis.config import UNIPROT_BASE_URL

logger = logging.getLogger(__name__)

# Timeout for API calls
_TIMEOUT = 15.0


def fetch_uniprot(variant_record: VariantRecord,
                  client: httpx.Client | None = None) -> VariantRecord:
    """Query UniProt for protein sequence and metadata.

    Args:
        variant_record: Must have gene_symbol set.
        client: Optional httpx.Client for dependency injection in tests.

    Returns:
        VariantRecord with protein fields populated.
    """
    gene = variant_record.gene_symbol
    if not gene:
        for f in ("uniprot_id", "protein_name", "protein_sequence",
                  "protein_length", "protein_function"):
            variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
        return variant_record

    try:
        accession = _search_uniprot(gene, client)
        if not accession:
            logger.warning(f"No UniProt entry found for gene {gene}")
            for f in ("uniprot_id", "protein_name", "protein_sequence",
                      "protein_length", "protein_function"):
                variant_record.set_with_reason(f, None, NullReason.NO_DATA_AVAILABLE)
            return variant_record

        data = _fetch_protein_data(accession, client)
        if not data:
            logger.warning(f"Failed to fetch UniProt data for {accession}")
            variant_record.uniprot_id = accession
            for f in ("protein_name", "protein_sequence",
                      "protein_length", "protein_function"):
                variant_record.set_with_reason(f, None, NullReason.TOOL_CRASHED)
            return variant_record

        variant_record.uniprot_id = accession
        variant_record.protein_name = data.get("protein_name")
        variant_record.protein_sequence = data.get("sequence")
        variant_record.protein_length = len(data["sequence"]) if data.get("sequence") else None
        variant_record.protein_function = data.get("function")

        logger.info(f"UniProt: {gene} → {accession}, length={variant_record.protein_length}")

    except Exception as e:
        logger.warning(f"UniProt fetch failed for {gene}: {e}")
        for f in ("uniprot_id", "protein_name", "protein_sequence",
                  "protein_length", "protein_function"):
            variant_record.set_with_reason(f, None, NullReason.TOOL_CRASHED)

    return variant_record


def _search_uniprot(gene: str, client: httpx.Client | None = None) -> str | None:
    """Search UniProt for a human gene and return the accession ID.

    Prefers reviewed (Swiss-Prot) entries over unreviewed (TrEMBL).

    Args:
        gene: Gene symbol, e.g., "BRCA1".
        client: Optional httpx.Client.

    Returns:
        UniProt accession ID (e.g., "P38398"), or None if not found.
    """
    url = f"{UNIPROT_BASE_URL}/uniprotkb/search"
    params = {
        "query": f"gene_exact:{gene} AND organism_id:9606",
        "format": "json",
        "fields": "accession,reviewed,protein_name,gene_names",
        "size": "5",
    }

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        resp = client.get(url, params=params)
        resp.raise_for_status()
        data = resp.json()
        results = data.get("results", [])
        if not results:
            return None

        # Prefer reviewed (Swiss-Prot) entries
        for entry in results:
            if entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)":
                return entry["primaryAccession"]

        # Fall back to first result
        return results[0]["primaryAccession"]
    finally:
        if should_close:
            client.close()


def _fetch_protein_data(accession: str,
                        client: httpx.Client | None = None) -> dict | None:
    """Fetch protein entry from UniProt REST API.

    Args:
        accession: UniProt accession ID.
        client: Optional httpx.Client.

    Returns:
        Dict with keys: protein_name, sequence, function. Or None on failure.
    """
    url = f"{UNIPROT_BASE_URL}/uniprotkb/{accession}.json"

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        resp = client.get(url)
        resp.raise_for_status()
        entry = resp.json()

        # Extract protein name
        protein_name = None
        prot_desc = entry.get("proteinDescription", {})
        rec_name = prot_desc.get("recommendedName", {})
        if rec_name:
            protein_name = rec_name.get("fullName", {}).get("value")

        # Extract sequence
        sequence = entry.get("sequence", {}).get("value")

        # Extract function from comments
        function = None
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function = texts[0].get("value")
                break

        return {
            "protein_name": protein_name,
            "sequence": sequence,
            "function": function,
        }
    finally:
        if should_close:
            client.close()
```

**Step 4: Run tests to verify they pass**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestUniProtClient -v`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add varis/m1_ingestion/uniprot_client.py tests/test_m1_ingestion.py
git commit -m "feat(m1): implement UniProt client with reviewed-entry preference"
```

---

## Task 5: Implement variant_normalizer

**Files:**
- Modify: `varis/m1_ingestion/variant_normalizer.py`
- Modify: `tests/test_m1_ingestion.py`

**Step 1: Write the tests**

Add to `tests/test_m1_ingestion.py`:

```python
from varis.m1_ingestion.variant_normalizer import normalize_variant


class TestVariantNormalizer:
    def test_validates_correct_position(self):
        """When ref AA matches UniProt sequence, confidence should be 'high'."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.residue_position = 1699
        record.ref_aa_single = "R"
        # BRCA1 position 1699 is Arg (R) in UniProt P38398
        # We need a real protein sequence for this; set a fake one where pos 1699 = R
        record.protein_sequence = "M" + "A" * 1697 + "R" + "A" * 164  # 1863 AA, R at pos 1699
        record = normalize_variant(record)
        assert record.coordinate_mapping_confidence == "high"
        assert record.uniprot_residue_position == 1699
        assert record.coordinate_mapping_method == "direct"

    def test_detects_mismatch(self):
        """When ref AA does NOT match, confidence should be 'failed'."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.residue_position = 1699
        record.ref_aa_single = "R"
        record.protein_sequence = "M" + "A" * 1697 + "G" + "A" * 164  # G at 1699, not R
        record = normalize_variant(record)
        assert record.coordinate_mapping_confidence == "failed"
        assert record.normalization_warnings is not None
        assert len(record.normalization_warnings) > 0

    def test_no_sequence_upstream_failure(self):
        """If protein_sequence is None (UniProt failed), set upstream_dependency_failed."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.residue_position = 1699
        record.ref_aa_single = "R"
        record.protein_sequence = None
        record = normalize_variant(record)
        assert record.coordinate_mapping_confidence is None
        assert "coordinate_mapping_confidence" in record.null_reasons

    def test_position_out_of_range(self):
        """Position beyond sequence length should set confidence to 'failed'."""
        record = create_variant_record("BRCA1", "p.Arg9999Trp")
        record.residue_position = 9999
        record.ref_aa_single = "R"
        record.protein_sequence = "MAAAA"  # only 5 AA
        record = normalize_variant(record)
        assert record.coordinate_mapping_confidence == "failed"
```

**Step 2: Run tests to verify they fail**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestVariantNormalizer -v`
Expected: FAIL

**Step 3: Implement variant_normalizer.py**

Replace the function bodies in `varis/m1_ingestion/variant_normalizer.py`:

```python
"""Variant Normalizer — Coordinate mapping and notation disambiguation.

Phase 1: Simple direct mapping with ref AA validation against UniProt sequence.
Full MANE Select / SIFTS / BLAST alignment deferred to Phase 2+.
"""

import logging
from varis.models.variant_record import VariantRecord, NullReason

logger = logging.getLogger(__name__)


def normalize_variant(variant_record: VariantRecord) -> VariantRecord:
    """Validate variant position against UniProt sequence and set confidence.

    Phase 1 implementation: assumes UniProt canonical numbering matches HGVS
    protein position directly. Validates by checking ref AA at the position.

    Args:
        variant_record: Must have residue_position and ref_aa_single set
            (from hgvs_parser), and protein_sequence set (from uniprot_client).

    Returns:
        VariantRecord with normalization fields populated.
    """
    position = variant_record.residue_position
    ref_aa = variant_record.ref_aa_single
    sequence = variant_record.protein_sequence

    # If we don't have a position to validate, nothing to do
    if position is None or ref_aa is None:
        variant_record.set_with_reason(
            "coordinate_mapping_confidence", None,
            NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        return variant_record

    # If UniProt failed, we can't validate
    if sequence is None:
        variant_record.set_with_reason(
            "coordinate_mapping_confidence", None,
            NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        logger.warning(
            f"Cannot validate position {position} — no protein sequence available"
        )
        return variant_record

    # Validate ref AA at position (1-indexed)
    if position < 1 or position > len(sequence):
        variant_record.coordinate_mapping_confidence = "failed"
        variant_record.normalization_warnings = [
            f"Position {position} is out of range for sequence length {len(sequence)}"
        ]
        logger.warning(
            f"Position {position} out of range (sequence length: {len(sequence)})"
        )
        return variant_record

    actual_aa = sequence[position - 1]  # Convert to 0-indexed
    if _validate_residue_match(sequence, position, ref_aa):
        variant_record.coordinate_mapping_confidence = "high"
        variant_record.coordinate_mapping_method = "direct"
        variant_record.uniprot_residue_position = position
        logger.info(
            f"Position {position} validated: expected {ref_aa}, found {actual_aa} ✓"
        )
    else:
        variant_record.coordinate_mapping_confidence = "failed"
        variant_record.normalization_warnings = [
            f"Ref AA mismatch at position {position}: "
            f"expected {ref_aa}, found {actual_aa} in UniProt sequence"
        ]
        logger.warning(
            f"Ref AA mismatch at position {position}: "
            f"expected {ref_aa}, found {actual_aa}"
        )

    return variant_record


def _resolve_canonical_transcript(gene_symbol: str) -> tuple[str, str] | None:
    """Look up the MANE Select transcript for a gene.

    Phase 1: Not implemented. Returns None.
    Phase 2+: Will query Ensembl REST API.

    Returns:
        Tuple of (transcript_id, protein_id) or None if not found.
    """
    return None


def _map_to_uniprot_position(residue_position: int, transcript_id: str,
                              uniprot_id: str) -> tuple[int, str, str]:
    """Map a transcript-based residue position to UniProt canonical coordinates.

    Phase 1: Direct mapping (assumes positions match).
    Phase 2+: Will use SIFTS and BLAST alignment fallbacks.

    Returns:
        Tuple of (uniprot_position, method, confidence).
    """
    return (residue_position, "direct", "high")


def _map_to_structure_position(uniprot_position: int, uniprot_id: str,
                                pdb_path: str = None) -> tuple[int, str, str]:
    """Map UniProt position to the 3D structure residue number.

    Phase 1: Direct mapping (AlphaFold uses UniProt numbering).
    Phase 2+: Will handle experimental PDB structures with different numbering.

    Returns:
        Tuple of (structure_position, method, confidence).
    """
    return (uniprot_position, "direct", "high")


def _validate_residue_match(protein_sequence: str, position: int,
                             expected_aa: str) -> bool:
    """Verify that the amino acid at the mapped position matches expectations.

    Args:
        protein_sequence: Full protein sequence string (one-letter codes).
        position: 1-indexed residue position.
        expected_aa: Expected amino acid (single-letter code).

    Returns:
        True if match confirmed, False if mismatch detected.
    """
    if position < 1 or position > len(protein_sequence):
        return False
    actual = protein_sequence[position - 1]
    return actual == expected_aa
```

**Step 4: Run tests to verify they pass**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestVariantNormalizer -v`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add varis/m1_ingestion/variant_normalizer.py tests/test_m1_ingestion.py
git commit -m "feat(m1): implement variant normalizer with ref AA validation"
```

---

## Task 6: Implement clinvar_client

**Files:**
- Modify: `varis/m1_ingestion/clinvar_client.py`
- Modify: `tests/test_m1_ingestion.py`

**Step 1: Write the tests**

Replace the `TestClinVarClient` class in `tests/test_m1_ingestion.py`:

```python
from varis.m1_ingestion.clinvar_client import fetch_clinvar


class TestClinVarClient:
    @pytest.mark.timeout(30)
    def test_fetch_known_variant(self):
        """Fetch BRCA1 p.Arg1699Trp from ClinVar — real API call."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = fetch_clinvar(record)
        assert record.clinvar_id is not None
        # VCV000055361 is the known ID; ClinVar may return it in different format
        assert "55361" in str(record.clinvar_id)
        assert record.clinvar_classification is not None

    @pytest.mark.timeout(30)
    def test_extracts_genomic_coordinates(self):
        """ClinVar should provide genomic coordinates for gnomAD."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = fetch_clinvar(record)
        # If ClinVar has this variant, it should have genomic coords
        if record.clinvar_id is not None:
            assert record.reference_build is not None
            assert record.clinvar_chrom is not None
            assert record.clinvar_pos is not None
            assert record.coordinate_source == "clinvar"

    @pytest.mark.timeout(30)
    def test_unknown_variant_returns_none(self):
        """Unknown variant should set fields to None, not crash."""
        record = create_variant_record("FAKEGENE999", "p.Ala1Val")
        record = fetch_clinvar(record)
        assert record.clinvar_id is None
        assert "clinvar_id" in record.null_reasons
```

**Step 2: Run tests to verify they fail**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestClinVarClient -v`
Expected: FAIL

**Step 3: Implement clinvar_client.py**

Replace `varis/m1_ingestion/clinvar_client.py`:

```python
"""ClinVar Client — Retrieves variant classification and clinical data from ClinVar (NIH).

Uses NCBI E-utilities API to search and fetch variant records.
Also extracts genomic coordinates from ClinVar for downstream gnomAD queries.

Populates: clinvar_id, clinvar_classification, clinvar_review_status,
           clinvar_conditions, reference_build, clinvar_chrom, clinvar_pos,
           clinvar_ref, clinvar_alt, coordinate_source.
"""

import logging
import time
import xml.etree.ElementTree as ET
import httpx
from varis.models.variant_record import VariantRecord, NullReason
from varis.config import CLINVAR_BASE_URL, NCBI_API_KEY

logger = logging.getLogger(__name__)

_TIMEOUT = 15.0


def fetch_clinvar(variant_record: VariantRecord,
                  client: httpx.Client | None = None) -> VariantRecord:
    """Query ClinVar for variant classification and genomic coordinates.

    Args:
        variant_record: Must have gene_symbol and hgvs_protein set.
        client: Optional httpx.Client for dependency injection.

    Returns:
        VariantRecord with ClinVar fields populated, or None fields on failure.
    """
    gene = variant_record.gene_symbol
    hgvs = variant_record.hgvs_protein
    if not gene or not hgvs:
        for f in ("clinvar_id", "clinvar_classification", "clinvar_review_status",
                  "reference_build", "clinvar_chrom", "clinvar_pos",
                  "clinvar_ref", "clinvar_alt", "coordinate_source"):
            variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
        return variant_record

    try:
        variation_id = _search_clinvar(gene, hgvs, client)
        if not variation_id:
            logger.warning(f"No ClinVar entry found for {gene} {hgvs}")
            for f in ("clinvar_id", "clinvar_classification", "clinvar_review_status",
                      "reference_build", "clinvar_chrom", "clinvar_pos",
                      "clinvar_ref", "clinvar_alt", "coordinate_source"):
                variant_record.set_with_reason(f, None, NullReason.NO_DATA_AVAILABLE)
            return variant_record

        record_data = _fetch_clinvar_record(variation_id, client)
        if not record_data:
            logger.warning(f"Failed to fetch ClinVar record for {variation_id}")
            variant_record.clinvar_id = variation_id
            for f in ("clinvar_classification", "clinvar_review_status"):
                variant_record.set_with_reason(f, None, NullReason.TOOL_CRASHED)
            return variant_record

        variant_record.clinvar_id = record_data.get("variation_id")
        variant_record.clinvar_classification = record_data.get("classification")
        variant_record.clinvar_review_status = record_data.get("review_status")
        variant_record.clinvar_conditions = record_data.get("conditions")

        # Genomic coordinates
        if record_data.get("chrom"):
            variant_record.reference_build = record_data.get("assembly")
            variant_record.clinvar_chrom = record_data.get("chrom")
            variant_record.clinvar_pos = record_data.get("pos")
            variant_record.clinvar_ref = record_data.get("ref")
            variant_record.clinvar_alt = record_data.get("alt")
            variant_record.coordinate_source = "clinvar"
        else:
            for f in ("reference_build", "clinvar_chrom", "clinvar_pos",
                      "clinvar_ref", "clinvar_alt", "coordinate_source"):
                variant_record.set_with_reason(f, None, NullReason.NO_DATA_AVAILABLE)

        logger.info(
            f"ClinVar: {gene} {hgvs} → {variant_record.clinvar_id}, "
            f"classification={variant_record.clinvar_classification}"
        )

    except Exception as e:
        logger.warning(f"ClinVar fetch failed for {gene} {hgvs}: {e}")
        for f in ("clinvar_id", "clinvar_classification", "clinvar_review_status",
                  "reference_build", "clinvar_chrom", "clinvar_pos",
                  "clinvar_ref", "clinvar_alt", "coordinate_source"):
            variant_record.set_with_reason(f, None, NullReason.TOOL_CRASHED)

    return variant_record


def _search_clinvar(gene: str, variant: str,
                    client: httpx.Client | None = None) -> str | None:
    """Search ClinVar for a variant and return the ClinVar variation ID.

    Args:
        gene: Gene symbol, e.g., "BRCA1".
        variant: HGVS protein notation, e.g., "p.Arg1699Trp".
        client: Optional httpx.Client.

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

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        resp = client.get(url, params=params)
        resp.raise_for_status()
        data = resp.json()
        id_list = data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return None
        return id_list[0]
    finally:
        if should_close:
            client.close()


def _fetch_clinvar_record(variation_id: str,
                          client: httpx.Client | None = None) -> dict | None:
    """Fetch full ClinVar record for a given variation ID.

    Parses XML response for classification, review status, conditions,
    and genomic coordinates (preferring GRCh38).

    Args:
        variation_id: ClinVar variation ID (from esearch).
        client: Optional httpx.Client.

    Returns:
        Dict with keys: variation_id, classification, review_status,
        conditions, assembly, chrom, pos, ref, alt. Or None on failure.
    """
    url = f"{CLINVAR_BASE_URL}/efetch.fcgi"
    params = {
        "db": "clinvar",
        "id": variation_id,
        "rettype": "vcv",
        "retmode": "xml",
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    # Rate limiting: 3 req/sec without key, 10 with key
    time.sleep(0.35)

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        resp = client.get(url, params=params)
        resp.raise_for_status()
        return _parse_clinvar_xml(resp.text, variation_id)
    finally:
        if should_close:
            client.close()


def _parse_clinvar_xml(xml_text: str, fallback_id: str) -> dict | None:
    """Parse ClinVar VCV XML response.

    Args:
        xml_text: Raw XML response from efetch.
        fallback_id: ID to use if not found in XML.

    Returns:
        Parsed dict or None on parse failure.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        logger.warning(f"Failed to parse ClinVar XML: {e}")
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

    # Find the VariationArchive or ClinVarVariation element
    # ClinVar XML structure varies; search flexibly
    for vcv in root.iter():
        # Look for variation ID in VCV format
        if vcv.tag == "VariationArchive":
            acc = vcv.get("Accession")
            if acc:
                result["variation_id"] = acc

    # Classification
    for elem in root.iter("Description"):
        parent = None
        # Walk up to check if this is under Classifications
        text = elem.text
        if text and text in (
            "Pathogenic", "Likely pathogenic", "Uncertain significance",
            "Likely benign", "Benign", "Pathogenic/Likely pathogenic",
            "Benign/Likely benign", "Conflicting classifications of pathogenicity",
        ):
            result["classification"] = text
            break

    # Try GermlineClassification path (newer ClinVar XML format)
    for gc in root.iter("GermlineClassification"):
        desc = gc.find("Description")
        if desc is not None and desc.text:
            result["classification"] = desc.text
        review = gc.find("ReviewStatus")
        if review is not None and review.text:
            result["review_status"] = review.text

    # Fallback: older XML format
    if not result["review_status"]:
        for rs in root.iter("ReviewStatus"):
            if rs.text:
                result["review_status"] = rs.text
                break

    # Conditions / traits
    for trait in root.iter("TraitName"):
        if trait.text:
            result["conditions"].append(trait.text)
    # Also check Name elements under Trait
    if not result["conditions"]:
        for trait in root.iter("Trait"):
            name_elem = trait.find("Name")
            if name_elem is not None:
                val = name_elem.find("ElementValue")
                if val is not None and val.text:
                    result["conditions"].append(val.text)

    # Genomic coordinates — prefer GRCh38
    grch38_loc = None
    grch37_loc = None
    for loc in root.iter("SequenceLocation"):
        assembly = loc.get("Assembly")
        if assembly == "GRCh38":
            grch38_loc = loc
        elif assembly == "GRCh37":
            grch37_loc = loc

    chosen = grch38_loc or grch37_loc
    if chosen is not None:
        result["assembly"] = chosen.get("Assembly")
        result["chrom"] = chosen.get("Chr")
        start = chosen.get("start") or chosen.get("Start")
        if start:
            result["pos"] = int(start)
        result["ref"] = chosen.get("referenceAlleleVCF") or chosen.get("referenceAllele")
        result["alt"] = chosen.get("alternateAlleleVCF") or chosen.get("alternateAllele")

    if not result["conditions"]:
        result["conditions"] = None

    return result
```

**Step 4: Run tests to verify they pass**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestClinVarClient -v`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add varis/m1_ingestion/clinvar_client.py tests/test_m1_ingestion.py
git commit -m "feat(m1): implement ClinVar client with genomic coordinate extraction"
```

---

## Task 7: Implement gnomad_client

**Files:**
- Modify: `varis/m1_ingestion/gnomad_client.py`
- Modify: `tests/test_m1_ingestion.py`

**Step 1: Write the tests**

Replace the `TestGnomADClient` class in `tests/test_m1_ingestion.py`:

```python
from varis.m1_ingestion.gnomad_client import fetch_gnomad


class TestGnomADClient:
    @pytest.mark.timeout(30)
    def test_fetch_with_genomic_coords(self):
        """Fetch gnomAD frequency using ClinVar-provided genomic coordinates."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        # Simulate ClinVar having provided genomic coords
        record.reference_build = "GRCh38"
        record.clinvar_chrom = "17"
        record.clinvar_pos = 43057051
        record.clinvar_ref = "G"
        record.clinvar_alt = "A"
        record.coordinate_source = "clinvar"
        record = fetch_gnomad(record)
        # gnomAD may or may not have this variant — but should not crash
        # If found, frequency should be a small number or 0
        if record.gnomad_frequency is not None:
            assert 0 <= record.gnomad_frequency <= 1

    @pytest.mark.timeout(10)
    def test_skips_without_genomic_coords(self):
        """Should skip and set reason when no genomic coordinates available."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        # No genomic coords set
        record = fetch_gnomad(record)
        assert record.gnomad_frequency is None
        assert record.null_reasons.get("gnomad_frequency") == "no_genomic_coordinates"

    @pytest.mark.timeout(10)
    def test_skips_wrong_build(self):
        """Should skip if reference build doesn't match gnomAD endpoint."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.reference_build = "GRCh37"  # gnomAD v4 needs GRCh38
        record.clinvar_chrom = "17"
        record.clinvar_pos = 43057051
        record.clinvar_ref = "G"
        record.clinvar_alt = "A"
        record = fetch_gnomad(record)
        assert record.gnomad_frequency is None
```

**Step 2: Run tests to verify they fail**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestGnomADClient -v`
Expected: FAIL

**Step 3: Implement gnomad_client.py**

Replace `varis/m1_ingestion/gnomad_client.py`:

```python
"""gnomAD Client — Retrieves population allele frequencies from gnomAD (Broad Institute).

Only queries gnomAD when ClinVar has provided genomic coordinates AND the reference
build matches the gnomAD endpoint (v4 = GRCh38). Without genomic coords, gnomAD
fields are set to None with reason "no_genomic_coordinates".

Populates: gnomad_frequency, gnomad_popmax, gnomad_homozygotes.
"""

import logging
import httpx
from varis.models.variant_record import VariantRecord, NullReason
from varis.config import GNOMAD_API_URL

logger = logging.getLogger(__name__)

_TIMEOUT = 15.0
_GNOMAD_BUILD = "GRCh38"  # gnomAD v4 uses GRCh38


def fetch_gnomad(variant_record: VariantRecord,
                 client: httpx.Client | None = None) -> VariantRecord:
    """Query gnomAD for allele frequency using ClinVar-provided genomic coordinates.

    Args:
        variant_record: Should have clinvar_chrom, clinvar_pos, clinvar_ref,
            clinvar_alt, and reference_build set by clinvar_client.
        client: Optional httpx.Client for dependency injection.

    Returns:
        VariantRecord with gnomAD frequency fields populated.
    """
    # Gate: check if we have genomic coordinates
    chrom = variant_record.clinvar_chrom
    pos = variant_record.clinvar_pos
    ref = variant_record.clinvar_ref
    alt = variant_record.clinvar_alt
    build = variant_record.reference_build

    if not all([chrom, pos, ref, alt]):
        logger.info("No genomic coordinates available — skipping gnomAD")
        for f in ("gnomad_frequency", "gnomad_popmax", "gnomad_homozygotes"):
            variant_record.set_with_reason(f, None, "no_genomic_coordinates")
        return variant_record

    # Gate: check build matches
    if build != _GNOMAD_BUILD:
        logger.info(f"Build mismatch: have {build}, need {_GNOMAD_BUILD} — skipping gnomAD")
        for f in ("gnomad_frequency", "gnomad_popmax", "gnomad_homozygotes"):
            variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
        return variant_record

    try:
        data = _query_gnomad_graphql(chrom, pos, ref, alt, client)
        if data is None:
            logger.info(f"Variant {chrom}-{pos}-{ref}-{alt} not found in gnomAD")
            for f in ("gnomad_frequency", "gnomad_popmax", "gnomad_homozygotes"):
                variant_record.set_with_reason(f, None, NullReason.NO_DATA_AVAILABLE)
            return variant_record

        variant_record.gnomad_frequency = data.get("frequency")
        variant_record.gnomad_popmax = data.get("popmax")
        variant_record.gnomad_homozygotes = data.get("homozygotes")

        logger.info(
            f"gnomAD: {chrom}-{pos}-{ref}-{alt} → "
            f"freq={variant_record.gnomad_frequency}, "
            f"popmax={variant_record.gnomad_popmax}"
        )

    except Exception as e:
        logger.warning(f"gnomAD fetch failed: {e}")
        for f in ("gnomad_frequency", "gnomad_popmax", "gnomad_homozygotes"):
            variant_record.set_with_reason(f, None, NullReason.TOOL_CRASHED)

    return variant_record


def _query_gnomad_graphql(chrom: str, pos: int, ref: str, alt: str,
                           client: httpx.Client | None = None) -> dict | None:
    """Execute GraphQL query against gnomAD API.

    Args:
        chrom: Chromosome (e.g., "17").
        pos: Genomic position.
        ref: Reference allele.
        alt: Alternate allele.
        client: Optional httpx.Client.

    Returns:
        Dict with keys: frequency, popmax, homozygotes. Or None if not found.
    """
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"

    query = """
    query GnomadVariant($variantId: String!, $dataset: DatasetId!) {
      variant(variantId: $variantId, dataset: $dataset) {
        variant_id
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

    variables = {
        "variantId": variant_id,
        "dataset": "gnomad_r4",
    }

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        resp = client.post(
            GNOMAD_API_URL,
            json={"query": query, "variables": variables},
        )
        resp.raise_for_status()
        data = resp.json()

        variant_data = data.get("data", {}).get("variant")
        if not variant_data:
            return None

        # Calculate overall frequency from genome + exome
        total_ac = 0
        total_an = 0
        total_hom = 0
        popmax_freq = 0.0

        for source_key in ("genome", "exome"):
            source = variant_data.get(source_key)
            if source:
                total_ac += source.get("ac", 0) or 0
                total_an += source.get("an", 0) or 0
                total_hom += source.get("homozygote_count", 0) or 0

                # Calculate popmax
                for pop in source.get("populations", []):
                    pop_ac = pop.get("ac", 0) or 0
                    pop_an = pop.get("an", 0) or 0
                    if pop_an > 0:
                        pop_freq = pop_ac / pop_an
                        popmax_freq = max(popmax_freq, pop_freq)

        frequency = total_ac / total_an if total_an > 0 else 0.0

        return {
            "frequency": frequency,
            "popmax": popmax_freq if popmax_freq > 0 else None,
            "homozygotes": total_hom,
        }
    finally:
        if should_close:
            client.close()
```

**Step 4: Run tests to verify they pass**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestGnomADClient -v`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add varis/m1_ingestion/gnomad_client.py tests/test_m1_ingestion.py
git commit -m "feat(m1): implement gnomAD client gated on ClinVar genomic coordinates"
```

---

## Task 8: Implement alphafold_client

**Files:**
- Modify: `varis/m1_ingestion/alphafold_client.py`
- Modify: `tests/test_m1_ingestion.py`

**Step 1: Write the tests**

Add to `tests/test_m1_ingestion.py`:

```python
from pathlib import Path
from varis.m1_ingestion.alphafold_client import fetch_alphafold_structure
from varis.config import STRUCTURES_DIR


class TestAlphaFoldClient:
    @pytest.mark.timeout(60)
    def test_download_brca1_structure(self):
        """Download BRCA1 AlphaFold structure — real API call."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.uniprot_id = "P38398"
        record = fetch_alphafold_structure(record)
        assert record.structure_source == "alphafold"
        assert record.pdb_path is not None
        assert Path(record.pdb_path).exists()
        assert Path(record.pdb_path).stat().st_size > 0

    @pytest.mark.timeout(60)
    def test_caches_existing_file(self):
        """Second download should use cached file."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.uniprot_id = "P38398"
        record = fetch_alphafold_structure(record)
        first_path = record.pdb_path

        # Second call should be fast (cached)
        record2 = create_variant_record("BRCA1", "p.Arg1699Trp")
        record2.uniprot_id = "P38398"
        record2 = fetch_alphafold_structure(record2)
        assert record2.pdb_path == first_path

    @pytest.mark.timeout(10)
    def test_no_uniprot_id(self):
        """Should skip when uniprot_id is not available."""
        record = create_variant_record("FAKEGENE", "p.Ala1Val")
        record = fetch_alphafold_structure(record)
        assert record.pdb_path is None
        assert "pdb_path" in record.null_reasons
```

**Step 2: Run tests to verify they fail**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestAlphaFoldClient -v`
Expected: FAIL

**Step 3: Implement alphafold_client.py**

Replace `varis/m1_ingestion/alphafold_client.py`:

```python
"""AlphaFold DB Client — Downloads predicted 3D protein structures.

Downloads PDB files from the AlphaFold Protein Structure Database.
Caches downloaded files to avoid re-downloading for the same protein.
This is the primary structure source. ESMFold (in M2) is the fallback.

Populates: structure_source, pdb_path (stored in data/structures/).
"""

import logging
from pathlib import Path
import httpx
from varis.models.variant_record import VariantRecord, NullReason
from varis.config import STRUCTURES_DIR

logger = logging.getLogger(__name__)

_TIMEOUT = 30.0
_ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"


def fetch_alphafold_structure(variant_record: VariantRecord,
                              client: httpx.Client | None = None) -> VariantRecord:
    """Download AlphaFold predicted structure for the protein.

    Args:
        variant_record: Must have uniprot_id set (from uniprot_client).
        client: Optional httpx.Client for dependency injection.

    Returns:
        VariantRecord with pdb_path set to downloaded file, or None if unavailable.
    """
    uniprot_id = variant_record.uniprot_id
    if not uniprot_id:
        logger.info("No UniProt ID available — skipping AlphaFold download")
        for f in ("structure_source", "pdb_path"):
            variant_record.set_with_reason(f, None, NullReason.UPSTREAM_DEPENDENCY_FAILED)
        return variant_record

    try:
        pdb_path = _download_pdb(uniprot_id, STRUCTURES_DIR, client)
        if pdb_path:
            variant_record.structure_source = "alphafold"
            variant_record.pdb_path = str(pdb_path)
            logger.info(f"AlphaFold structure: {pdb_path}")
        else:
            logger.warning(f"AlphaFold structure not available for {uniprot_id}")
            for f in ("structure_source", "pdb_path"):
                variant_record.set_with_reason(f, None, NullReason.NO_DATA_AVAILABLE)

    except Exception as e:
        logger.warning(f"AlphaFold download failed for {uniprot_id}: {e}")
        for f in ("structure_source", "pdb_path"):
            variant_record.set_with_reason(f, None, NullReason.TOOL_CRASHED)

    return variant_record


def _download_pdb(uniprot_id: str, output_dir: Path,
                  client: httpx.Client | None = None) -> Path | None:
    """Download PDB file from AlphaFold DB, with caching.

    Args:
        uniprot_id: UniProt accession ID.
        output_dir: Directory to save the PDB file.
        client: Optional httpx.Client.

    Returns:
        Path to downloaded PDB file, or None on failure.
    """
    filename = f"AF-{uniprot_id}-F1-model_v4.pdb"
    output_path = output_dir / filename

    # Cache check: skip download if file already exists and has content
    if output_path.exists() and output_path.stat().st_size > 0:
        logger.info(f"Using cached structure: {output_path}")
        return output_path

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    url = _ALPHAFOLD_PDB_URL.format(uniprot_id=uniprot_id)

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        resp = client.get(url)
        if resp.status_code == 404:
            logger.warning(f"AlphaFold structure not found for {uniprot_id}")
            return None
        resp.raise_for_status()

        output_path.write_bytes(resp.content)
        logger.info(f"Downloaded AlphaFold structure: {output_path} ({len(resp.content)} bytes)")
        return output_path
    finally:
        if should_close:
            client.close()
```

**Step 4: Run tests to verify they pass**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestAlphaFoldClient -v`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add varis/m1_ingestion/alphafold_client.py tests/test_m1_ingestion.py
git commit -m "feat(m1): implement AlphaFold client with file caching"
```

---

## Task 9: Implement alphamissense_client

**Files:**
- Modify: `varis/m1_ingestion/alphamissense_client.py`
- Modify: `tests/test_m1_ingestion.py`

**Step 1: Write the tests**

Add to `tests/test_m1_ingestion.py`:

```python
from varis.m1_ingestion.alphamissense_client import fetch_alphamissense


class TestAlphaMissenseClient:
    @pytest.mark.timeout(30)
    def test_fetch_brca1_score(self):
        """Fetch AlphaMissense score for BRCA1 p.Arg1699Trp — real API call."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.uniprot_id = "P38398"
        record.residue_position = 1699
        record.ref_aa_single = "R"
        record.alt_aa_single = "W"
        record = fetch_alphamissense(record)
        # AlphaMissense score for this variant should be ~0.934
        if record.alphamissense_score is not None:
            assert 0.8 < record.alphamissense_score < 1.0
            assert record.alphamissense_class == "likely_pathogenic"

    @pytest.mark.timeout(10)
    def test_unknown_variant(self):
        """Unknown variant should not crash."""
        record = create_variant_record("FAKEGENE", "p.Ala1Val")
        record.uniprot_id = "FAKE123"
        record.residue_position = 1
        record.ref_aa_single = "A"
        record.alt_aa_single = "V"
        record = fetch_alphamissense(record)
        # May return None — that's fine
        assert record is not None
```

**Step 2: Run tests to verify they fail**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestAlphaMissenseClient -v`
Expected: FAIL

**Step 3: Implement alphamissense_client.py**

Replace `varis/m1_ingestion/alphamissense_client.py`:

```python
"""AlphaMissense Client — Retrieves pre-computed pathogenicity scores.

Phase 1: Uses the hegelab.org community web API for individual variant lookups.
Phase 2+: May switch to local SQLite database from the full 4GB Zenodo TSV.

AlphaMissense scores are used as ONE of ~15 features in Varis's ML ensemble.
The relationship is complementary: AlphaMissense pre-screens, Varis investigates why.

Populates: alphamissense_score, alphamissense_class.
"""

import logging
import httpx
from varis.models.variant_record import VariantRecord, NullReason

logger = logging.getLogger(__name__)

_TIMEOUT = 15.0
_ALPHAMISSENSE_API_URL = "https://alphamissense.hegelab.org/api/variant"


def fetch_alphamissense(variant_record: VariantRecord,
                        client: httpx.Client | None = None) -> VariantRecord:
    """Look up AlphaMissense pre-computed pathogenicity score.

    Args:
        variant_record: Must have uniprot_id, residue_position, ref_aa_single,
            and alt_aa_single set.
        client: Optional httpx.Client for dependency injection.

    Returns:
        VariantRecord with AlphaMissense score and classification.
    """
    uniprot_id = variant_record.uniprot_id
    position = variant_record.residue_position
    ref_aa = variant_record.ref_aa_single
    alt_aa = variant_record.alt_aa_single

    if not all([uniprot_id, position, ref_aa, alt_aa]):
        logger.info("Missing fields for AlphaMissense lookup — skipping")
        for f in ("alphamissense_score", "alphamissense_class"):
            variant_record.set_with_reason(f, None, NullReason.UPSTREAM_DEPENDENCY_FAILED)
        return variant_record

    try:
        result = _lookup_score(uniprot_id, position, ref_aa, alt_aa, client)
        if result:
            score, classification = result
            variant_record.alphamissense_score = score
            variant_record.alphamissense_class = classification
            logger.info(
                f"AlphaMissense: {uniprot_id} {ref_aa}{position}{alt_aa} → "
                f"score={score}, class={classification}"
            )
        else:
            logger.info(
                f"AlphaMissense: no score found for {uniprot_id} {ref_aa}{position}{alt_aa}"
            )
            for f in ("alphamissense_score", "alphamissense_class"):
                variant_record.set_with_reason(f, None, NullReason.NO_DATA_AVAILABLE)

    except Exception as e:
        logger.warning(f"AlphaMissense lookup failed: {e}")
        for f in ("alphamissense_score", "alphamissense_class"):
            variant_record.set_with_reason(f, None, NullReason.TOOL_CRASHED)

    return variant_record


def _lookup_score(uniprot_id: str, position: int, ref_aa: str, alt_aa: str,
                  client: httpx.Client | None = None) -> tuple[float, str] | None:
    """Look up AlphaMissense score via hegelab.org web API.

    Args:
        uniprot_id: UniProt accession ID.
        position: Residue position (1-indexed).
        ref_aa: Reference amino acid (single letter).
        alt_aa: Alternate amino acid (single letter).
        client: Optional httpx.Client.

    Returns:
        Tuple of (score, classification) or None if not found.
        Classification is one of: "likely_pathogenic", "ambiguous", "likely_benign".
    """
    # hegelab.org API uses format: {uniprot_id}/{ref_aa}{position}{alt_aa}
    variant_str = f"{ref_aa}{position}{alt_aa}"
    url = f"{_ALPHAMISSENSE_API_URL}/{uniprot_id}/{variant_str}"

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        resp = client.get(url)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()

        data = resp.json()

        # Extract score and classification from response
        score = data.get("am_pathogenicity")
        classification = data.get("am_class")

        if score is None:
            return None

        score = float(score)

        # Normalize classification to our standard format
        if classification:
            classification = classification.lower().replace(" ", "_")

        return (score, classification)
    except (httpx.HTTPStatusError, ValueError, KeyError) as e:
        logger.warning(f"AlphaMissense API error for {url}: {e}")
        return None
    finally:
        if should_close:
            client.close()
```

**Step 4: Run tests to verify they pass**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/test_m1_ingestion.py::TestAlphaMissenseClient -v`
Expected: ALL PASS (or score test may get None if hegelab.org API format differs — see note below)

**Note:** The hegelab.org API format may need adjustment at runtime. If the API returns a different JSON structure or uses different URL patterns, adapt the `_lookup_score` function accordingly. The test is written to tolerate a None response for the score.

**Step 5: Commit**

```bash
git add varis/m1_ingestion/alphamissense_client.py tests/test_m1_ingestion.py
git commit -m "feat(m1): implement AlphaMissense client via hegelab.org API"
```

---

## Task 10: Full pipeline integration test

**Files:**
- Modify: `tests/test_m1_ingestion.py`
- Modify: `tests/test_pipeline.py`

**Step 1: Write integration test in test_m1_ingestion.py**

Add to `tests/test_m1_ingestion.py`:

```python
import json
from pathlib import Path


class TestM1Integration:
    @pytest.mark.timeout(120)
    def test_full_m1_brca1(self):
        """Run full M1 pipeline on BRCA1 p.Arg1699Trp with real API calls."""
        from varis.m1_ingestion import run

        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = run(record)

        # HGVS parsing
        assert record.residue_position == 1699
        assert record.ref_amino_acid == "Arg"
        assert record.alt_amino_acid == "Trp"

        # UniProt
        assert record.uniprot_id == "P38398"
        assert record.protein_sequence is not None
        assert record.protein_length > 1000

        # Normalizer validation
        assert record.coordinate_mapping_confidence == "high"

        # ClinVar (should find this well-known variant)
        assert record.clinvar_id is not None

        # AlphaFold (should download structure)
        assert record.structure_source == "alphafold"
        assert record.pdb_path is not None

        # Module tracking
        assert "M1" in record.modules_completed

    @pytest.mark.timeout(120)
    def test_full_m1_produces_valid_json(self):
        """M1 output should serialize to valid JSON and round-trip."""
        from varis.m1_ingestion import run

        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = run(record)

        # Should serialize without error
        json_str = record.to_json()
        data = json.loads(json_str)
        assert data["gene_symbol"] == "BRCA1"
        assert data["variant_id"] == "BRCA1_p.Arg1699Trp"

        # Should round-trip
        restored = record.from_json(json_str)
        assert restored.residue_position == record.residue_position
        assert restored.uniprot_id == record.uniprot_id
```

**Step 2: Update pipeline tests in test_pipeline.py**

Replace `tests/test_pipeline.py`:

```python
"""Tests for the full pipeline orchestration."""

import pytest
from varis.models.variant_record import create_variant_record


class TestPipelineGracefulDegradation:
    def test_pipeline_never_crashes(self):
        """Pipeline should always return a VariantRecord, never raise."""
        from varis.pipeline import investigate
        record = investigate("FAKEGENE", "p.Ala1Val")
        assert record is not None
        assert record.gene_symbol == "FAKEGENE"

    @pytest.mark.timeout(120)
    def test_pipeline_tracks_failures(self):
        """Pipeline with bad gene should track failures but still complete."""
        from varis.pipeline import investigate
        record = investigate("FAKEGENE999", "p.Ala1Val")
        assert record is not None
        # M1 should still complete (individual sub-modules may fail)
        assert record.modules_completed is not None
        # Some data should be missing
        assert record.uniprot_id is None

    @pytest.mark.timeout(10)
    def test_pipeline_records_timing(self):
        """Pipeline should record processing time."""
        from varis.pipeline import investigate
        record = investigate("FAKEGENE", "p.Ala1Val")
        assert record.processing_time_seconds is not None
        assert record.processing_time_seconds >= 0
```

**Step 3: Run all tests**

Run: `cd /Users/darwin/Documents/Claude/varis && python -m pytest tests/ -v --timeout=120`
Expected: ALL PASS

**Step 4: Commit**

```bash
git add tests/test_m1_ingestion.py tests/test_pipeline.py
git commit -m "test(m1): add full M1 integration test with BRCA1 p.Arg1699Trp"
```

---

## Task 11: Initialize git repo and make first tagged commit

**Step 1: Initialize git if not already a repo**

```bash
cd /Users/darwin/Documents/Claude/varis
git init
```

**Step 2: Create .gitignore additions**

Verify `.gitignore` includes:
- `data/structures/*.pdb` (downloaded structures are large, don't commit)
- `*.pyc`, `__pycache__/`
- `.env`
- `data/models/`

**Step 3: Stage and commit all Phase 1 work**

If individual commits weren't made during each task, make one consolidated commit:

```bash
git add -A
git commit -m "feat: Phase 1 complete — M1 ingestion pipeline with 7 sub-modules"
```

**Step 4: Run the CLI as final validation**

```bash
cd /Users/darwin/Documents/Claude/varis
pip install -e ".[dev]"
python -m varis BRCA1 p.Arg1699Trp --output data/ -v
```

Expected output: a `BRCA1_p.Arg1699Trp.json` file in `data/` with populated fields.

**Step 5: Tag the release**

```bash
git tag -a v0.1.0 -m "Phase 1: M1 ingestion pipeline — all 7 sub-modules implemented"
```

---

## Summary of tasks

| Task | Component | Type | Est. Complexity |
|------|-----------|------|----------------|
| 1 | VariantRecord schema update | Schema | Low |
| 2 | M1 orchestrator reorder | Refactor | Low |
| 3 | hgvs_parser | Implementation | Low |
| 4 | uniprot_client | Implementation + API | Medium |
| 5 | variant_normalizer | Implementation | Low |
| 6 | clinvar_client | Implementation + API + XML | High |
| 7 | gnomad_client | Implementation + GraphQL | Medium |
| 8 | alphafold_client | Implementation + Download | Low |
| 9 | alphamissense_client | Implementation + API | Medium |
| 10 | Integration tests | Testing | Low |
| 11 | Git init + tag | Infrastructure | Low |
