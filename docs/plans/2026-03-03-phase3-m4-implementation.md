# Phase 3: M4 Conservation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement evolutionary conservation analysis (M4) that fetches orthologs, aligns them, computes Shannon entropy, and falls back to ConSurf when the primary pipeline fails.

**Architecture:** Three-stage pipeline (UniProt orthologs → Clustal Omega MSA → entropy scoring) with ConSurf fallback and protein-level caching. Each stage is a thin wrapper in its own file, matching the M1/M2/M3 pattern.

**Tech Stack:** httpx (UniProt API, EBI Clustal Omega API, ConSurf API), BioPython (alignment parsing), math (entropy), pytest (testing)

---

### Task 1: Update VariantRecord Schema to v1.3.0

**Files:**
- Modify: `varis/models/variant_record.py`
- Modify: `tests/conftest.py`
- Create: `tests/test_m4_conservation.py`

**Step 1: Write failing test**

Create `tests/test_m4_conservation.py`:

```python
"""Tests for M4: Conservation Analysis."""
import pytest
import math
from unittest.mock import MagicMock, patch
from varis.models.variant_record import (
    VariantRecord, create_variant_record, NullReason, RECORD_SCHEMA_VERSION,
)


class TestSchemaV130:
    """Verify schema v1.3.0 fields exist on VariantRecord."""

    def test_schema_version_is_1_3_0(self):
        assert RECORD_SCHEMA_VERSION == "1.3.0"

    def test_new_msa_fields_exist(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "msa_num_sequences")
        assert hasattr(record, "msa_gap_fraction_at_site")
        assert hasattr(record, "msa_column_index")

    def test_insufficient_data_reason_exists(self):
        assert hasattr(NullReason, "INSUFFICIENT_DATA")
        assert NullReason.INSUFFICIENT_DATA == "insufficient_data"
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_m4_conservation.py::TestSchemaV130 -v`
Expected: FAIL — schema version is "1.2.0", new fields don't exist

**Step 3: Update VariantRecord**

In `varis/models/variant_record.py`:

1. Change `RECORD_SCHEMA_VERSION = "1.2.0"` → `"1.3.0"`

2. Add `INSUFFICIENT_DATA = "insufficient_data"` to the `NullReason` class.

3. Add new MSA fields after the existing CONSERVATION section:
```python
    # CONSERVATION — Evolutionary analysis (populated by M4)
    conservation_score: Optional[float] = None
    conservation_method: Optional[str] = None
    num_orthologs: Optional[int] = None
    position_entropy: Optional[float] = None
    conserved_across_mammals: Optional[bool] = None
    msa_num_sequences: Optional[int] = None
    msa_gap_fraction_at_site: Optional[float] = None
    msa_column_index: Optional[int] = None
```

4. Update `tests/conftest.py` `fully_populated_record` fixture — add:
```python
    record.msa_num_sequences = 46
    record.msa_gap_fraction_at_site = 0.02
    record.msa_column_index = 1699
```

**Step 4: Run tests**

Run: `pytest tests/test_m4_conservation.py::TestSchemaV130 -v`
Expected: PASS

Run: `pytest tests/ -v --tb=short`
Expected: ALL PASS

**Step 5: Commit**

```
refactor(schema): update VariantRecord to v1.3.0 for Phase 3

Add MSA metadata fields: msa_num_sequences, msa_gap_fraction_at_site,
msa_column_index. Add NullReason.INSUFFICIENT_DATA.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 2: Conservation Scorer (Pure Logic, No APIs)

**Files:**
- Replace: `varis/m4_conservation/conservation_scorer.py`
- Modify: `tests/test_m4_conservation.py`

This is the core logic with no external dependencies — implement it first so other modules have something to feed into.

**Step 1: Write failing tests**

Add to `tests/test_m4_conservation.py`:

```python
class TestConservationScorer:
    """Tests for conservation_scorer.py — entropy and position mapping."""

    def test_entropy_fully_conserved(self):
        """All same AA → entropy=0, score=1.0."""
        from varis.m4_conservation.conservation_scorer import _shannon_entropy
        column = ["R"] * 20
        assert _shannon_entropy(column) == pytest.approx(0.0, abs=0.001)

    def test_entropy_maximally_variable(self):
        """All different AAs → entropy ≈ log2(20), score ≈ 0.0."""
        from varis.m4_conservation.conservation_scorer import _shannon_entropy
        aas = list("ACDEFGHIKLMNPQRSTVWY")
        assert _shannon_entropy(aas) == pytest.approx(math.log2(20), abs=0.01)

    def test_entropy_ignores_gaps(self):
        """Gaps excluded from frequency, but present in input."""
        from varis.m4_conservation.conservation_scorer import _shannon_entropy
        column = ["R"] * 18 + ["-", "-"]
        # Should compute entropy from 18 R's only → fully conserved
        assert _shannon_entropy(column) == pytest.approx(0.0, abs=0.001)

    def test_position_mapping_no_gaps(self):
        """Direct mapping when query has no gaps."""
        from varis.m4_conservation.conservation_scorer import _map_position_to_column
        # Query row: "MKRST" (no gaps), position 3 (K=1, R=2, S=3 → 0-indexed col 2)
        query_row = "MKRST"
        col = _map_position_to_column(query_row, position=3)
        assert col == 2  # 0-indexed: M=0, K=1, R=2 is pos 3... wait
        # Actually position is 1-indexed: M=pos1, K=pos2, R=pos3
        # So column index = 2 (0-indexed)
        assert col == 2

    def test_position_mapping_with_gaps(self):
        """Gaps in query shift the column index."""
        from varis.m4_conservation.conservation_scorer import _map_position_to_column
        # Query row: "M-KR-ST" → M=pos1(col0), K=pos2(col2), R=pos3(col3), S=pos4(col5), T=pos5(col6)
        query_row = "M-KR-ST"
        col = _map_position_to_column(query_row, position=3)
        assert col == 3  # R is at position 3, alignment column 3

    def test_position_mapping_validates_ref_aa(self):
        """Returns None if ref AA doesn't match."""
        from varis.m4_conservation.conservation_scorer import _map_position_to_column
        query_row = "MKRST"
        col = _map_position_to_column(query_row, position=3, expected_aa="W")
        assert col is None  # R != W → mismatch

    def test_score_conservation_full(self, m1_completed_record):
        """Full scorer with a small test alignment."""
        from varis.m4_conservation.conservation_scorer import score_conservation
        # Build a simple alignment: query + 10 orthologs, all have R at position 3
        alignment = {
            "sequences": {
                "query": "MKRST",
                "orth1": "MKRST", "orth2": "MKRST", "orth3": "MKRST",
                "orth4": "MKRST", "orth5": "MKRST", "orth6": "MKRST",
                "orth7": "MKRST", "orth8": "MKRST", "orth9": "MKRST",
                "orth10": "MKRST",
            },
            "query_id": "query",
            "taxonomy": {
                "orth1": 9606, "orth2": 9615, "orth3": 10090,
                "orth4": 9913, "orth5": 9823,  # 5 mammals
                "orth6": 7955, "orth7": 8364, "orth8": 9031,
                "orth9": 28377, "orth10": 13616,
            },
        }
        m1_completed_record.residue_position = 3
        m1_completed_record.ref_aa_single = "R"
        result = score_conservation(m1_completed_record, alignment)
        assert result.conservation_score == pytest.approx(1.0, abs=0.01)
        assert result.position_entropy == pytest.approx(0.0, abs=0.01)
        assert result.msa_column_index == 2
        assert result.conservation_available is True

    def test_mammal_conservation_threshold(self, m1_completed_record):
        """≥90% mammals have ref AA → True."""
        from varis.m4_conservation.conservation_scorer import score_conservation
        alignment = {
            "sequences": {
                "query": "MKRST",
                "m1": "MKRST", "m2": "MKRST", "m3": "MKRST",
                "m4": "MKRST", "m5": "MKRST",  # All 5 mammals have R
                "o1": "MKWST",  # non-mammal with different AA
            },
            "query_id": "query",
            "taxonomy": {
                "m1": 9606, "m2": 9615, "m3": 10090, "m4": 9913, "m5": 9823,
                "o1": 7955,
            },
        }
        m1_completed_record.residue_position = 3
        m1_completed_record.ref_aa_single = "R"
        result = score_conservation(m1_completed_record, alignment)
        assert result.conserved_across_mammals is True

    def test_mammal_conservation_too_few(self, m1_completed_record):
        """<5 mammals → conserved_across_mammals=None."""
        from varis.m4_conservation.conservation_scorer import score_conservation
        alignment = {
            "sequences": {
                "query": "MKRST",
                "m1": "MKRST", "m2": "MKRST",  # Only 2 mammals
                "o1": "MKRST", "o2": "MKRST", "o3": "MKRST",
                "o4": "MKRST", "o5": "MKRST", "o6": "MKRST",
                "o7": "MKRST", "o8": "MKRST",
            },
            "query_id": "query",
            "taxonomy": {"m1": 9606, "m2": 9615, "o1": 7955, "o2": 8364,
                          "o3": 9031, "o4": 28377, "o5": 13616, "o6": 7719,
                          "o7": 6239, "o8": 7227},
        }
        m1_completed_record.residue_position = 3
        m1_completed_record.ref_aa_single = "R"
        result = score_conservation(m1_completed_record, alignment)
        assert result.conserved_across_mammals is None
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m4_conservation.py::TestConservationScorer -v`
Expected: FAIL

**Step 3: Implement conservation_scorer.py**

Replace `varis/m4_conservation/conservation_scorer.py`:

```python
"""Conservation Scorer — Computes per-position conservation from MSA.

Calculates Shannon entropy at the mutation position and determines
mammalian conservation. Entropy is normalized by log2(20) so that
conservation_score is comparable across proteins.

Definitions:
  H = -Σ p(aa) * log2(p(aa)) over 20 standard amino acids (gaps excluded)
  conservation_score = 1 - (H / log2(20)), range [0.0, 1.0]
  conserved_across_mammals: ref AA in ≥90% mammalian sequences (≥5 required)

Populates: conservation_score, position_entropy, conservation_method,
           msa_column_index, msa_gap_fraction_at_site, msa_num_sequences,
           conserved_across_mammals, num_orthologs.
Sets: conservation_available, conservation_missing_reason.
"""

import logging
import math
from collections import Counter

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

_STANDARD_AAS = set("ACDEFGHIKLMNPQRSTVWY")
_MAX_ENTROPY = math.log2(20)  # ≈ 4.322

# NCBI taxonomy IDs for Mammalia (class-level taxon 40674)
# This is a representative set of common mammalian taxon IDs.
# In production, check if taxon_id descends from 40674.
_MAMMAL_TAXON_IDS = {
    9606,   # Homo sapiens
    9598,   # Pan troglodytes
    9544,   # Macaca mulatta
    9615,   # Canis lupus familiaris
    9913,   # Bos taurus
    9823,   # Sus scrofa
    10090,  # Mus musculus
    10116,  # Rattus norvegicus
    9986,   # Oryctolagus cuniculus
    9685,   # Felis catus
    9796,   # Equus caballus
    13616,  # Monodelphis domestica
    9258,   # Ornithorhynchus anatinus
    9669,   # Mustela putorius furo
    9940,   # Ovis aries
    30611,  # Pan paniscus
    9601,   # Pongo abelii
}

_MIN_MAMMALS = 5
_MAMMAL_CONSERVATION_THRESHOLD = 0.9


def score_conservation(variant_record: VariantRecord,
                       alignment: dict) -> VariantRecord:
    """Calculate conservation score at the mutation position from MSA.

    Args:
        variant_record: Must have residue_position and ref_aa_single set.
        alignment: Dict with keys:
            "sequences": {id: sequence_string, ...} (aligned, with gaps)
            "query_id": str identifying the query sequence
            "taxonomy": {id: taxon_id, ...} (optional, for mammal check)

    Returns:
        VariantRecord with conservation fields populated.
    """
    position = variant_record.residue_position
    ref_aa = variant_record.ref_aa_single
    sequences = alignment.get("sequences", {})
    query_id = alignment.get("query_id", "query")
    taxonomy = alignment.get("taxonomy", {})

    if not sequences or query_id not in sequences:
        logger.warning("No alignment data or query not found in alignment")
        variant_record.set_feature_status(
            "conservation", False, NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        return variant_record

    try:
        query_row = sequences[query_id]
        num_seqs = len(sequences)
        variant_record.msa_num_sequences = num_seqs
        variant_record.num_orthologs = num_seqs - 1  # exclude query

        # Map protein position to MSA column
        col_idx = _map_position_to_column(query_row, position, expected_aa=ref_aa)
        if col_idx is None:
            logger.warning(
                "Could not map position %d to MSA column (ref AA mismatch or out of range)",
                position,
            )
            variant_record.set_feature_status(
                "conservation", False, NullReason.COORDINATE_MAPPING_FAILED
            )
            return variant_record

        variant_record.msa_column_index = col_idx

        # Extract column from all sequences
        column = []
        for seq_id, seq in sequences.items():
            if col_idx < len(seq):
                column.append(seq[col_idx])
            else:
                column.append("-")

        # Gap fraction
        gap_count = sum(1 for c in column if c == "-" or c == ".")
        variant_record.msa_gap_fraction_at_site = round(gap_count / len(column), 4)

        # Shannon entropy (gaps excluded)
        entropy = _shannon_entropy(column)
        variant_record.position_entropy = round(entropy, 6)
        variant_record.conservation_score = round(
            1.0 - (entropy / _MAX_ENTROPY), 4
        )
        variant_record.conservation_method = "clustal_omega"

        # Mammal conservation
        variant_record.conserved_across_mammals = _mammal_conservation(
            sequences, taxonomy, col_idx, ref_aa
        )

        variant_record.set_feature_status("conservation", True)
        logger.info(
            "Conservation: position %d entropy=%.3f score=%.3f mammals=%s",
            position, entropy, variant_record.conservation_score,
            variant_record.conserved_across_mammals,
        )

    except Exception as e:
        logger.warning("Conservation scoring failed: %s", e)
        variant_record.set_feature_status(
            "conservation", False, NullReason.TOOL_CRASHED
        )

    return variant_record


def _map_position_to_column(query_row: str, position: int,
                             expected_aa: str | None = None) -> int | None:
    """Map a 1-indexed protein position to a 0-indexed MSA column.

    Walks the query row, counting non-gap characters. When count equals
    position, returns that column index. Optionally validates the amino acid.

    Args:
        query_row: The aligned query sequence (with gaps).
        position: 1-indexed residue position in the unaligned protein.
        expected_aa: If provided, verify the AA at the mapped position matches.

    Returns:
        0-indexed column index, or None if position out of range or AA mismatch.
    """
    residue_count = 0
    for col_idx, char in enumerate(query_row):
        if char != "-" and char != ".":
            residue_count += 1
            if residue_count == position:
                if expected_aa and char.upper() != expected_aa.upper():
                    logger.warning(
                        "Ref AA mismatch at position %d: expected %s, found %s",
                        position, expected_aa, char,
                    )
                    return None
                return col_idx
    return None  # Position beyond sequence length


def _shannon_entropy(column: list[str]) -> float:
    """Calculate Shannon entropy for an MSA column.

    Gaps ("-" and ".") are excluded from frequency calculation.
    Entropy is in bits (log base 2). Range: [0, log2(20)].

    Args:
        column: List of single-character amino acids (may include gaps).

    Returns:
        Shannon entropy in bits. 0.0 = fully conserved.
    """
    # Filter out gaps
    residues = [c.upper() for c in column if c in _STANDARD_AAS or c.upper() in _STANDARD_AAS]
    if not residues:
        return 0.0

    total = len(residues)
    counts = Counter(residues)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


def _mammal_conservation(sequences: dict, taxonomy: dict,
                          col_idx: int, ref_aa: str) -> bool | None:
    """Check if the reference AA is conserved across mammalian orthologs.

    Args:
        sequences: {seq_id: aligned_sequence, ...}
        taxonomy: {seq_id: ncbi_taxon_id, ...}
        col_idx: 0-indexed MSA column to check.
        ref_aa: Expected reference amino acid (single letter).

    Returns:
        True if ≥90% of mammals have ref AA, False if not, None if <5 mammals.
    """
    mammal_aas = []
    for seq_id, seq in sequences.items():
        taxon = taxonomy.get(seq_id)
        if taxon and taxon in _MAMMAL_TAXON_IDS:
            if col_idx < len(seq):
                aa = seq[col_idx]
                if aa != "-" and aa != ".":  # Ignore gaps
                    mammal_aas.append(aa.upper())

    if len(mammal_aas) < _MIN_MAMMALS:
        return None

    ref_count = sum(1 for aa in mammal_aas if aa == ref_aa.upper())
    fraction = ref_count / len(mammal_aas)
    return fraction >= _MAMMAL_CONSERVATION_THRESHOLD
```

**Step 4: Run tests**

Run: `pytest tests/test_m4_conservation.py::TestConservationScorer -v`
Expected: PASS

**Step 5: Commit**

```
feat(m4): implement conservation scorer with entropy and position mapping

Shannon entropy normalized by log2(20). Gap-aware MSA position mapping.
Mammal conservation with ≥5 threshold and ≥90% identity requirement.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 3: UniProt Orthologs Client

**Files:**
- Create: `varis/m4_conservation/uniprot_orthologs.py`
- Modify: `tests/test_m4_conservation.py`

**Step 1: Write failing tests**

Add to `tests/test_m4_conservation.py`:

```python
class TestUniProtOrthologs:
    """Tests for uniprot_orthologs.py — fetch ortholog sequences."""

    def test_fetch_orthologs_brca1(self, m1_completed_record):
        """Mocked: returns sequences with taxonomy."""
        from varis.m4_conservation.uniprot_orthologs import fetch_orthologs
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = (
            ">sp|P38398|BRCA1_HUMAN\nMKRST\n"
            ">sp|Q9GKK4|BRCA1_PANTR\nMKRST\n"
            * 10  # 10 orthologs + query
        )
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        record, orthologs = fetch_orthologs(m1_completed_record, client=mock_client)
        assert orthologs is not None
        assert len(orthologs) >= 10

    def test_fetch_orthologs_too_few(self, m1_completed_record):
        """<10 sequences → returns None."""
        from varis.m4_conservation.uniprot_orthologs import fetch_orthologs
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = ">sp|P38398|BRCA1_HUMAN\nMKRST\n"  # Only 1
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        record, orthologs = fetch_orthologs(m1_completed_record, client=mock_client)
        assert orthologs is None

    def test_fetch_orthologs_no_uniprot(self, m1_completed_record):
        """No uniprot_id → skip."""
        from varis.m4_conservation.uniprot_orthologs import fetch_orthologs
        m1_completed_record.uniprot_id = None
        record, orthologs = fetch_orthologs(m1_completed_record)
        assert orthologs is None

    def test_fetch_orthologs_caps_at_100(self, m1_completed_record):
        """Max 100 sequences returned."""
        from varis.m4_conservation.uniprot_orthologs import fetch_orthologs
        # Build 150 FASTA entries
        fasta = "".join(f">sp|P{i:05d}|ORTH{i}\nMKRST\n" for i in range(150))
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = fasta
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        record, orthologs = fetch_orthologs(m1_completed_record, client=mock_client)
        assert orthologs is not None
        assert len(orthologs["sequences"]) <= 101  # 100 orthologs + query
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m4_conservation.py::TestUniProtOrthologs -v`
Expected: FAIL

**Step 3: Implement uniprot_orthologs.py**

Create `varis/m4_conservation/uniprot_orthologs.py`. Key points:
- Query UniProt REST API: `https://rest.uniprot.org/uniref/stream?query=uniprot_id:{uniprot_id}&format=fasta`
- Or use UniProt's cluster endpoint to get UniRef50/90 members
- Parse FASTA response, extract sequence IDs and taxonomy from headers
- Cap at 100 sequences, require minimum 10
- Accept optional `client: httpx.Client | None` for DI
- Return `(variant_record, orthologs_dict | None)` where orthologs_dict has keys: `"sequences"`, `"query_id"`, `"taxonomy"`
- Cache to `data/conservation/{uniprot_id}_orthologs.json`

**Step 4: Run tests**

Run: `pytest tests/test_m4_conservation.py::TestUniProtOrthologs -v`
Expected: PASS

**Step 5: Commit**

```
feat(m4): implement UniProt ortholog fetcher

Queries UniProt REST API for ortholog sequences. Caps at 100,
requires minimum 10. Extracts taxonomy from FASTA headers.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 4: Clustal Omega Client

**Files:**
- Replace: `varis/m4_conservation/clustal_client.py` (rename from alignment.py)
- Modify: `tests/test_m4_conservation.py`

**Step 1: Write failing tests**

Add to `tests/test_m4_conservation.py`:

```python
class TestClustalClient:
    """Tests for clustal_client.py — EBI Clustal Omega API."""

    def test_clustal_alignment(self, m1_completed_record):
        """Mocked: submit + poll returns valid alignment."""
        from varis.m4_conservation.clustal_client import run_alignment
        orthologs = {
            "sequences": {
                "query": "MKRST",
                "orth1": "MKRST",
                "orth2": "MKRAT",
            },
            "query_id": "query",
            "taxonomy": {"orth1": 9606, "orth2": 10090},
        }
        # Mock submit response (job ID)
        mock_submit = MagicMock()
        mock_submit.status_code = 200
        mock_submit.text = "clustalo-R20260303-123456"
        # Mock status response (FINISHED)
        mock_status = MagicMock()
        mock_status.status_code = 200
        mock_status.text = "FINISHED"
        # Mock result response (aligned FASTA)
        mock_result = MagicMock()
        mock_result.status_code = 200
        mock_result.text = ">query\nMKRST\n>orth1\nMKRST\n>orth2\nMKRAT\n"

        mock_client = MagicMock()
        mock_client.post.return_value = mock_submit
        mock_client.get.side_effect = [mock_status, mock_result]

        record, alignment = run_alignment(m1_completed_record, orthologs, client=mock_client)
        assert alignment is not None
        assert "sequences" in alignment
        assert len(alignment["sequences"]) == 3

    def test_clustal_timeout(self, m1_completed_record):
        """Poll exceeds max retries → TIMED_OUT."""
        from varis.m4_conservation.clustal_client import run_alignment
        orthologs = {
            "sequences": {"query": "MKRST", "orth1": "MKRST"},
            "query_id": "query",
            "taxonomy": {},
        }
        mock_submit = MagicMock()
        mock_submit.status_code = 200
        mock_submit.text = "clustalo-R20260303-123456"
        mock_status = MagicMock()
        mock_status.status_code = 200
        mock_status.text = "RUNNING"  # Never finishes

        mock_client = MagicMock()
        mock_client.post.return_value = mock_submit
        mock_client.get.return_value = mock_status

        record, alignment = run_alignment(
            m1_completed_record, orthologs, client=mock_client, max_polls=2, poll_interval=0
        )
        assert alignment is None

    def test_clustal_no_sequences(self, m1_completed_record):
        """No sequences → skip."""
        from varis.m4_conservation.clustal_client import run_alignment
        record, alignment = run_alignment(m1_completed_record, None)
        assert alignment is None
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m4_conservation.py::TestClustalClient -v`
Expected: FAIL

**Step 3: Implement clustal_client.py**

Create `varis/m4_conservation/clustal_client.py`. Key points:
- EBI Clustal Omega API: submit sequences as FASTA to `https://www.ebi.ac.uk/Tools/services/rest/clustalo/run`
- Poll status at `https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}`
- Get result from `https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-fasta`
- Accept `max_polls` (default 30) and `poll_interval` (default 5s) parameters
- Parse aligned FASTA result back into the orthologs dict format (update sequences with aligned versions)
- Preserve taxonomy from input orthologs dict
- Accept optional `client: httpx.Client | None` for DI

Also delete the old `varis/m4_conservation/alignment.py` stub (replaced by clustal_client.py).

**Step 4: Run tests**

Run: `pytest tests/test_m4_conservation.py::TestClustalClient -v`
Expected: PASS

**Step 5: Commit**

```
feat(m4): implement Clustal Omega client via EBI REST API

Submit sequences, poll for completion, parse aligned FASTA.
Replaces alignment.py stub with clustal_client.py.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 5: ConSurf Fallback

**Files:**
- Replace: `varis/m4_conservation/consurf_fallback.py`
- Modify: `tests/test_m4_conservation.py`

**Step 1: Write failing tests**

Add to `tests/test_m4_conservation.py`:

```python
class TestConSurfFallback:
    """Tests for consurf_fallback.py — pre-computed conservation."""

    def test_consurf_known_protein(self, m1_completed_record):
        """Mocked: returns grade, mapped to score."""
        from varis.m4_conservation.consurf_fallback import fetch_consurf
        # ConSurf grades: 1=variable, 9=conserved
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "grades": {str(i): {"grade": 5} for i in range(1, 1864)},
        }
        # Override position 1699 to grade 9 (conserved)
        mock_response.json.return_value["grades"]["1699"] = {"grade": 9}
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        result = fetch_consurf(m1_completed_record, client=mock_client)
        assert result.conservation_score is not None
        assert result.conservation_score > 0.8  # Grade 9 → high score
        assert result.conservation_method == "consurf"

    def test_consurf_unknown_protein(self, m1_completed_record):
        """Mocked: 404 → conservation_available=False."""
        from varis.m4_conservation.consurf_fallback import fetch_consurf
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        result = fetch_consurf(m1_completed_record, client=mock_client)
        assert result.conservation_available is False

    def test_consurf_no_uniprot_id(self, m1_completed_record):
        """No uniprot_id → skip."""
        from varis.m4_conservation.consurf_fallback import fetch_consurf
        m1_completed_record.uniprot_id = None
        result = fetch_consurf(m1_completed_record)
        assert result.conservation_available is False
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m4_conservation.py::TestConSurfFallback -v`
Expected: FAIL

**Step 3: Implement consurf_fallback.py**

Replace `varis/m4_conservation/consurf_fallback.py`. Key points:
- Query ConSurf DB API (best-effort — the actual API endpoint may vary)
- Map ConSurf grade (1-9) to conservation_score: `score = (grade - 1) / 8.0`
  - Grade 1 → 0.0 (variable), Grade 9 → 1.0 (conserved)
- Set `conservation_method = "consurf"`
- Accept optional `client: httpx.Client | None` for DI
- If API unavailable (404, timeout, etc.): set_feature_status with appropriate reason, return

**Step 4: Run tests**

Run: `pytest tests/test_m4_conservation.py::TestConSurfFallback -v`
Expected: PASS

**Step 5: Commit**

```
feat(m4): implement ConSurf fallback for pre-computed conservation

Best-effort ConSurf DB lookup. Maps grade (1-9) to score (0.0-1.0).
Falls back gracefully when unavailable.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 6: M4 Orchestrator with Caching

**Files:**
- Modify: `varis/m4_conservation/__init__.py`
- Modify: `tests/test_m4_conservation.py`

**Step 1: Write failing tests**

Add to `tests/test_m4_conservation.py`:

```python
class TestM4Orchestrator:
    """Tests for M4 orchestration — caching, fallback, integration."""

    def test_m4_no_sequence(self, m1_completed_record):
        """No protein_sequence → M4 fails gracefully."""
        from varis.m4_conservation import run
        m1_completed_record.protein_sequence = None
        m1_completed_record.uniprot_id = None
        result = run(m1_completed_record)
        assert "M4" in result.modules_failed

    def test_m4_fallback_to_consurf(self, m1_completed_record):
        """Primary fails → ConSurf runs."""
        from varis.m4_conservation import run
        # Patch orthologs to fail, consurf to succeed
        with patch("varis.m4_conservation.uniprot_orthologs.fetch_orthologs",
                    return_value=(m1_completed_record, None)):
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.json.return_value = {
                "grades": {"1699": {"grade": 8}},
            }
            mock_client = MagicMock()
            mock_client.get.return_value = mock_response
            with patch("varis.m4_conservation.consurf_fallback.httpx.Client",
                        return_value=mock_client):
                result = run(m1_completed_record)
        assert result.conservation_score is not None
        assert result.conservation_method == "consurf"

    def test_m4_integration(self, m1_completed_record):
        """Full pipeline with mocks — golden record check."""
        from varis.m4_conservation import run
        # Mock orthologs
        orthologs = {
            "sequences": {
                "query": "MKRST" * 340,  # ~1700 chars
                **{f"orth{i}": "MKRST" * 340 for i in range(15)},
            },
            "query_id": "query",
            "taxonomy": {f"orth{i}": 9606 + i for i in range(15)},
        }
        with patch("varis.m4_conservation.uniprot_orthologs.fetch_orthologs",
                    return_value=(m1_completed_record, orthologs)):
            # Mock Clustal to return same sequences (already "aligned")
            with patch("varis.m4_conservation.clustal_client.run_alignment",
                        return_value=(m1_completed_record, orthologs)):
                result = run(m1_completed_record)
        assert "M4" in result.modules_completed
        assert result.conservation_available is True
        assert result.conservation_score is not None
        assert 0.0 <= result.conservation_score <= 1.0
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_m4_conservation.py::TestM4Orchestrator -v`
Expected: FAIL

**Step 3: Update M4 orchestrator**

Replace `varis/m4_conservation/__init__.py`:

```python
"""M4: Conservation Engine — Calculates evolutionary conservation.

CRITICAL INDEPENDENCE: M4 depends on M1 (needs protein sequence), NOT on M2 or M3.
If all structural analysis fails, conservation still works.

Pipeline: UniProt orthologs → Clustal Omega MSA → Shannon entropy scoring.
Fallback: ConSurf pre-computed grades.
Caching: Conservation is protein-level, cached by uniprot_id.

Depends on: M1 (needs uniprot_id, protein_sequence, residue_position)
Populates: conservation_score, conservation_method, num_orthologs,
           position_entropy, conserved_across_mammals, msa_* fields
"""

import json
import logging
from pathlib import Path

from varis.config import DATA_DIR

logger = logging.getLogger(__name__)

_CACHE_DIR = DATA_DIR / "conservation"


def run(variant_record):
    """Execute M4: evolutionary conservation analysis.

    Order: check cache → fetch orthologs → align → score → (consurf fallback).
    Results are cached by uniprot_id for reuse across variants.

    Args:
        variant_record: VariantRecord with M1 fields populated.

    Returns:
        VariantRecord with conservation fields populated (or None with reasons).
    """
    from varis.m4_conservation.clustal_client import run_alignment
    from varis.m4_conservation.conservation_scorer import score_conservation
    from varis.m4_conservation.consurf_fallback import fetch_consurf
    from varis.m4_conservation.uniprot_orthologs import fetch_orthologs

    uniprot_id = variant_record.uniprot_id

    # Check cache first
    if uniprot_id:
        cached = _load_cache(uniprot_id, variant_record)
        if cached:
            if variant_record.conservation_score is not None:
                variant_record.mark_module_completed("M4")
                return variant_record

    # Step 1: Fetch orthologs from UniProt
    orthologs = None
    try:
        variant_record, orthologs = fetch_orthologs(variant_record)
    except Exception as e:
        logger.warning("Ortholog fetch failed: %s", e)
        variant_record.mark_module_failed("M4.orthologs")

    # Step 2: Multiple sequence alignment via Clustal Omega
    alignment = None
    if orthologs is not None:
        try:
            variant_record, alignment = run_alignment(variant_record, orthologs)
        except Exception as e:
            logger.warning("Alignment failed: %s", e)
            variant_record.mark_module_failed("M4.alignment")

    # Step 3: Score conservation from alignment
    if alignment is not None:
        try:
            variant_record = score_conservation(variant_record, alignment)
        except Exception as e:
            logger.warning("Conservation scoring failed: %s", e)
            variant_record.mark_module_failed("M4.scoring")

    # Fallback: ConSurf pre-computed scores
    if variant_record.conservation_score is None:
        try:
            variant_record = fetch_consurf(variant_record)
        except Exception as e:
            logger.warning("ConSurf fallback failed: %s", e)
            variant_record.mark_module_failed("M4.consurf")

    # Save to cache if we got a score
    if uniprot_id and variant_record.conservation_score is not None:
        _save_cache(uniprot_id, variant_record)

    if variant_record.conservation_score is not None:
        variant_record.mark_module_completed("M4")
    else:
        variant_record.mark_module_failed("M4")

    return variant_record


def _load_cache(uniprot_id: str, variant_record) -> bool:
    """Load cached conservation scores for a protein.

    Args:
        uniprot_id: UniProt accession.
        variant_record: VariantRecord to populate from cache.

    Returns:
        True if cache hit, False if miss.
    """
    cache_path = _CACHE_DIR / f"{uniprot_id}_scores.json"
    if not cache_path.exists():
        return False

    try:
        with open(cache_path) as f:
            data = json.load(f)

        position = str(variant_record.residue_position)
        if position not in data.get("positions", {}):
            return False

        pos_data = data["positions"][position]
        variant_record.conservation_score = pos_data.get("conservation_score")
        variant_record.position_entropy = pos_data.get("position_entropy")
        variant_record.msa_column_index = pos_data.get("msa_column_index")
        variant_record.msa_gap_fraction_at_site = pos_data.get("msa_gap_fraction_at_site")
        variant_record.conserved_across_mammals = pos_data.get("conserved_across_mammals")
        variant_record.conservation_method = data.get("method", "clustal_omega")
        variant_record.num_orthologs = data.get("num_orthologs")
        variant_record.msa_num_sequences = data.get("msa_num_sequences")
        variant_record.set_feature_status("conservation", True)

        logger.info("Cache hit for %s position %s", uniprot_id, position)
        return True

    except Exception as e:
        logger.warning("Cache load failed: %s", e)
        return False


def _save_cache(uniprot_id: str, variant_record) -> None:
    """Save conservation scores to cache.

    Args:
        uniprot_id: UniProt accession.
        variant_record: VariantRecord with conservation fields populated.
    """
    try:
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache_path = _CACHE_DIR / f"{uniprot_id}_scores.json"

        # Load existing cache or create new
        data = {}
        if cache_path.exists():
            with open(cache_path) as f:
                data = json.load(f)

        data.setdefault("positions", {})
        data["method"] = variant_record.conservation_method
        data["num_orthologs"] = variant_record.num_orthologs
        data["msa_num_sequences"] = variant_record.msa_num_sequences

        position = str(variant_record.residue_position)
        data["positions"][position] = {
            "conservation_score": variant_record.conservation_score,
            "position_entropy": variant_record.position_entropy,
            "msa_column_index": variant_record.msa_column_index,
            "msa_gap_fraction_at_site": variant_record.msa_gap_fraction_at_site,
            "conserved_across_mammals": variant_record.conserved_across_mammals,
        }

        with open(cache_path, "w") as f:
            json.dump(data, f, indent=2)

        logger.info("Cached conservation for %s position %s", uniprot_id, position)

    except Exception as e:
        logger.warning("Cache save failed: %s", e)
```

**Step 4: Run tests**

Run: `pytest tests/test_m4_conservation.py -v`
Expected: PASS

Run: `pytest tests/ -v --tb=short`
Expected: ALL PASS

**Step 5: Remove old alignment.py stub**

```bash
git rm varis/m4_conservation/alignment.py
```

**Step 6: Commit**

```
feat(m4): implement M4 orchestrator with caching and fallback

Orchestrates: cache check → UniProt orthologs → Clustal Omega →
entropy scoring → ConSurf fallback. Caches by uniprot_id.
Removes alignment.py stub (replaced by clustal_client.py).

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 7: Update Config and Dependencies

**Files:**
- Modify: `varis/config.py`

**Step 1: Add conservation-related URLs to config**

Add to `varis/config.py` API ENDPOINTS section:

```python
CONSURF_API_URL = "https://consurf.tau.ac.il/api"
CLUSTAL_OMEGA_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
```

**Step 2: Run full test suite**

Run: `pytest tests/ -v --tb=short`
Expected: ALL PASS

**Step 3: Commit**

```
chore: add ConSurf and Clustal Omega API URLs to config

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 8: Final Verification and Cleanup

**Step 1: Search for stale references**

Search for `run_blast` imports (should only be in blast_client.py itself and old __init__.py imports that were replaced). Search for `alignment.py` references.

**Step 2: Verify no import cycles**

```
python -c "from varis.m4_conservation import run"
```

**Step 3: Run full test suite**

```
pytest tests/ -v --tb=short
```

**Step 4: Commit any cleanup**

```
chore: Phase 3 cleanup — remove stale M4 references

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```
