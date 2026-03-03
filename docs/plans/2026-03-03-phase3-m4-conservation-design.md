# Phase 3 Design: M4 — Conservation Analysis

**Date:** 2026-03-03
**Status:** Approved
**Depends on:** Phase 1 (M1 ingestion) — complete
**Independent of:** Phase 2 (M2/M3 structural analysis)

## Decisions Made

- **UniProt REST API** for ortholog retrieval (fast, no BLAST polling)
- **EBI Clustal Omega API** for multiple sequence alignment
- **ConSurf DB API** as best-effort fallback
- **Approach A: Compute-first** — full pipeline primary, ConSurf fallback
- **Caching by uniprot_id** — conservation is protein-level, reused across variants
- **BLAST stub kept** — not deleted, available as future alternative

## Module Architecture

### Files

| File | Responsibility |
|------|---------------|
| `uniprot_orthologs.py` | Fetch ortholog sequences from UniProt REST API. Cap 10-100 sequences. Filter taxonomy. Record selection rules. |
| `clustal_client.py` | Submit sequences to EBI Clustal Omega API, poll for MSA result. |
| `conservation_scorer.py` | Build sequence-index → column map (gap-aware). Shannon entropy normalized by log2(20). Mammal conservation with ≥5 threshold. |
| `consurf_fallback.py` | Best-effort ConSurf DB lookup by UniProt ID. Fail gracefully if unavailable. |
| `__init__.py` (update) | Orchestrate with caching: check cache → orthologs → align → score → (consurf fallback). |
| `blast_client.py` (keep) | Stub retained for future BLAST-based ortholog retrieval. |

### Data Flow

```
M1 output (uniprot_id, protein_sequence, residue_position)
  ↓
Cache check: data/conservation/{uniprot_id}_scores.json exists?
  → Yes: load cached scores, extract at position, done
  → No: continue pipeline
  ↓
uniprot_orthologs: fetch ortholog sequences + taxonomy from UniProt
  - Min 10, max 100 sequences
  - Filter by identity (30-90% to capture useful divergence)
  - Record taxonomy IDs for mammal classification
  - If <10 orthologs: conservation_available=False, reason=insufficient_orthologs
  ↓
clustal_client: submit query + orthologs to EBI Clustal Omega
  - Poll for completion (~30s typical)
  - Parse FASTA alignment result
  ↓
conservation_scorer: extract column at mutation position
  - Build sequence-index → alignment-column map (count gaps)
  - Compute Shannon entropy at column
  - Compute mammal conservation percentage
  - Cache all per-position scores
  ↓
(if any step fails) consurf_fallback: query ConSurf DB
  - Pre-computed conservation grade (1-9)
  - Map grade to conservation_score approximation
```

## Fields Written

```
# Primary conservation outputs
conservation_score             — float: 1 - (H / log2(20)), range [0.0, 1.0]
                                  1.0 = perfectly conserved, 0.0 = maximally variable
conservation_method            — str: "clustal_omega" or "consurf"
num_orthologs                  — int: number of ortholog sequences used
position_entropy               — float: raw Shannon entropy at position (bits)
conserved_across_mammals       — bool: ref AA in ≥90% mammalian orthologs

# MSA metadata (debugging + ML features)
msa_num_sequences              — int: total sequences in alignment (including query)
msa_gap_fraction_at_site       — float: fraction of sequences with gap at mutation column
msa_column_index               — int: alignment column index for mutation position

# Feature availability (already in schema)
conservation_available         — bool
conservation_missing_reason    — str
```

## Key Definitions

### Shannon Entropy

```
H = -Σ p(aa) * log2(p(aa))
```

- Computed over 20 standard amino acids only
- Gaps **excluded** from frequency calculation (but gap fraction stored separately)
- `msa_gap_fraction_at_site` recorded for transparency
- Normalization: `conservation_score = 1 - (H / log2(20))`
- log2(20) ≈ 4.322 is the maximum possible entropy (all 20 AAs equally frequent)

### Mammal Conservation

- **Mammal**: NCBI taxonomy class Mammalia (taxon ID 40674)
- **Threshold**: ref AA present in ≥90% of mammalian ortholog sequences
- **Minimum**: Require ≥5 mammalian sequences; if fewer, `conserved_across_mammals=None`
- **Gap handling**: Ignore sequences with gaps at the mutation position when counting

### Ortholog Selection

- **Source**: UniProt REST API, cross-references or sequence clusters
- **Identity range**: 30-90% to the query (too similar = uninformative, too distant = noise)
- **Count**: min 10, max 100 sequences
- **Below minimum**: `conservation_available=False`, reason `insufficient_orthologs`
- **Record**: `ortholog_selection_rule` and `taxonomy_filter_applied` in cache metadata

### MSA Position Mapping

The mutation position in the protein sequence does NOT directly map to the MSA column
due to gaps. The scorer must:

1. Find the query sequence row in the MSA
2. Walk the query row, counting non-gap characters
3. When count == residue_position, that's the alignment column
4. Verify: the character at that position should match ref_amino_acid
5. Store `msa_column_index` for debugging

## Caching Strategy

Conservation is protein-level — every variant in the same protein reuses the same
alignment and per-position scores. Only the column extraction runs per variant.

```
data/conservation/{uniprot_id}_orthologs.json    — ortholog sequences + metadata
data/conservation/{uniprot_id}_msa.fasta          — Clustal Omega alignment
data/conservation/{uniprot_id}_scores.json        — per-position: entropy, gap_fraction, etc.
```

Cache is keyed by `uniprot_id`. Cache invalidation: manual or by schema version bump.
No TTL in Phase 3 — protein orthologs don't change frequently.

## Schema Changes

Add to VariantRecord (bump to v1.3.0):

```python
# New MSA metadata fields
msa_num_sequences: Optional[int] = None
msa_gap_fraction_at_site: Optional[float] = None
msa_column_index: Optional[int] = None
```

`conservation_score`, `conservation_method`, `num_orthologs`, `position_entropy`,
`conserved_across_mammals`, `conservation_available`, `conservation_missing_reason`
already exist in v1.2.0.

Also add `NullReason.INSUFFICIENT_DATA = "insufficient_data"` for the <10 orthologs case.

## Testing Strategy

### Test Markers

- Default: offline, mocked API responses
- `@pytest.mark.integration`: real API calls (UniProt, Clustal Omega, ConSurf)

### Tests (`tests/test_m4_conservation.py`)

| Test | Assertion Style |
|------|----------------|
| **UniProt Orthologs** | |
| `test_fetch_orthologs_brca1` | Mocked: returns ≥10 sequences with taxonomy (mocked) |
| `test_fetch_orthologs_caps` | Max 100 sequences returned even if more available |
| `test_fetch_orthologs_too_few` | <10 sequences → returns None, reason recorded |
| `test_fetch_orthologs_no_uniprot` | No uniprot_id → skip with upstream dependency reason |
| **Clustal Omega** | |
| `test_clustal_alignment` | Mocked: submit + poll returns valid FASTA alignment |
| `test_clustal_timeout` | Mocked: poll exceeds max retries → TIMED_OUT |
| `test_clustal_no_sequences` | No sequences → skip |
| **Conservation Scorer** | |
| `test_entropy_fully_conserved` | All same AA → entropy=0, score=1.0 |
| `test_entropy_maximally_variable` | All different AAs → entropy≈log2(20), score≈0.0 |
| `test_entropy_ignores_gaps` | Gaps excluded from frequency, gap_fraction stored |
| `test_position_mapping_with_gaps` | Known MSA with gaps → correct column extracted |
| `test_position_mapping_validates_ref_aa` | Mismatch detected and logged |
| `test_mammal_conservation_threshold` | ≥90% mammals have ref AA → True |
| `test_mammal_conservation_too_few` | <5 mammals → conserved_across_mammals=None |
| **ConSurf Fallback** | |
| `test_consurf_known_protein` | Mocked: returns grade, mapped to score |
| `test_consurf_unknown_protein` | Mocked: 404 → conservation_available=False |
| **Orchestrator** | |
| `test_m4_caching` | Second call uses cached result |
| `test_m4_fallback_to_consurf` | Primary fails → ConSurf runs |
| `test_m4_no_sequence` | No protein_sequence → M4 fails gracefully |
| `test_m4_integration` | Full pipeline with mocks, golden record check |

## Dependencies

No new Python packages needed — uses httpx (already installed) for all API calls.
BioPython's `AlignIO` for parsing FASTA alignments (already installed).
