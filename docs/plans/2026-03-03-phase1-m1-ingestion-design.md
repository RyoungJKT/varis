# Phase 1 Design: M1 Ingestion

**Date:** 2026-03-03
**Status:** Approved
**Scope:** Implement all 7 M1 sub-modules so `python -m varis BRCA1 p.Arg1699Trp` produces a complete VariantRecord JSON.

---

## Data Flow

```
Input: gene="BRCA1", hgvs="p.Arg1699Trp"
  |
  1. hgvs_parser          -> position, ref/alt AA, charge change
  2. uniprot_client        -> uniprot_id, sequence, protein metadata
  3. variant_normalizer    -> validate ref AA against UniProt sequence, set confidence
  4. clinvar_client        -> clinvar_id, classification, genomic coordinates
  5. gnomad_client         -> frequency (only if ClinVar provided genomic coords + build matches)
  6. alphafold_client      -> download PDB file using uniprot_id
  7. alphamissense_client  -> pathogenicity score via hegelab.org web lookup
  |
Output: VariantRecord JSON -> {gene}_{variant}.json
```

Each step is wrapped in try/except by the existing M1 orchestrator (`m1_ingestion/__init__.py`). If any step fails, its fields stay None with a reason code, and the pipeline continues.

## Key Design Decisions

### 1. ClinVar as genomic coordinate source (not Ensembl variant_recoder)

Protein HGVS (p.Arg1699Trp) does not uniquely define genomic coordinates due to codon degeneracy. Instead of guessing nucleotide changes, we extract genomic coordinates from ClinVar records, which contain exact chromosome/position/ref/alt.

- If variant is in ClinVar -> we have unambiguous genomic coords -> gnomAD query is safe
- If variant is NOT in ClinVar -> gnomAD fields set to null with reason "no_genomic_coordinates"
- This is an acceptable Phase 1 limitation (novel variants skip gnomAD)

### 2. UniProt before normalizer (dependency ordering)

The normalizer validates ref AA against the UniProt protein sequence. UniProt must be fetched first. This prevents stamping "high confidence" before we have data to validate against.

### 3. Simple direct position mapping (Phase 1)

Assume UniProt canonical numbering matches HGVS protein position directly. Validate with ref AA check. Full MANE Select / SIFTS / BLAST alignment deferred to Phase 2+.

### 4. ClinVar disambiguation rule

When ClinVar search returns multiple records for the same protein change:
1. Prefer exact protein HGVS match (same ref/alt and position)
2. Among exact matches, prefer highest review stars
3. If still ambiguous, take first match but set coordinate_mapping_confidence = "low" and skip gnomAD genomic lookup

### 5. Assembly-aware gnomAD queries

gnomAD v4 uses GRCh38, v2 uses GRCh37. The gnomad_client only queries if the ClinVar genomic build matches the target gnomAD endpoint. Mismatched builds -> null with reason.

---

## New VariantRecord Fields

Added to `varis/models/variant_record.py`:

```python
# Genomic coordinates (populated by clinvar_client)
reference_build: str | None = None       # "GRCh37" or "GRCh38"
clinvar_chrom: str | None = None         # e.g. "17"
clinvar_pos: int | None = None           # genomic position
clinvar_ref: str | None = None           # reference nucleotide allele
clinvar_alt: str | None = None           # alternate nucleotide allele
coordinate_source: str | None = None     # "clinvar" for Phase 1
```

Also add to schema/variant_record_schema.json.

---

## Sub-Module Specifications

### 1. hgvs_parser.py

Pure string parsing, no API calls.

- Two regex patterns (already defined in stub): three-letter and single-letter HGVS protein notation
- Extracts: residue_position, ref_amino_acid, alt_amino_acid, ref_aa_single, alt_aa_single
- Calculates charge_change using config.AA_CHARGE mapping
- Invalid input -> set fields to None with reason "validation_failed"

### 2. uniprot_client.py

REST API via httpx.

- Search: `GET https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606&fields=accession,protein_name,sequence,gene_names`
- Take first reviewed (Swiss-Prot) result over unreviewed (TrEMBL)
- Populate: uniprot_id, protein_name, protein_sequence, protein_length, protein_function
- Accept optional httpx.Client parameter for dependency injection in tests

### 3. variant_normalizer.py

Validation against UniProt sequence. No API calls.

- Read protein_sequence from record (set by uniprot_client)
- Check if amino acid at residue_position matches ref_aa_single
- Match -> coordinate_mapping_confidence = "high", uniprot_residue_position = residue_position
- Mismatch -> coordinate_mapping_confidence = "failed", add warning to normalization_warnings
- If no protein_sequence available -> coordinate_mapping_confidence = None, reason "upstream_dependency_failed"

### 4. clinvar_client.py

NCBI E-utilities via httpx.

- Search: `esearch.fcgi?db=clinvar&term={gene}[gene]+AND+{hgvs}[variant+name]`
- Fetch: `efetch.fcgi?db=clinvar&id={id}&rettype=vcv&retmode=xml`
- Parse XML (xml.etree.ElementTree) for: variation ID, classification, review status, conditions
- Extract genomic coordinates from `<SequenceLocation>` element, prefer GRCh38
- Populate: clinvar_id, clinvar_classification, clinvar_review_status, reference_build, clinvar_chrom, clinvar_pos, clinvar_ref, clinvar_alt, coordinate_source
- Disambiguation: exact protein match -> highest review stars -> first match with low confidence
- Rate limit: respect 3 req/sec without API key, 10 with NCBI_API_KEY

### 5. gnomad_client.py

GraphQL API via httpx.

- Gate: only run if clinvar_chrom, clinvar_pos, clinvar_ref, clinvar_alt are all populated AND reference_build matches gnomAD endpoint (v4 = GRCh38)
- Variant ID format: `{chrom}-{pos}-{ref}-{alt}`
- GraphQL POST to https://gnomad.broadinstitute.org/api
- Extract: allele_frequency, popmax, homozygote_count
- No genomic coords -> null with reason "no_genomic_coordinates"
- Build mismatch -> null with reason "validation_failed" (assembly mismatch)

### 6. alphafold_client.py

Direct file download via httpx.

- URL: `https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb`
- Save to: data/structures/AF-{uniprot_id}-F1-model_v4.pdb
- Cache: skip download if file already exists at that path
- Set structure_source = "alphafold", pdb_path = <path>
- No uniprot_id -> null with reason "upstream_dependency_failed"

### 7. alphamissense_client.py

Web lookup via httpx.

- Query hegelab.org API with UniProt ID + protein variant notation
- Extract: am_pathogenicity score, am_class (likely_pathogenic/ambiguous/likely_benign)
- Populate: alphamissense_score, alphamissense_class
- Lookup failure -> null with reason "no_data_available"

---

## M1 Orchestrator Change

Update execution order in `m1_ingestion/__init__.py`:

```
Current:  hgvs_parser -> normalizer -> clinvar -> gnomad -> uniprot -> alphafold -> alphamissense
Revised:  hgvs_parser -> uniprot -> normalizer -> clinvar -> gnomad -> alphafold -> alphamissense
```

---

## Error Handling

Consistent pattern across all sub-modules (per CLAUDE.md):

```python
try:
    # do work
    variant_record.field = result
except Exception as e:
    logger.warning(f"SubModule failed for {variant_record.variant_id}: {e}")
    variant_record.set_with_reason("field", None, "tool_crashed")
return variant_record
```

---

## Testing

Integration tests with real API calls against BRCA1 p.Arg1699Trp.

- Each sub-module gets its own test class in test_m1_ingestion.py
- Assertions validate against known values from data/validation/BRCA1_p.Arg1699Trp_expected.json:
  - ClinVar ID = VCV000055361
  - AlphaMissense score ~ 0.934
  - UniProt ID = P38398
  - Structure source = alphafold
- @pytest.mark.timeout(30) on API-calling tests
- Pipeline integration test: full flow produces valid VariantRecord JSON with expected fields populated

---

## Phase 1 Limitations (explicit, by design)

1. Variants not in ClinVar will have no gnomAD frequency data
2. Position mapping is direct (no MANE Select / SIFTS / BLAST alignment)
3. AlphaMissense uses web lookup (not local 4GB TSV)
4. No genomic coordinate inference from protein HGVS (deferred to Phase 2+)
5. No liftover between GRCh37 and GRCh38
