# Phase 2 Design: M2 (Structure) + M3 (Structural Analysis)

**Date:** 2026-03-03
**Status:** Approved
**Depends on:** Phase 1 (M1 ingestion) — complete

## Decisions Made

- **Free tools only** — FoldX/PyRosetta wrappers as optional stubs that activate when binary detected
- **ESM Metagenomic Atlas API** for structure fallback (sequence ≤400aa only)
- **BioPython DSSP** with `mkdssp` binary for secondary structure
- **InterPro REST API** for domain identification (replaces HMMER, no local Pfam DB)
- **Approach A: Thin wrappers, flat orchestration** — matches M1 pattern

## M2: Structure Retrieval & Preparation

### Purpose

Validate the structure from M1, optionally retrieve a fallback structure, conditionally
repair it, and extract quality metrics. Prepare the structure for M3 analysis.

### Files

| File | Responsibility |
|------|---------------|
| `structure_validator.py` | Validate PDB exists, check mutation residue is present, extract pLDDT (AlphaFold/ESMFold only), detect missing residues near site |
| `esmfold_predictor.py` | Optional fallback: ESM Atlas API. Skip with `sequence_too_long` if >400aa |
| `pdb_fixer.py` | Conditional: only runs when validation detects missing heavy atoms. Adds atoms/hydrogens. Never reconstructs missing residues |
| `run.py` | Orchestrates: validate → (esmfold if no structure) → (fix if needed) → extract quality |

### Fields Written

```
pdb_path                        — str: path to usable structure (from M1 or ESMFold)
pdb_fixed_path                  — str: path to repaired PDB (only if fixer ran)
pdb_hash                        — str: SHA256 of PDB for provenance
structure_source_url             — str: URL structure was fetched from

mutation_site_present            — bool: is the residue in the structure?
mutation_site_plddt              — float: pLDDT at mutation residue (null if not predicted)
plddt_mean                       — float: mean pLDDT across protein (null if not predicted)
plddt_available                  — bool: whether pLDDT is meaningful for this source

mutation_site_confidence_bucket  — "high" / "medium" / "low" (based on pLDDT thresholds)
numbering_scheme                 — "uniprot_canonical" (explicit, for future PDB support)

structure_quality_summary        — dict: {plddt_mean, plddt_site, percent_low_confidence}
preparation_steps                — list[str]: e.g. ["validated", "added_hydrogens"]
```

### Key Logic

- pLDDT extraction gated: `if structure_source in ("alphafold", "esmfold")`
- If `mutation_site_present=False` → `coordinate_mapping_confidence="failed"`, M3 skips site-dependent features
- If pLDDT at site < 70 → `mutation_site_confidence_bucket="low"`, M3 runs but results carry lower weight in M5
- PDBFixer only runs if validation detects missing heavy atoms — never reconstructs missing residues/loops
- PDBFixer does NOT "fill gaps" — missing residues near mutation site → `mutation_site_present=False` with reason `missing_residues_near_site`

### Design Rationale

- **No gap filling:** Reconstructing missing loops is modeling, not fixing. Creates false certainty, violates project honesty principle.
- **pLDDT source-awareness:** B-factor column means pLDDT only in AlphaFold/ESMFold. Experimental PDB B-factors are thermal displacement — different quantity entirely.
- **ESMFold 400aa limit:** BRCA1 is 1863aa — fallback won't work for many clinically relevant proteins. Handled with explicit `sequence_too_long` reason code.
- **Explicit numbering scheme:** AlphaFold uses UniProt canonical numbering. Making this explicit prevents bugs when PDB structures are added later.
- **Conditional PDBFixer:** AlphaFold PDBs are clean enough for SASA/DSSP. Running fixer by default adds latency for no benefit.
- **`structure_resolution` dropped:** Predicted models don't have experimental resolution. Quality is represented by pLDDT-based metrics instead.

## M3: Structural Analysis

### Purpose

Extract structural features from the prepared 3D structure. Each feature is independently
computed — if one fails, others still run.

### Files

| File | Priority | Responsibility | External Dependency |
|------|----------|---------------|-------------------|
| `freesasa_wrapper.py` | 1 | Relative SASA at mutation site, burial classification (core/surface) | `freesasa` library |
| `dssp_wrapper.py` | 2 | Secondary structure assignment (helix/sheet/coil) via BioPython DSSP | `mkdssp` binary |
| `biopython_contacts.py` | 3 | WT local environment: heavy-atom contacts within 4.5Å, H-bonds, packing density | BioPython |
| `interpro_client.py` | 4 | Pfam domain at mutation position via InterPro REST API | httpx |
| `foldx_wrapper.py` | 5 | Stub: detects FoldX binary, runs ΔΔG if available | FoldX (optional) |
| `pyrosetta_wrapper.py` | 6 | Stub: same pattern as FoldX, fallback ΔΔG | PyRosetta (optional) |
| `run.py` | — | Orchestrates all wrappers. Skips site-dependent tools if `mutation_site_present=False` | — |

### Fields Written

```
# FreeSASA (Priority 1)
solvent_accessibility_relative   — float: relative SASA in [0.0, 1.0] normalized by residue-type max
burial_category                  — "core" / "surface"
sasa_available                   — bool

# DSSP (Priority 2)
secondary_structure              — str: DSSP code ("H", "E", "C", etc.)
secondary_structure_name         — "helix" / "sheet" / "coil"
dssp_available                   — bool

# Contacts (Priority 3) — WT environment only
contacts_wt                      — int: heavy-atom contacts within 4.5Å
hbonds_wt                        — int: H-bonds at position in WT
packing_density                  — float: local packing metric
contacts_available               — bool

# InterPro/Pfam (Priority 4)
domain_name                      — str: Pfam domain name (e.g., "BRCT")
domain_id                        — str: Pfam ID (e.g., "PF00533")
domain_start                     — int: domain boundary start position
domain_end                       — int: domain boundary end position
domain_criticality               — "critical" / "important" / "peripheral"
domain_available                 — bool

# FoldX/PyRosetta (Priority 5-6, optional stubs)
ddg_foldx                        — float: null until licensed
ddg_pyrosetta                    — float: null until licensed
ddg_mean                         — float: null until licensed
ddg_available                    — bool (False in Phase 2)
```

### Key Design Decisions

- **No "interface" burial category:** AlphaFold structures are monomers — no binding partners means no true interface detection. Reserved for future phases with complex prediction.
- **No `helix_disruption`:** Dropped to avoid overclaiming. M5 can learn this pattern from raw features (secondary_structure + amino acid substitution) if it matters.
- **WT environment only:** Without a mutant structure, we cannot measure "changed contacts" or "lost H-bonds." Phase 2 measures the wild-type local environment only. Mutant modeling deferred until FoldX/PyRosetta available.
- **InterPro stores domain boundaries:** `domain_start` and `domain_end` stored for transparency. If mapping confidence isn't high, `domain_available=False` with reason `mapping_uncertain`.
- **Site-dependent gating:** If `mutation_site_present=False`, FreeSASA/DSSP/Contacts are all skipped. InterPro still runs (uses sequence position, not structure).
- **SASA definition:** Relative SASA (%) normalized by residue-type maximum, range [0.0, 1.0]. More comparable across residue types than absolute Ų.
- **Contact definition:** Heavy atoms only within 4.5Å. Hydrogen positions unreliable in predicted structures.

### Fallback Chain (Structure Retrieval)

```
AlphaFold DB (M1 already downloads)
  → ESMFold API (if no AlphaFold, sequence ≤400aa)
    → Skip structure (set all structural fields to null, clear reason codes)
```

## Testing Strategy

### Test Markers

- Default test suite: fully offline, deterministic
- `@pytest.mark.integration`: real API calls (ESMFold, InterPro live), excluded by default

### M2 Tests (`tests/test_m2_structure.py`)

| Test | Assertion Style |
|------|----------------|
| `test_structure_sanity_fixture` | PDB has ≥1699 residues, numbering starts at 1, continuous |
| `test_validate_existing_pdb` | Finds residue 1699, `mutation_site_present=True` |
| `test_validate_missing_residue` | `mutation_site_present=False`, reason code set |
| `test_site_out_of_range` | Position > length → `mutation_site_present=False`, reason `site_out_of_range` |
| `test_plddt_only_for_predicted` | Skipped when source is not alphafold/esmfold |
| `test_plddt_range_valid` | All values in [0, 100] |
| `test_plddt_site_bfactor_value` | Extracted value matches known B-factor at one residue |
| `test_esmfold_sequence_too_long` | >400aa → skip, reason `sequence_too_long` (mocked) |
| `test_esmfold_short_sequence` | Mocked API success for short sequence |
| `test_pdb_fixer_adds_hydrogens` | `preparation_steps` contains `"added_hydrogens"`, output parses, site preserved |
| `test_pdb_fixer_skipped_when_clean` | Clean PDB → `preparation_steps=["validated"]` |
| `test_m2_integration` | Full pipeline on BRCA1, golden record key check |

### M3 Tests (`tests/test_m3_structural.py`)

| Test | Assertion Style |
|------|----------------|
| `test_freesasa_computes_sasa` | Returns float, `sasa_available=True` |
| `test_freesasa_relative_bounds` | Value in [0.0, 1.0] |
| `test_freesasa_burial_category` | Returns "core" or "surface" only |
| `test_freesasa_no_structure` | `sasa_available=False` with reason |
| `test_dssp_returns_valid_code` | Result in allowed DSSP code set |
| `test_dssp_mkdssp_missing` | `dssp_available=False`, reason `mkdssp_missing` |
| `test_contacts_wt_valid` | Integers ≥ 0, nonzero for folded region |
| `test_contacts_heavy_atoms_only` | Confirms heavy-atom-only counting |
| `test_interpro_brca1_domain` | Returns domain with `domain_start`/`domain_end` (mocked) |
| `test_interpro_position_outside_domain` | `domain_available=False` (mocked) |
| `test_foldx_stub_no_binary` | `ddg_available=False`, reason `tool_missing` |
| `test_m3_skips_when_site_absent` | All site-dependent features skipped, InterPro still runs |
| `test_m3_integration_golden_record` | Full pipeline, snapshot keys + reason codes |

## Dependencies to Add

```toml
# In pyproject.toml [project.optional-dependencies.structure]
# Already present: freesasa, mdanalysis, gemmi, openbabel-wheel
# Need to add:
"pdbfixer>=1.9",       # Structure repair (requires OpenMM)
"mdtraj>=1.9",         # Backup DSSP if needed later
```

`mkdssp` binary installed separately via `brew install dssp` or `conda install -c conda-forge dssp`.

## Schema Changes

The VariantRecord schema (`varis/models/variant_record.py`) will need new fields added for:
- `mutation_site_present`, `mutation_site_plddt`, `plddt_available`, `mutation_site_confidence_bucket`
- `numbering_scheme`, `structure_quality_summary`, `preparation_steps`
- `pdb_hash`, `structure_source_url`
- `solvent_accessibility_relative` (replaces `solvent_accessibility`)
- `burial_category` restricted to `"core"` / `"surface"` (no `"interface"`)
- `contacts_wt`, `hbonds_wt`, `packing_density`
- `domain_start`, `domain_end`

Remove from schema (deferred):
- `helix_disruption`
- `hbonds_lost`, `contacts_changed` (require mutant structure)
- `structure_resolution` (replaced by `structure_quality_summary`)

Schema version bump: 1.1.0 → 1.2.0
