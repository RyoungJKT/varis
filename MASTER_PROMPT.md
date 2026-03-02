# MASTER_PROMPT.md — Domain Knowledge & Scientific Context for Varis

> **Author and sole developer:** Russell Young
> **Organization:** Russell Genetics (founded by Russell Young)
> **Date:** March 2026

This document provides the scientific domain knowledge, strategic reasoning, and
implementation guidance that Claude Code needs to build Varis correctly. Read this
alongside CLAUDE.md (which covers coding standards and architecture rules).

CLAUDE.md tells you HOW to code. This document tells you WHAT you are building and WHY.

---

## 1. Who Is Building This

Russell Young is a Year 12 IB Diploma student at British School Jakarta. He is the
sole founder, developer, architect, and maintainer of Russell Genetics and all its
products. There is no team. There are no co-founders. Every line of code, every
design decision, every document, and every scientific choice is Russell's.

When writing comments, documentation, commit messages, or any text that references
the creator of this project, use "Russell Young" or "Russell." Never invent
co-authors, contributors, or team members. The README, pyproject.toml, LICENSE,
and any publication-facing text should credit Russell Young as sole author.

---

## 2. What Varis Is (The Scientific Mission)

### The Problem: The VUS Crisis

When a patient receives a genetic test, roughly 40% of results come back as a
**Variant of Uncertain Significance (VUS)** — a mutation that might cause disease,
but nobody can confirm it. There are ~2.2 million VUS in ClinVar (NIH's central
variant database), and ~400,000 new ones are added every year.

For each VUS, a patient may wait 2–5+ years for reclassification. Functional
studies to resolve a single VUS cost $5,000–$50,000. Roughly 300 million people
globally have a rare disease, and many are stuck in "diagnostic limbo" because
their variant is a VUS.

### What Varis Does

Varis is a **structural investigator** for genetic variants. It determines whether
a patient's missense mutation causes disease by investigating the structural
evidence — not by classifying with a black-box score.

A **missense mutation** changes one amino acid in a protein. Whether this causes
disease depends on where that amino acid sits in the protein's 3D structure:
- A change **buried in the hydrophobic core** is devastating (destabilizes the fold)
- A change **on a flexible surface loop** may be tolerated
- A change **near the active site** disrupts function
- A change **in a critical domain** (like BRCT in BRCA1) is likely damaging

Varis exploits the AlphaFold revolution — free, accurate 3D structure predictions
for ~20,000 human proteins — to investigate these questions computationally.

### What Varis Is NOT

- Varis does NOT discover which disease a patient has
- Varis does NOT replace genetic counselors or clinicians
- Varis does NOT do experimental/wet-lab structural biology
- Varis provides **structural evidence to confirm a diagnosis doctors already suspect**

### The Clinical Workflow

1. Patient presents with symptoms → Doctor suspects a genetic disease
2. Genetic test reveals a missense variant in a relevant gene
3. Lab classifies it as VUS (uncertain)
4. **Varis investigates**: analyzes the variant's effect on protein structure
5. Varis produces structural evidence: ΔΔG, burial depth, domain context, ACMG codes
6. Genetic counselor uses this evidence to reclassify the VUS
7. Patient gets a diagnosis and treatment begins

---

## 3. The Key Differentiator: Structural Investigator vs. Classifier

### AlphaMissense (Google DeepMind, 2023)

AlphaMissense pre-computed pathogenicity predictions for all 71 million possible
human missense variants. It gives clinicians **a single score (0–1)**. No explanation.
No structural mechanism. No clinical report. 11% of variants are left as "ambiguous."

The model is a deep learning transformer trained on population frequency data
(weak labels). The architecture is opaque — no one can explain why a specific
variant got its score. The model is **frozen at 2023** and will not be updated
(DeepMind's GitHub explicitly states the repository will not be actively maintained).

### Varis (Russell Genetics, 2026)

Varis investigates WHY a mutation is damaging. For each variant, it produces:
- **ΔΔG calculation**: How much the protein is destabilized (in kcal/mol)
- **3D visualization**: Where the damage sits in the protein's architecture
- **SHAP explanation**: Which features drove the ML prediction
- **ACMG evidence codes**: PM1, PM5, PP3, PP2, PS3-proxy
- **Clinical report**: PDF with all evidence for the genetic counselor

The ML model is an **interpretable ensemble** (CatBoost + XGBoost + LightGBM) trained
on expert-reviewed ClinVar classifications (strong labels), with full SHAP explainability.

### The Relationship Is Complementary, Not Competitive

Varis **uses** AlphaMissense as one of ~15 input features in its ML ensemble.
AlphaMissense pre-screens; Varis investigates why. Varis incorporates the
AlphaMissense score alongside its own structural features. The tools are
complementary — AlphaMissense is fast and broad, Varis is deep and explanatory.

When implementing the AlphaMissense client (M1), always treat it as "one feature
among many," never as a competing tool to be excluded.

---

## 4. Domain Knowledge: Protein Structure

### What You Need to Know About Protein Structure

**Amino acids** are the building blocks of proteins. There are 20 standard amino acids,
each with different properties (size, charge, hydrophobicity). A missense mutation
swaps one amino acid for another.

**Protein folding**: The amino acid chain folds into a 3D structure. The folded
structure determines function. The interior ("hydrophobic core") is tightly packed
with nonpolar residues. The surface is exposed to water.

**Secondary structure** (assigned by DSSP):
- **H** = alpha-helix (rigid, rod-like)
- **E** = beta-sheet (flat, extended)
- **C/T/S** = coil/turn/bend (flexible)
- A mutation that breaks a helix or sheet is structurally catastrophic

**ΔΔG (delta-delta-G)**: The change in folding free energy caused by a mutation.
- ΔΔG ≈ 0 kcal/mol → neutral (protein stability unchanged)
- ΔΔG > 2 kcal/mol → destabilizing (protein is significantly less stable)
- ΔΔG > 4 kcal/mol → severely destabilizing (protein likely unfolds)
- This is the single most powerful structural feature for pathogenicity prediction

**Solvent Accessible Surface Area (SASA)**: How exposed a residue is to water.
- < 10% exposed → buried in core (mutations here are devastating)
- 10–40% → partially buried/interface
- > 40% → surface-exposed (mutations may be tolerated)

**pLDDT**: AlphaFold's per-residue confidence score (0–100).
- > 90 → very high confidence in predicted structure
- 70–90 → confident
- 50–70 → low confidence
- < 50 → disordered region (structural predictions unreliable)
- When pLDDT < 70 at the mutation site, flag that structural features may be unreliable

**Conservation**: If an amino acid is the same across many species (mammals, vertebrates),
evolution has preserved it because mutations there are harmful. 100% conservation
across mammals strongly suggests the position is structurally or functionally critical.

### ACMG Evidence Framework

The **American College of Medical Genetics (ACMG)** framework is how clinicians
classify variants. Evidence codes are combined to reach a classification:

| Code | Evidence Type | How Varis Provides It |
|------|--------------|----------------------|
| PP3 | Computational prediction supports pathogenic | Ensemble score > 0.8 |
| PM1 | Located in critical functional domain | HMMER domain identification |
| PM5 | Known pathogenic variant at same position | ClinVar cross-reference |
| PS3-proxy | Functional damage evidence | ΔΔG > 2 kcal/mol |
| PP2 | Low rate of benign missense in this gene | gnomAD analysis |

AlphaMissense provides only PP3. Varis provides PP3 + PM1 + PM5 + PS3-proxy + PP2.
This is a core differentiator — always implement ACMG mapping carefully.

**Critical framing rule**: Varis SUGGESTS evidence codes to support professional
review. It does NOT replace ACMG adjudication by qualified variant curation teams.
All evidence codes — particularly PS3-proxy, which is computational structural
evidence rather than wet-lab functional data — must be clearly labelled as
computational suggestions. Every report, API response, and user-facing output that
includes ACMG codes MUST include a disclaimer to this effect. Never present Varis's
codes as final clinical classifications.

---

## 5. Data Sources and API Guidance

All data comes from publicly funded, open-access databases. Varis requires
**zero access to patients, hospitals, or medical records.**

### ClinVar (NIH) — via NCBI E-utilities
- ~2.2M VUS + pathogenic/benign training labels
- API: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`
- Use `esearch.fcgi` to find variants, `efetch.fcgi` to retrieve records
- Rate limit: 3 requests/sec without API key, 10 with key (set NCBI_API_KEY)
- Return format: XML (parse with xml.etree.ElementTree)
- Search by gene + protein change: `{gene}[gene] AND {hgvs}[variant name]`
- ClinVar IDs look like: VCV000055361 (variation), RCV000074440 (record)

### gnomAD (Broad Institute) — via GraphQL API
- Allele frequencies across populations
- API: `https://gnomad.broadinstitute.org/api`
- Uses GraphQL — send POST with query string
- Key fields: `allele_frequency`, `popmax`, `homozygote_count`
- A variant with frequency < 0.0001 is rare (consistent with pathogenic)
- Absence from gnomAD is also informative (very rare or de novo)

### UniProt — via REST API
- Protein sequences, functional annotations, domain boundaries
- API: `https://rest.uniprot.org/uniprotkb/search`
- Search: `gene_exact:{gene} AND organism_id:9606` (9606 = human)
- Returns JSON with sequence, function description, domains, GO terms
- The `features` array contains domain annotations with start/end positions

### AlphaFold DB — via REST API
- Predicted 3D structures for ~20,000 human proteins
- API: `https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}`
- Downloads PDB file from: `https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb`
- pLDDT scores are stored in the B-factor column of the PDB file
- If no AlphaFold structure exists, fall back to ESMFold

### ESMFold (Meta) — via ESM Atlas API
- Fallback structure prediction when AlphaFold DB doesn't have the protein
- API: `https://api.esmatlas.com/foldSequence/v1/pdb/`
- POST with raw amino acid sequence, returns PDB
- Slower and sometimes less accurate than AlphaFold, but covers proteins AlphaFold doesn't
- For local GPU prediction: `pip install fair-esm` (requires GPU)

### AlphaMissense — via pre-computed lookup
- Pre-computed scores for all 71M human missense variants
- Downloaded from Zenodo as a large TSV file (~4GB)
- Columns: `uniprot_id`, `protein_variant`, `am_pathogenicity`, `am_class`
- For individual lookups, use the community web resource at alphamissense.hegelab.org
- Or load the TSV into a local SQLite database for batch queries

### BLAST — via NCBI web API
- Find orthologous sequences across species for conservation analysis
- API: `https://blast.ncbi.nlm.nih.gov/Blast.cgi`
- Submit: PUT request with `PROGRAM=blastp`, `DATABASE=swissprot`, `QUERY={sequence}`
- Poll for results using the RID (Request ID)
- Parse results to extract top ~50 ortholog sequences
- This is the one API that requires polling/waiting (BLAST runs take 30–120 seconds)

---

### Variant Normalization and Coordinate Mapping (Critical M1 Sub-Module)

Variant normalization is one of the hardest and most error-prone steps in the
pipeline. The same mutation can have different notations across databases, and
residue numbering varies between transcripts, isoforms, and structures. **A wrong
position means silently wrong ΔΔG, wrong SASA, wrong conservation score —
everything downstream is garbage.**

This is implemented as `variant_normalizer.py` under M1. It runs AFTER hgvs_parser
and BEFORE any database lookups.

**Common traps** (build tests for all of these):

1. **Transcript differences**: BRCA1 has NM_007294.4 and NM_007297.4 with different
   exon structures. The same genomic position maps to different protein positions.
   Solution: always resolve to the MANE Select transcript first.

2. **Multiple notations**: p.Arg1699Trp = p.R1699W = c.5095C>T = NC_000017.11:g.43057051G>A.
   Solution: normalize to canonical protein HGVS and store the mapping.

3. **Coordinate system mismatches**: ClinVar uses genomic coordinates, UniProt
   uses canonical isoform numbering, AlphaFold structures use their own residue
   numbering (usually but not always matching UniProt). Solution: map through
   SIFTS (PDBe structure integration) or pairwise BLAST alignment, and validate
   by checking the reference amino acid at the mapped position.

4. **Isoform confusion**: TP53 has multiple isoforms (p53alpha, p53beta) with
   different domain compositions. The wrong isoform = wrong domain annotation.
   Solution: always use UniProt canonical isoform and flag any discrepancy.

**Validation rule**: After mapping, ALWAYS check that the reference amino acid at
the structure position matches the expected residue from the HGVS notation. If it
doesn't match, the mapping is wrong — set coordinate_mapping_confidence to "failed"
and record the warning. Do NOT proceed with structural analysis on a mismatched residue.

---

## 6. Tool-Specific Implementation Notes

### FoldX (ΔΔG calculation — the most important structural feature)

FoldX requires a license (free for academic use, request at foldxsuite.crg.eu).
Binary is not pip-installable — it's a standalone executable.

**Critical implementation steps:**
1. **RepairPDB first**: AlphaFold structures lack hydrogens and have suboptimal
   rotamers. ALWAYS run `FoldX --command=RepairPdb` before `BuildModel`.
2. **BuildModel command**: `FoldX --command=BuildModel --pdb=repaired.pdb --mutant-file=individual_list.txt`
3. **Mutant file format**: `RA1699W;` (one-letter ref AA, chain A, position, one-letter alt AA, semicolon)
4. **Output**: Look for `Dif_` file — the ΔΔG value is in the `total energy` column.
5. **Interpretation**: ΔΔG > 0 = destabilizing, > 2 = significantly destabilizing.
6. **Runtime**: ~2–5 minutes per mutation on CPU.
7. **If FoldX fails**: Fall back to PyRosetta, then DDGun web server, then skip ΔΔG.

### PyRosetta (Independent second ΔΔG)

Requires a license (free for academic, request at pyrosetta.org).
Install via: `pip install pyrosetta` (after license approval).

**Key difference from FoldX**: PyRosetta performs backbone relaxation after
introducing the mutation, allowing the protein to adjust. This gives a more
realistic ΔΔG estimate for some mutations. When FoldX and PyRosetta agree,
confidence in the ΔΔG value is high. When they disagree, flag for review.

### FreeSASA (Solvent accessibility — easiest structural tool)

`pip install freesasa` — always works, no license needed.

```python
import freesasa
structure = freesasa.Structure(pdb_path)
result = freesasa.calc(structure)
# Get per-residue SASA
residue_sasa = result.residueAreas()
```

Classify: < 10% relative SASA → "core", 10–40% → "interface", > 40% → "surface"

### DSSP (Secondary structure assignment)

Install `mkdssp` via system package or `pip install dssp`. BioPython has a DSSP wrapper.

```python
from Bio.PDB import PDBParser, DSSP
parser = PDBParser()
structure = parser.get_structure("protein", pdb_path)
dssp = DSSP(structure[0], pdb_path)
# dssp[(chain, (' ', position, ' '))] → (index, aa, ss, ...)
```

Secondary structure codes: H=helix, E=sheet, C=coil, T=turn, S=bend, G=3-10 helix.

### HMMER (Domain identification)

`pip install pyhmmer` for the Python bindings, or use the HMMER web server.
Requires Pfam database download (~2GB). Identifies which functional domain
(BRCT, kinase, SH3, etc.) the mutation position falls in.

Fallback: Use InterProScan web API, or simply look up UniProt domain annotations.

### CatBoost / XGBoost / LightGBM (The ML ensemble)

**Why three models**: Ensembles almost always outperform single models. Each model
has different strengths:
- **CatBoost**: Best on mixed feature types (numerical + categorical like domain_name).
  Handles missing values natively. Least hyperparameter tuning needed. **Use as primary.**
- **XGBoost**: Industry standard. Best SHAP integration. Most widely recognized.
- **LightGBM**: Fastest training. Often highest raw accuracy. Good for iteration.

**Why NOT deep learning**: Varis's core argument against AlphaMissense is that
black-box predictions aren't clinically trustworthy. Using a deep learning model
would undermine that narrative. Tree-based models with SHAP are interpretable —
a genetic counselor can see exactly which features drove the prediction. Research
(Grinsztajn et al., 2022) confirms gradient-boosted trees consistently outperform
deep learning on structured tabular data.

**Handling missing features**: All three models handle None/NaN natively. If FoldX
failed for a variant, ddg_foldx is None, and the models learn to use remaining
features. This is critical for graceful degradation — never impute missing values
with zeros or means for the final model.

**Feature availability indicators**: Always include feature_available_* flags
(ddg_available, sasa_available, etc.) as additional ML features alongside the
structural features. These let the model distinguish "tool was unavailable" from
"scientifically difficult variant." Use `record.get_ml_features()` which returns
both structural features and availability flags.

**Simulated missingness during training**: To prevent the model from learning
shortcuts from systematic missingness (e.g., FoldX failing on certain protein
families that happen to correlate with pathogenicity), simulate missingness during
training:
  1. For each training fold, randomly drop entire feature blocks (all structural,
     or all conservation, or ΔΔG only) for 10–20% of samples
  2. Set the corresponding feature_available_* flags to False for dropped features
  3. This forces the model to make good predictions even when features are absent
  4. After training, audit: does the model's performance on "feature-absent" subsets
     match its performance on "feature-present" subsets? Large gaps indicate the
     model over-relies on missingness patterns.

**Training data**: Expert-reviewed ClinVar variants labeled as Pathogenic (1) or
Benign (0). Start with 1,000–2,000, scale to 10K–200K. Use stratified k-fold
cross-validation and Optuna for hyperparameter tuning.

### Evaluation Protocol (Critical for Credibility)

The ROC-AUC numbers in the build plan (0.80–0.93+) are TARGETS, not expectations.
Making them credible requires a rigorous evaluation protocol:

**Split strategies** (implement all three, report the hardest):

1. **Gene-stratified split**: Hold out entire genes from training. If BRCA1 is in
   the test set, NO BRCA1 variants appear in training. This tests generalization
   across genes. The hardest and most honest evaluation — expect 5–10 points lower
   than random split.

2. **Time-split validation**: Train on variants classified before a cutoff date,
   test on variants classified after. This simulates real-world use (predicting
   new variants). Use ClinVar submission dates for the split.

3. **Random stratified split**: Standard k-fold with class balance. The easiest
   evaluation — include for comparison but do NOT lead with these numbers.

**Metrics to report** (minimum set):
- ROC-AUC (discrimination)
- PR-AUC (performance on imbalanced classes — pathogenic variants are rarer)
- Precision and recall at the classification threshold
- Calibration plot (does score=0.8 really mean 80% probability?)
- Performance by feature availability (full features vs. lite features vs. minimal)

**Failure modes table**: For every evaluation, also report:
- Performance on "hard genes" (genes with few training examples)
- Performance on low-confidence structures (pLDDT < 70)
- Performance when ΔΔG is unavailable vs. available
- Variants where the ensemble disagrees (model_agreement = "low")

Never present random-split ROC-AUC as the headline number. Lead with
gene-stratified results — they're lower but honest, and reviewers will respect it.

**Model disagreement**: When the three models' predictions differ by more than 0.2,
flag the variant for closer review. This is a novel contribution — disagreement
detection signals genuinely uncertain cases.

---

## 7. The Self-Evolution Architecture

### Loop 1: Auto-Retraining (Fully Autonomous)

Monthly job:
1. Scan ClinVar for newly expert-reviewed variants (VUS → Pathogenic or Benign)
2. Run each new labeled variant through M1–M4 to extract structural features
3. Add to training dataset
4. Retrain all three models with Optuna tuning
5. Benchmark new model vs. current model on held-out test set
6. Run automated regression tests on fixed benchmark variant set (see below)
7. Deploy ONLY if ALL key metrics improve (ROC-AUC, PR-AUC, precision, recall)
   AND regression tests pass
8. If new model is worse → keep current model and log why
9. Tag and archive the model version (deployed or rejected)
10. Record everything in the Evolution Log

**Key implementation detail**: The decision to deploy must be conservative. A model
that improves ROC-AUC but drops precision should NOT be deployed — it would increase
false positives, which in clinical genetics means telling a patient their variant is
pathogenic when it isn't.

### Release Governance (Critical for Clinical Trust)

Varis follows a strict release discipline to prevent regressions from reaching users:

**Model versioning**: Every trained model is tagged with a date stamp (e.g.
v2026.03, v2026.04). Previous versions are preserved in a version archive. The
current production model, the candidate model, and at least 6 months of history
must always be available.

**Regression test suite**: A fixed benchmark set of well-characterised variants
(minimum 50, spanning pathogenic, benign, and known edge cases across diverse
genes) is re-evaluated with every candidate model. Any unexpected reclassification
of a benchmark variant triggers manual review before deployment. The benchmark set
is stored in `tests/benchmark_variants.json` and must never be included in
training data.

**Rollback protocol**: If a deployed model produces unexpected results in
production (flagged by user reports or automated monitoring), the previous version
can be restored immediately from the version archive. The rollback is recorded in
the Evolution Log with the reason.

**Changelog with metric deltas**: Every deployment entry in the Evolution Log must
include: old model version, new model version, metric comparisons (with exact
numbers), number of new training variants added, regression test results, and
deploy/reject decision with rationale.

Implementation note: Build the version archive and regression test infrastructure
BEFORE enabling auto-retrain. The governance framework is what makes monthly
retraining safe enough for clinical-adjacent use.

### Loop 2: Tool Discovery (Semi-Autonomous)

A scout subsystem monitors arXiv, GitHub, HuggingFace, PyPI, and BioRxiv for new
tools. When it finds a candidate, it sends the paper abstract + Varis's pipeline
description to an LLM (Claude API) and asks: "Would this tool improve any stage of
Varis's pipeline? If so, which stage and why?"

The LLM generates a structured proposal. This is a stretch goal (Priority 5).

### Loop 3: Auto-Integration (Partially Autonomous)

When a proposed tool passes review: install → generate wrapper → benchmark → integrate.
Works autonomously for ~60–70% of standard Python tools. The rest need human help.
This is also a stretch goal (Priority 6).

### The Evolution Log

Every action Varis takes on itself is recorded in a structured, timestamped, public
log. Event types: AUTO_RETRAIN, TOOL_DISCOVERY, LLM_ASSESSMENT, TOOL_INTEGRATION,
MODEL_DEPLOY. Published on VarisDB for full transparency.

**This is unprecedented for a bioinformatics tool** — no other tool publicly logs
every decision it makes about its own improvement. Implement the Evolution Log
BEFORE the auto-retrain loop, because the log records everything the loop does.

---

## 8. The Validation Variant: BRCA1 p.Arg1699Trp

This is the canonical test case used throughout development. Every module should
be tested against this variant because the expected results are well-documented:

| Property | Expected Value | Source |
|----------|---------------|--------|
| Gene | BRCA1 | — |
| Position | 1699 | — |
| Ref AA | Arg (R) — positively charged | — |
| Alt AA | Trp (W) — neutral, bulky | — |
| ClinVar ID | VCV000055361 | ClinVar |
| ClinVar status | Uncertain significance (historically) | ClinVar |
| AlphaMissense | 0.934 (likely_pathogenic) | Cheng et al., 2023 |
| Domain | BRCT (phosphopeptide binding) | Pfam PF00533 |
| ΔΔG | ~3.5–4.5 kcal/mol (destabilizing) | Williams et al., 2004 |
| Burial | Core (~3% exposed) | FreeSASA |
| Secondary structure | Helix | DSSP |
| Conservation | 100% across mammals | BLAST + Clustal Omega |
| Charge change | Positive → neutral | Sequence |
| Expected classification | Likely Pathogenic | All evidence combined |
| ACMG codes expected | PM1 (BRCT domain), PP3 (computational), PS3-proxy (ΔΔG > 2) | — |

If your implementation produces results significantly different from these values,
something is wrong. Debug before moving on.

---

## 9. Open Science Philosophy

Varis follows a radical openness strategy. Everything is published:

- **VarisDB**: Free public database at varisdb.russellgenetics.org
- **ML ensemble**: Weights on HuggingFace under Russell Genetics
- **Structural feature dataset**: 200K+ variants × 15 features on Zenodo with DOI
- **Pipeline code**: Complete codebase on GitHub under Russell Genetics
- **Educational notebooks**: 5 Google Colab tutorials
- **Evolution Log**: Public, auditable record of every change
- **ClinVar submissions**: Structural evidence contributed to the global variant record

**The competitive advantages that can't be copied** (even when everything is public):
1. First-mover data lead (VarisDB grows daily)
2. Community trust and reputation (takes years to build)
3. Continuous retraining (static copies degrade over time)
4. ClinVar submitter status (requires track record and NIH review)
5. Russell's domain expertise (in his head, not in the code)
6. Russell Genetics brand (the field remembers who built the tool)

Never implement anything that contradicts this openness. No paywalls, no proprietary
features, no secret sauce. The science is free. The infrastructure is sustainable
through grants and institutional partnerships.

### Licensing Clarification (Important for Documentation and README)

Varis's own code is released under the MIT license — fully open, no restrictions.
However, some structural analysis dependencies have their own licensing terms:

- **FoldX**: Free academic license required from CRG Barcelona. Non-redistributable.
  Users must obtain their own license at https://foldxsuite.crg.eu/
- **PyRosetta**: Free academic license required from RosettaCommons. Non-redistributable.
  Users must apply at https://www.rosettacommons.org/software/license-and-download

These tools CANNOT be bundled or redistributed with Varis. The README, installation
docs, and any public-facing materials must clearly state this. When writing
documentation, always say "Varis's code is open source under MIT license. Some
structural analysis tools require separate free academic licenses."

Varis's fallback architecture (M3 has multiple tools per analysis step) means the
pipeline still functions if a particular tool is unavailable, though with reduced
feature depth. This is by design — graceful degradation applies to licensing
constraints as well as runtime failures.

---

## 10. What "Done" Looks Like at Each Phase

### Phase 1 (Weeks 1–3): Data Foundation
**Done when:** `python -m varis BRCA1 p.Arg1699Trp` produces a JSON file with gene
info, protein sequence, ClinVar status, gnomAD frequency, AlphaMissense score, and
a downloaded PDB file. All fields populated from real API calls.

### Phase 2 (Weeks 3–6): Structural Analysis
**Done when:** The JSON now includes ΔΔG (from FoldX or PyRosetta or both), solvent
accessibility (FreeSASA), secondary structure (DSSP), contacts (BioPython), and
domain identification (HMMER). For BRCA1 p.Arg1699Trp, ΔΔG should be ~3.5–4.5.

### Phase 3 (Weeks 5–7): Conservation
**Done when:** The JSON includes conservation score (close to 1.0 for BRCA1 position
1699) from BLAST + Clustal Omega alignment. Independently implementable alongside Phase 2.

### Phase 4 (Weeks 6–9): ML Ensemble
**Done when:** Three models trained on 1,000+ ClinVar variants, SHAP working, ensemble
predicts BRCA1 p.Arg1699Trp as likely_pathogenic with high confidence. Benchmark
comparison against AlphaMissense on held-out test set.

### Phase 5 (Weeks 8–11): VarisDB
**Done when:** A web interface at localhost where Russell can search for a variant and
see the full investigation report with charts, SHAP plot, and (ideally) 3D structure.

#### VarisDB Investigation UI Specification

The VarisDB frontend is the public face of every investigation. It must be **honest**
— every visual element must reflect real computed data, not decorative placeholders.

**Architecture (Backend to Frontend contract)**

The backend (FastAPI) serves a single `/api/v1/investigate/{variant_id}` endpoint
that returns a JSON payload with four sections:

- `structure`: source (alphafold/esmfold/experimental), url_or_data, chain, residue_index,
  ref_aa, alt_aa, plddt_at_residue, coordinate_mapping_confidence (exact/high/low/failed),
  normalization_warnings
- `features`: list of { name, value, units, evidence_tag, available } — from
  get_structural_features() plus feature_available flags
- `prediction`: score, classification, confidence_lower, confidence_upper, model_agreement,
  individual_scores (catboost/xgboost/lightgbm), uncertainty_flags
- `explanation`: list of { feature, value, shap } sorted by |shap| descending — computed
  server-side by the ACTUAL trained CatBoost/XGBoost/LightGBM models

CRITICAL: The `explanation` array MUST come from actual SHAP values computed server-side.
The React frontend is a pure visualization layer — no ML computation in the browser. A
pretty waterfall disconnected from the computed features will backfire in a demo.

**Mapping from VariantRecord to API response:**

The `/investigate` endpoint reads directly from the VariantRecord and maps fields:
- `structure` section ← `normalization.*` fields (canonical_transcript,
  structure_residue_position, coordinate_mapping_confidence, normalization_warnings)
  plus `structural_features.plddt_score`
- `features` section ← `get_structural_features()` output + `feature_availability.*` flags
- `prediction` section ← `scoring.*` fields (ensemble_score, classification,
  confidence_lower, confidence_upper, model_agreement, individual model scores)
- `explanation` section ← `scoring.shap_top_features` (list of dicts with feature,
  value, shap — sorted by |shap| descending). This is computed by the SHAP explainer
  in M5 using the actual trained CatBoost/XGBoost/LightGBM models.

**Frontend Layout (React + Tailwind)**

The investigation page has three zones:

**Left panel: Mol* 3D Protein Viewer**
- Embed pdbe-molstar (the viewer used by PDB and AlphaFold DB)
- On load: color entire structure by pLDDT confidence
- Highlight the mutation site residue in red/orange
- Hover tooltip: "p.Arg1699Trp (Chain A: 1699) — pLDDT: 92.4"
- Mapping robustness (critical):
  - If coordinate_mapping_confidence is "exact" or "high": highlight residue normally
  - If "low": show yellow banner "Residue mapping uncertain (isoform mismatch).
    Showing closest match via sequence alignment." Highlight with dashed outline.
  - If "failed": show red banner "Coordinate mapping failed. 3D context unavailable."
    Do NOT highlight any residue. Never silently highlight the wrong residue.

**Right panel: Evidence and Prediction**
- **Reliability Strip** (compact bar, always visible at top):
  - Structure confidence badge: green (pLDDT > 90), yellow (70-90), red (< 70)
  - Feature availability indicators: icons per block (Structure, Conservation, DDG)
    with null_reason on hover (e.g., "tool_missing", "license_unavailable")
  - Ensemble agreement: green check or orange warning
  - This strip prevents overclaiming and is the UI expression of null reason codes
    and feature_available flags from the Variant Record schema

- **SHAP Waterfall Chart** (animated):
  - Start from baseline (population average), show each feature pushing the score
    toward pathogenic (right/red) or benign (left/blue)
  - Animate features appearing one by one (top SHAP magnitude first)
  - Each bar labeled: feature name + actual value + SHAP contribution
  - Final score shown at bottom with confidence interval band
  - Implementation: Recharts first (stacked bar with invisible baselines). Upgrade
    to D3 only if Recharts rendering is insufficient.

- **ACMG Evidence Panel**:
  - Show suggested evidence codes (PM1, PP3, PS3-proxy, etc.)
  - Header says "Suggested Evidence Codes" — never "Assigned" or "Determined"
  - Each code shows the supporting measurement (e.g., PM1: "BRCT domain, Pfam PF00533")

### Phase 6 (Weeks 10–14): Self-Evolution
**Done when:** Auto-retrain loop runs monthly, Evolution Log records events, ClinVar
submission formatter produces valid submissions, codebase is packaged on GitHub.

---

## 11. Common Mistakes to Avoid

1. **Using ΔΔG = 0 for missing values**: Use None. The tree models handle it natively.
   Imputing 0 tells the model "this mutation has no effect on stability," which is wrong.

2. **Treating AlphaMissense as ground truth**: It's one feature among 15. The ensemble
   should be able to disagree with AlphaMissense when structural evidence says otherwise.

3. **Ignoring pLDDT**: If pLDDT < 70 at the mutation site, all structural features
   from that region are unreliable. Flag this in the output.

4. **Running FoldX without RepairPDB first**: AlphaFold structures need hydrogen atoms
   and rotamer optimization before FoldX can calculate accurate ΔΔG.

5. **Hardcoding API responses for tests**: Use real API calls for integration tests,
   mock responses for unit tests. The conftest.py fixtures provide pre-populated
   VariantRecords for unit testing without API calls.

6. **Making M3 depend on M4 or vice versa**: They are completely independent. M3 uses
   the 3D structure (from M2). M4 uses the protein sequence (from M1). Never import
   between them.

7. **Using deep learning anywhere**: The entire narrative is built on interpretability.
   A single neural network layer anywhere in the pipeline undermines the argument
   against AlphaMissense's black-box approach.

8. **Forgetting to log failures**: Every tool failure must be recorded in
   modules_failed on the VariantRecord. This is how graceful degradation is auditable.

9. **Setting fields to None without a reason code**: Every null field MUST have a
   corresponding entry in `null_reasons` (e.g., `"ddg_foldx": "tool_missing"`).
   Use `record.set_with_reason(field, value, reason)` instead of direct assignment.
   Use `record.set_feature_status(group, available, reason)` for feature availability
   flags. Without reason codes, you cannot distinguish a pipeline bug from a
   scientifically difficult variant — and the ML model may learn artifacts.

10. **Mixing records from different schema versions**: Every VariantRecord is stamped
    with `record_schema_version`. When adding fields to the schema, increment the
    version. Before training, verify all records in the training set have the same
    schema version. Records from older schemas must be migrated or excluded.

11. **Skipping coordinate validation after mapping**: After mapping residue positions
    between databases, ALWAYS verify the reference amino acid matches. If UniProt says
    position 1699 is Arg but the AlphaFold structure has Gly at that position, the
    mapping is wrong. Set confidence to "failed" and do NOT run structural analysis.

12. **Spending more than 3 days on one tool**: Switch to the fallback. Document what
    failed and why. Move on. A working system with 8 features beats a broken system
    that was supposed to have 15.

13. **Skipping the Evolution Log**: Build it before auto-retrain. The log is what
    makes self-evolution auditable and transparent. Without it, auto-retrain is just
    a cron job — with it, it's unprecedented scientific infrastructure.

---

## 12. The Narrative in One Paragraph

Varis is an open-source structural investigator for genetic variants, built by
Russell Young as a Russell Genetics product. While Google DeepMind's AlphaMissense
gives clinicians a single frozen score from a black-box model, Varis investigates
why a mutation is damaging — orchestrating a 41-tool pipeline to calculate protein
destabilization, map functional domains, measure evolutionary conservation, and
generate multiple ACMG evidence lines. The ML ensemble (CatBoost + XGBoost +
LightGBM) is fully interpretable via SHAP. The platform auto-retrains on new ClinVar
data, discovers new tools by scanning arXiv, and logs every decision in a public
Evolution Log. Everything is open-source — the model, the data, the pipeline, the
database — because the fastest way to help rare disease patients is to put the tools
in every lab's hands.
