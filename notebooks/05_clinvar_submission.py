# %% [markdown]
# # Notebook 5: Submitting Structural Evidence to ClinVar
#
# Format Varis investigation results into ClinVar submission format
# and contribute to the global variant record.
#
# ## Overview
#
# ClinVar is NIH's central database for genotype-phenotype relationships.
# Varis can format its structural investigation results into ClinVar
# submissions to contribute computational evidence.
#
# **Important:** All submissions are labelled as computational suggestions.
# Varis does NOT make clinical adjudications.

# %%
from varis.models.variant_record import create_variant_record
from varis.m6_platform.api.clinvar_submitter import (
    format_clinvar_submission,
    submit_to_clinvar,
)

# %% [markdown]
# ## 1. Load a completed investigation

# %%
# Example: load a completed VariantRecord
# record = VariantRecord.load("data/training/variants/BRCA1_p.Arg1699Trp.json")

# For demo, create a mock record with typical pathogenic results
record = create_variant_record("BRCA1", "p.Arg1699Trp")
record.score_ensemble = 0.92
record.classification = "likely_pathogenic"
record.confidence_lower = 0.85
record.confidence_upper = 0.96
record.model_agreement = "high"
record.evidence_tags = ["computational_support", "energetics_support", "domain_context"]
record.evidence_computational_support = True
record.evidence_energetics = True
record.evidence_domain_context = True
record.ddg_mean = 4.2
record.burial_category = "core"
record.domain_name = "BRCT"
record.domain_id = "PF00533"
record.conservation_score = 0.98
record.modules_completed = ["M1", "M2", "M3", "M4", "M5"]
record.ensemble_version = "v2026.03"
record.canonical_transcript = "NM_007294.4"
record.hgvs_coding = "c.5095C>T"

# %% [markdown]
# ## 2. Format the submission

# %%
submission = format_clinvar_submission(record)
if submission:
    print("Eligible for ClinVar submission!")
    print(f"Classification: {submission['clinical_significance']}")
    print(f"Gene: {submission['gene_symbol']}")
    print(f"HGVS: {submission['hgvs_coding']}")
    print(f"\nEvidence comment:\n{submission['comment']}")
else:
    print("Not eligible for ClinVar submission")
    print(f"Classification: {record.classification}")
    print(f"Model agreement: {record.model_agreement}")

# %% [markdown]
# ## 3. Dry-run submission (does NOT call ClinVar API)

# %%
if submission:
    result = submit_to_clinvar(submission, dry_run=True)
    print(f"Status: {result['status']}")
    print(f"Variant: {result['variant_id']}")

# %% [markdown]
# ## 4. Real submission (requires CLINVAR_API_KEY)
#
# To submit for real:
# 1. Set CLINVAR_API_KEY environment variable
# 2. Register Russell Genetics as a ClinVar submitter
# 3. Use use_test_endpoint=True first to validate
# 4. Then use_test_endpoint=False for production
#
# ```python
# result = submit_to_clinvar(submission, dry_run=False, use_test_endpoint=True)
# ```
