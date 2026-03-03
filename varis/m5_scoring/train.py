"""M5 Training Pipeline — ClinVar download, feature computation, and model training.

Three-phase pipeline:
  Phase A (select):  Download ClinVar, select balanced missense variants per gene.
  Phase B (compute): Run M1-M4 on each variant, cache VariantRecords. Resumable.
  Phase C (train):   Build feature matrix, gene-stratified CV, train final ensemble.

Usage:
    python -m varis.m5_scoring select          # Phase A only
    python -m varis.m5_scoring compute         # Phase B only (resumable)
    python -m varis.m5_scoring train-only      # Phase C only (from cached features)
    python -m varis.m5_scoring train           # All phases
    python -m varis.m5_scoring status          # Show progress
"""

import argparse
import copy
import gzip
import json
import logging
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import httpx
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

from varis.config import (
    AA_ONE_TO_THREE,
    DATA_DIR,
    MODELS_DIR,
    PIPELINE_VERSION,
)
from varis.models.variant_record import NullReason, VariantRecord, create_variant_record

logger = logging.getLogger(__name__)

# =============================================================================
# CONSTANTS
# =============================================================================

TRAINING_CACHE_DIR = DATA_DIR / "training"
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_LOCAL_PATH = DATA_DIR / "variant_summary.txt.gz"

BENCHMARK_VARIANTS = {
    ("BRCA1", "p.Arg1699Trp"),
    ("BRCA1", "p.Lys1183Arg"),
    ("TP53", "p.Arg175His"),
    ("BRCA2", "p.Asp2723His"),
    ("CFTR", "p.Gly551Asp"),
}

TARGET_GENES = [
    "BRCA1", "BRCA2", "TP53", "CFTR", "MSH2", "MLH1", "MSH6", "PMS2",
    "PTEN", "RB1", "APC", "VHL", "MEN1", "RET", "CDH1", "PALB2",
    "ATM", "CHEK2", "RAD51C", "RAD51D",
]

MAX_PER_CLASS_PER_GENE = 15
MIN_PER_CLASS_PER_GENE = 5
API_SLEEP_SECONDS = 2.0
MAX_CONSECUTIVE_FAILURES = 5

HGVS_PROTEIN_RE = re.compile(r'\(p\.([A-Z][a-z]{2}\d+[A-Z][a-z]{2})\)')
HGVS_SINGLE_RE = re.compile(r'p\.([A-Z]\d+[A-Z])')

VALID_AA_THREE = {
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
}

# Map ClinVar ReviewStatus strings to star counts
REVIEW_STATUS_STARS = {
    "no assertion criteria provided": 0,
    "no assertion for the individual variant": 0,
    "no assertion provided": 0,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}


# =============================================================================
# PHASE A: VARIANT SELECTION
# =============================================================================

def download_clinvar(force: bool = False) -> Path:
    """Download ClinVar variant_summary.txt.gz if not cached locally.

    Checks whether the local file exists and is less than 7 days old.
    If not, downloads a fresh copy from NCBI FTP via httpx.

    Args:
        force: If True, re-download even if a recent cached copy exists.

    Returns:
        Path to the local variant_summary.txt.gz file.
    """
    CLINVAR_LOCAL_PATH.parent.mkdir(parents=True, exist_ok=True)

    if not force and CLINVAR_LOCAL_PATH.exists():
        age_seconds = time.time() - CLINVAR_LOCAL_PATH.stat().st_mtime
        if age_seconds < 7 * 24 * 3600:
            logger.info(
                "ClinVar file is %.1f days old, using cached copy: %s",
                age_seconds / 86400,
                CLINVAR_LOCAL_PATH,
            )
            return CLINVAR_LOCAL_PATH

    logger.info("Downloading ClinVar variant_summary.txt.gz ...")
    try:
        with httpx.Client(timeout=300, follow_redirects=True) as client:
            with client.stream("GET", CLINVAR_URL) as response:
                response.raise_for_status()
                with open(CLINVAR_LOCAL_PATH, "wb") as f:
                    for chunk in response.iter_bytes(chunk_size=65536):
                        f.write(chunk)
        logger.info("Downloaded ClinVar to %s", CLINVAR_LOCAL_PATH)
    except Exception as e:
        logger.error("Failed to download ClinVar: %s", e)
        if CLINVAR_LOCAL_PATH.exists():
            logger.warning("Using stale cached copy")
        else:
            raise
    return CLINVAR_LOCAL_PATH


def _extract_hgvs_protein(name: str) -> Optional[str]:
    """Extract p.XxxNNNYyy from ClinVar Name column.

    Tries three-letter format first: (p.Arg1699Trp).
    Falls back to single-letter: p.R1699W and converts to three-letter.
    Validates that ref and alt are valid amino acid codes.

    Args:
        name: The ClinVar Name field value (e.g.,
              "NM_007294.4(BRCA1):c.5095C>T (p.Arg1699Trp)").

    Returns:
        Three-letter HGVS string like "p.Arg1699Trp", or None if unparseable.
    """
    if not isinstance(name, str):
        return None

    # Try three-letter: (p.Arg1699Trp)
    match = HGVS_PROTEIN_RE.search(name)
    if match:
        hgvs = match.group(1)
        # Parse ref and alt AA codes
        aa_match = re.match(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs)
        if aa_match:
            ref_aa, _, alt_aa = aa_match.group(1), aa_match.group(2), aa_match.group(3)
            if ref_aa in VALID_AA_THREE and alt_aa in VALID_AA_THREE:
                return f"p.{hgvs}"
        return None

    # Try single-letter: p.R1699W
    match = HGVS_SINGLE_RE.search(name)
    if match:
        hgvs_single = match.group(1)
        single_match = re.match(r'([A-Z])(\d+)([A-Z])', hgvs_single)
        if single_match:
            ref_single = single_match.group(1)
            position = single_match.group(2)
            alt_single = single_match.group(3)
            ref_three = AA_ONE_TO_THREE.get(ref_single)
            alt_three = AA_ONE_TO_THREE.get(alt_single)
            if ref_three and alt_three:
                return f"p.{ref_three}{position}{alt_three}"
        return None

    return None


def _slugify_variant(gene: str, hgvs: str) -> str:
    """Create filesystem-safe cache filename.

    Args:
        gene: Gene symbol (e.g., "BRCA1").
        hgvs: HGVS protein notation (e.g., "p.Arg1699Trp").

    Returns:
        Sanitized string safe for use as a filename.
    """
    raw = f"{gene}_{hgvs}"
    return re.sub(r'[^A-Za-z0-9._-]', '_', raw)


def _review_status_to_stars(review_status: str) -> int:
    """Convert ClinVar ReviewStatus string to star count.

    Args:
        review_status: ClinVar ReviewStatus field value.

    Returns:
        Integer star count (0-4). Defaults to 0 for unrecognised strings.
    """
    if not isinstance(review_status, str):
        return 0
    status_lower = review_status.strip().lower()
    for pattern, stars in REVIEW_STATUS_STARS.items():
        if pattern in status_lower:
            return stars
    return 0


def select_training_variants(
    clinvar_path: Optional[Path] = None,
    target_genes: Optional[list[str]] = None,
    max_per_class: int = MAX_PER_CLASS_PER_GENE,
    min_per_class: int = MIN_PER_CLASS_PER_GENE,
    min_stars: int = 1,
    seed: int = 42,
) -> pd.DataFrame:
    """Select balanced missense variants for training from ClinVar.

    Label hygiene rules:
      - Keep only: Pathogenic, Likely pathogenic, Benign, Likely benign
      - Exclude: VUS, conflicting interpretations
      - Prefer review status >= min_stars

    Args:
        clinvar_path: Path to variant_summary.txt.gz. Downloaded if None.
        target_genes: Gene list to select from. Defaults to TARGET_GENES.
        max_per_class: Max variants per class (pathogenic/benign) per gene.
        min_per_class: Min variants per class per gene. Gene skipped if unmet.
        min_stars: Minimum ClinVar review stars required.
        seed: Random seed for reproducible sampling.

    Returns:
        DataFrame with columns: gene, hgvs, label, clinvar_classification,
        clinvar_review_status, cache_path, computed, computed_at, failed_reason.
    """
    if target_genes is None:
        target_genes = TARGET_GENES

    if clinvar_path is None:
        clinvar_path = download_clinvar()

    TRAINING_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("Reading ClinVar file: %s", clinvar_path)
    try:
        df = pd.read_csv(clinvar_path, sep='\t', low_memory=False)
    except Exception as e:
        logger.error("Failed to read ClinVar file: %s", e)
        return pd.DataFrame()

    logger.info("ClinVar loaded: %d rows", len(df))

    # Filter to GRCh38, single nucleotide variant
    if "Assembly" in df.columns:
        df = df[df["Assembly"] == "GRCh38"]
    if "Type" in df.columns:
        df = df[df["Type"] == "single nucleotide variant"]

    # Filter clinical significance
    if "ClinicalSignificance" not in df.columns:
        logger.error("ClinicalSignificance column not found in ClinVar file")
        return pd.DataFrame()

    sig_col = df["ClinicalSignificance"].fillna("")

    # Pathogenic / Likely pathogenic (exclude conflicting)
    path_mask = (
        sig_col.str.contains("Pathogenic", case=False, na=False)
        & ~sig_col.str.contains("Conflicting", case=False, na=False)
        & ~sig_col.str.contains("Benign", case=False, na=False)
    )
    # Benign / Likely benign (exclude conflicting)
    benign_mask = (
        sig_col.str.contains("Benign", case=False, na=False)
        & ~sig_col.str.contains("Conflicting", case=False, na=False)
        & ~sig_col.str.contains("Pathogenic", case=False, na=False)
    )
    df = df[path_mask | benign_mask].copy()
    df["label"] = path_mask[df.index].astype(int)

    logger.info("After significance filter: %d variants", len(df))

    # Filter to target genes
    if "GeneSymbol" in df.columns:
        df = df[df["GeneSymbol"].isin(target_genes)]
    else:
        logger.error("GeneSymbol column not found")
        return pd.DataFrame()

    logger.info("After gene filter (%d genes): %d variants", len(target_genes), len(df))

    # Filter by review status stars
    if "ReviewStatus" in df.columns:
        df["_stars"] = df["ReviewStatus"].apply(_review_status_to_stars)
        df = df[df["_stars"] >= min_stars]
        logger.info("After review status filter (>= %d stars): %d variants", min_stars, len(df))

    # Extract HGVS protein from Name column
    if "Name" not in df.columns:
        logger.error("Name column not found in ClinVar file")
        return pd.DataFrame()

    df["hgvs"] = df["Name"].apply(_extract_hgvs_protein)
    n_before = len(df)
    df = df.dropna(subset=["hgvs"])
    logger.info(
        "Extracted HGVS protein: %d/%d variants had parseable protein notation",
        len(df), n_before,
    )

    # Deduplicate by (gene, hgvs) — keep highest star rating
    df = df.sort_values("_stars", ascending=False).drop_duplicates(
        subset=["GeneSymbol", "hgvs"], keep="first",
    )

    # Exclude benchmark variants
    df["_is_benchmark"] = df.apply(
        lambda row: (row["GeneSymbol"], row["hgvs"]) in BENCHMARK_VARIANTS,
        axis=1,
    )
    n_benchmarks = df["_is_benchmark"].sum()
    if n_benchmarks > 0:
        logger.info("Excluded %d benchmark variants", n_benchmarks)
    df = df[~df["_is_benchmark"]]

    # Select up to max_per_class per class per gene
    rng = np.random.RandomState(seed)
    selected_rows = []
    skipped_genes = []

    for gene in target_genes:
        gene_df = df[df["GeneSymbol"] == gene]
        pathogenic = gene_df[gene_df["label"] == 1]
        benign = gene_df[gene_df["label"] == 0]

        if len(pathogenic) < min_per_class or len(benign) < min_per_class:
            skipped_genes.append(gene)
            logger.info(
                "Skipping %s: pathogenic=%d, benign=%d (need >= %d each)",
                gene, len(pathogenic), len(benign), min_per_class,
            )
            continue

        # Sample up to max_per_class from each class
        n_path = min(max_per_class, len(pathogenic))
        n_benign = min(max_per_class, len(benign))

        path_sample = pathogenic.sample(n=n_path, random_state=rng)
        benign_sample = benign.sample(n=n_benign, random_state=rng)

        selected_rows.append(path_sample)
        selected_rows.append(benign_sample)
        logger.info(
            "Selected %s: %d pathogenic, %d benign",
            gene, n_path, n_benign,
        )

    if skipped_genes:
        logger.warning("Skipped %d genes: %s", len(skipped_genes), ", ".join(skipped_genes))

    if not selected_rows:
        logger.error("No variants selected! Check ClinVar file and gene list.")
        return pd.DataFrame()

    result = pd.concat(selected_rows, ignore_index=True)

    # Build output DataFrame
    cache_dir = TRAINING_CACHE_DIR / "variants"
    cache_dir.mkdir(parents=True, exist_ok=True)

    output = pd.DataFrame({
        "gene": result["GeneSymbol"].values,
        "hgvs": result["hgvs"].values,
        "label": result["label"].values,
        "clinvar_classification": result["ClinicalSignificance"].values,
        "clinvar_review_status": result["ReviewStatus"].values,
        "cache_path": [
            str(cache_dir / f"{_slugify_variant(g, h)}.json")
            for g, h in zip(result["GeneSymbol"], result["hgvs"])
        ],
        "computed": False,
        "computed_at": None,
        "failed_reason": None,
    })

    # Save manifest
    manifest_path = TRAINING_CACHE_DIR / "variants_manifest.csv"
    output.to_csv(manifest_path, index=False)
    logger.info("Saved manifest with %d variants to %s", len(output), manifest_path)

    # Save benchmark variants separately for reference
    benchmark_df = pd.DataFrame(
        list(BENCHMARK_VARIANTS), columns=["gene", "hgvs"],
    )
    benchmarks_path = TRAINING_CACHE_DIR / "benchmarks.csv"
    benchmark_df.to_csv(benchmarks_path, index=False)
    logger.info("Saved benchmark variants to %s", benchmarks_path)

    # Summary
    logger.info(
        "Selection summary: %d variants from %d genes "
        "(pathogenic=%d, benign=%d)",
        len(output),
        output["gene"].nunique(),
        (output["label"] == 1).sum(),
        (output["label"] == 0).sum(),
    )

    return output


# =============================================================================
# PHASE B: FEATURE COMPUTATION
# =============================================================================

def _run_m1_through_m4(
    gene: str,
    hgvs: str,
    skip_m4: bool = False,
) -> Optional[VariantRecord]:
    """Run M1-M4 pipeline for a single variant. Does NOT run M5.

    Follows the pipeline.py pattern with try/except around each module.
    Modules are imported lazily to avoid circular imports.

    Args:
        gene: Gene symbol (e.g., "BRCA1").
        hgvs: HGVS protein notation (e.g., "p.Arg1699Trp").
        skip_m4: If True, skip conservation analysis (faster).

    Returns:
        VariantRecord with M1-M4 fields populated, or None on total failure.
    """
    import varis.m1_ingestion as m1
    import varis.m2_structure as m2
    import varis.m3_structural_analysis as m3
    import varis.m4_conservation as m4

    record = create_variant_record(gene, hgvs)
    record.pipeline_version = PIPELINE_VERSION

    try:
        record = m1.run(record)
    except Exception as e:
        logger.warning("M1 failed for %s %s: %s", gene, hgvs, e)
        record.mark_module_failed("M1")

    try:
        record = m2.run(record)
    except Exception as e:
        logger.warning("M2 failed for %s %s: %s", gene, hgvs, e)
        record.mark_module_failed("M2")

    try:
        record = m3.run(record)
    except Exception as e:
        logger.warning("M3 failed for %s %s: %s", gene, hgvs, e)
        record.mark_module_failed("M3")

    if not skip_m4:
        try:
            record = m4.run(record)
        except Exception as e:
            logger.warning("M4 failed for %s %s: %s", gene, hgvs, e)
            record.mark_module_failed("M4")

    return record


def compute_features(
    manifest_path: Optional[Path] = None,
    sleep_seconds: float = API_SLEEP_SECONDS,
    max_variants: Optional[int] = None,
    skip_m4: bool = False,
) -> tuple[int, int, int, int]:
    """Run M1-M4 on each variant in manifest. Resumable.

    Reads the manifest CSV, processes each uncomputed variant, saves the
    VariantRecord JSON to its cache_path, and updates the manifest.
    Tracks consecutive failures per gene and skips a gene after
    MAX_CONSECUTIVE_FAILURES.

    Args:
        manifest_path: Path to variants_manifest.csv. Uses default if None.
        sleep_seconds: Seconds to sleep between API calls.
        max_variants: Maximum number of variants to compute. None = all.
        skip_m4: If True, skip conservation (M4) for faster processing.

    Returns:
        Tuple of (total, cached, computed, failed) counts.
    """
    if manifest_path is None:
        manifest_path = TRAINING_CACHE_DIR / "variants_manifest.csv"

    if not manifest_path.exists():
        logger.error("Manifest not found: %s. Run 'select' first.", manifest_path)
        return (0, 0, 0, 0)

    manifest = pd.read_csv(manifest_path)
    total = len(manifest)
    cached = int(manifest["computed"].fillna(False).sum())
    computed = 0
    failed = 0
    variants_processed = 0

    # Track consecutive failures per gene
    gene_failures: dict[str, int] = {}
    skipped_genes: set[str] = set()

    logger.info(
        "Starting feature computation: %d total, %d already cached, %d pending",
        total, cached, total - cached,
    )

    for idx, row in manifest.iterrows():
        if max_variants is not None and variants_processed >= max_variants:
            logger.info("Reached max_variants limit (%d), stopping.", max_variants)
            break

        # Skip already-computed variants
        if row.get("computed") is True or row.get("computed") == "True":
            continue

        gene = row["gene"]
        hgvs = row["hgvs"]

        # Skip genes with too many consecutive failures
        if gene in skipped_genes:
            continue

        # Check if cache file already exists (in case manifest is stale)
        cache_path = Path(row["cache_path"])
        if cache_path.exists():
            manifest.at[idx, "computed"] = True
            manifest.at[idx, "computed_at"] = datetime.now(timezone.utc).isoformat()
            cached += 1
            logger.info("Cache hit for %s %s", gene, hgvs)
            continue

        # Run M1-M4
        logger.info("Computing features for %s %s ...", gene, hgvs)
        try:
            record = _run_m1_through_m4(gene, hgvs, skip_m4=skip_m4)
            if record is None:
                raise RuntimeError("Pipeline returned None")

            # Save VariantRecord JSON
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            record.save(str(cache_path))

            manifest.at[idx, "computed"] = True
            manifest.at[idx, "computed_at"] = datetime.now(timezone.utc).isoformat()
            manifest.at[idx, "failed_reason"] = None
            computed += 1
            gene_failures[gene] = 0  # Reset consecutive failures

            logger.info(
                "Computed %s %s: %d features available",
                gene, hgvs, record.count_available_features(),
            )
        except Exception as e:
            manifest.at[idx, "computed"] = False
            manifest.at[idx, "failed_reason"] = str(e)
            failed += 1

            gene_failures[gene] = gene_failures.get(gene, 0) + 1
            if gene_failures[gene] >= MAX_CONSECUTIVE_FAILURES:
                skipped_genes.add(gene)
                logger.warning(
                    "Skipping gene %s after %d consecutive failures",
                    gene, MAX_CONSECUTIVE_FAILURES,
                )

            logger.warning("Failed %s %s: %s", gene, hgvs, e)

        variants_processed += 1

        # Log progress every 10 variants
        if variants_processed % 10 == 0:
            logger.info(
                "Progress: %d processed, %d computed, %d failed, %d cached",
                variants_processed, computed, failed, cached,
            )

        # Save manifest after each variant (enables resume)
        manifest.to_csv(manifest_path, index=False)

        # Rate-limit API calls
        if sleep_seconds > 0:
            time.sleep(sleep_seconds)

    # Final save
    manifest.to_csv(manifest_path, index=False)

    logger.info(
        "Feature computation complete: total=%d, cached=%d, computed=%d, failed=%d",
        total, cached, computed, failed,
    )
    return (total, cached, computed, failed)


# =============================================================================
# PHASE C: TRAINING & EVALUATION
# =============================================================================

def load_cached_records(
    manifest_path: Optional[Path] = None,
) -> tuple[list[VariantRecord], list[int], list[str]]:
    """Load cached VariantRecords and their labels from manifest.

    Only includes variants where computed=True and the cache file exists.

    Args:
        manifest_path: Path to variants_manifest.csv. Uses default if None.

    Returns:
        Tuple of (records, labels, genes) — parallel lists.

    Raises:
        FileNotFoundError: If manifest file does not exist.
    """
    if manifest_path is None:
        manifest_path = TRAINING_CACHE_DIR / "variants_manifest.csv"

    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    manifest = pd.read_csv(manifest_path)
    records: list[VariantRecord] = []
    labels: list[int] = []
    genes: list[str] = []

    computed_mask = manifest["computed"].fillna(False).astype(str).str.lower() == "true"
    computed_df = manifest[computed_mask]

    for _, row in computed_df.iterrows():
        cache_path = Path(row["cache_path"])
        if not cache_path.exists():
            logger.warning("Cache file missing for %s %s: %s", row["gene"], row["hgvs"], cache_path)
            continue
        try:
            record = VariantRecord.load(str(cache_path))
            records.append(record)
            labels.append(int(row["label"]))
            genes.append(row["gene"])
        except Exception as e:
            logger.warning(
                "Failed to load cached record for %s %s: %s",
                row["gene"], row["hgvs"], e,
            )

    logger.info("Loaded %d cached records from %d computed variants", len(records), len(computed_df))
    return records, labels, genes


def train_and_evaluate(
    records: list[VariantRecord],
    labels: list[int],
    genes: list[str],
    output_dir: Optional[Path] = None,
    n_cv_folds: int = 5,
    missingness_rate: float = 0.15,
    model_version: Optional[str] = None,
    seed: int = 42,
) -> dict:
    """Build features, cross-validate, train final model, save everything.

    Applies simulated missingness to a fraction of records, builds the
    feature matrix, runs gene-stratified cross-validation with fresh
    models per fold, then trains a final model on all data.

    Args:
        records: List of VariantRecords with features computed.
        labels: Parallel list of labels (0=benign, 1=pathogenic).
        genes: Parallel list of gene symbols.
        output_dir: Where to save trained models. Defaults to MODELS_DIR.
        n_cv_folds: Number of cross-validation folds.
        missingness_rate: Fraction of records to apply simulated missingness.
        model_version: Version tag for the model. Auto-generated if None.
        seed: Random seed for reproducibility.

    Returns:
        Dict of evaluation metrics including per-fold and aggregate results.
    """
    from varis.m5_scoring.data_loader import build_training_dataset, gene_stratified_split
    from varis.m5_scoring.ensemble import train_ensemble
    from varis.m5_scoring.feature_extractor import simulate_missingness
    from varis.m5_scoring.benchmarks import save_training_manifest

    if output_dir is None:
        output_dir = MODELS_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if model_version is None:
        model_version = datetime.now(timezone.utc).strftime("v%Y.%m")

    rng = np.random.RandomState(seed)

    # Apply simulated missingness to a subset of records
    logger.info(
        "Applying simulated missingness (rate=%.2f) to %d records ...",
        missingness_rate, len(records),
    )
    augmented_records = []
    for record in records:
        if rng.random() < missingness_rate:
            augmented_records.append(simulate_missingness(record, rate=0.5))
        else:
            augmented_records.append(record)

    # Build feature matrix
    logger.info("Building feature matrix from %d records ...", len(augmented_records))
    feature_df = build_training_dataset(augmented_records, labels)

    feature_cols = [c for c in feature_df.columns if c != "label"]
    X = feature_df[feature_cols]
    y = feature_df["label"]

    # Create DataFrame for gene-stratified splitting
    split_df = pd.DataFrame({
        "gene": genes,
        "label": labels,
    })

    # Gene-stratified cross-validation
    logger.info("Running %d-fold gene-stratified cross-validation ...", n_cv_folds)
    fold_metrics: list[dict] = []

    n_unique_genes = split_df["gene"].nunique()
    effective_folds = min(n_cv_folds, n_unique_genes)
    if effective_folds < n_cv_folds:
        logger.warning(
            "Only %d unique genes available, reducing folds from %d to %d",
            n_unique_genes, n_cv_folds, effective_folds,
        )

    for fold_idx, (train_idx, test_idx) in enumerate(
        gene_stratified_split(split_df, n_splits=effective_folds)
    ):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]
        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]

        train_genes = set(split_df.iloc[train_idx]["gene"])
        test_genes = set(split_df.iloc[test_idx]["gene"])

        logger.info(
            "Fold %d: train=%d (genes: %s), test=%d (genes: %s)",
            fold_idx + 1,
            len(train_idx),
            ", ".join(sorted(train_genes)),
            len(test_idx),
            ", ".join(sorted(test_genes)),
        )

        # Train fresh models for this fold
        import tempfile
        fold_dir = Path(tempfile.mkdtemp(prefix=f"varis_fold{fold_idx}_"))
        try:
            train_ensemble(X_train, y_train, output_dir=fold_dir)

            # Load fold models and predict
            from varis.m5_scoring.ensemble import load_ensemble, predict_from_models
            fold_models = load_ensemble(fold_dir)

            # Predict on test set
            y_scores = []
            for i in range(len(X_test)):
                sample = X_test.iloc[i].to_dict()
                try:
                    pred = predict_from_models(fold_models, sample)
                    y_scores.append(pred["score_ensemble"])
                except Exception as e:
                    logger.warning("Prediction failed for test sample %d: %s", i, e)
                    y_scores.append(0.5)

            y_scores_arr = np.array(y_scores)
            y_test_arr = y_test.values

            # Compute metrics
            fold_result = {
                "fold": fold_idx + 1,
                "train_size": len(train_idx),
                "test_size": len(test_idx),
                "train_genes": sorted(train_genes),
                "test_genes": sorted(test_genes),
            }

            # ROC-AUC (requires both classes in test set)
            if len(np.unique(y_test_arr)) > 1:
                fold_result["roc_auc"] = float(roc_auc_score(y_test_arr, y_scores_arr))
                fold_result["pr_auc"] = float(average_precision_score(y_test_arr, y_scores_arr))
            else:
                fold_result["roc_auc"] = None
                fold_result["pr_auc"] = None
                logger.warning("Fold %d has only one class in test set", fold_idx + 1)

            fold_metrics.append(fold_result)
            logger.info(
                "Fold %d: ROC-AUC=%.4f, PR-AUC=%.4f",
                fold_idx + 1,
                fold_result.get("roc_auc") or 0.0,
                fold_result.get("pr_auc") or 0.0,
            )
        except Exception as e:
            logger.warning("Fold %d training failed: %s", fold_idx + 1, e)
            fold_metrics.append({
                "fold": fold_idx + 1,
                "error": str(e),
                "roc_auc": None,
                "pr_auc": None,
            })
        finally:
            # Clean up temp fold directory
            import shutil
            try:
                shutil.rmtree(fold_dir)
            except OSError:
                pass

    # Aggregate CV metrics
    valid_roc = [f["roc_auc"] for f in fold_metrics if f.get("roc_auc") is not None]
    valid_pr = [f["pr_auc"] for f in fold_metrics if f.get("pr_auc") is not None]

    cv_results = {
        "n_folds": effective_folds,
        "folds": fold_metrics,
    }
    if valid_roc:
        cv_results["roc_auc_mean"] = float(np.mean(valid_roc))
        cv_results["roc_auc_std"] = float(np.std(valid_roc))
    if valid_pr:
        cv_results["pr_auc_mean"] = float(np.mean(valid_pr))
        cv_results["pr_auc_std"] = float(np.std(valid_pr))

    logger.info(
        "CV Results: ROC-AUC=%.4f +/- %.4f, PR-AUC=%.4f +/- %.4f",
        cv_results.get("roc_auc_mean", 0.0),
        cv_results.get("roc_auc_std", 0.0),
        cv_results.get("pr_auc_mean", 0.0),
        cv_results.get("pr_auc_std", 0.0),
    )

    # Train final model on ALL data
    logger.info("Training final model on all %d samples ...", len(X))
    train_ensemble(X, y, output_dir=output_dir)

    # Save training manifest
    metrics = {
        "cv_results": cv_results,
        "model_version": model_version,
        "pipeline_version": PIPELINE_VERSION,
        "n_training_samples": len(X),
        "n_features": len(feature_cols),
        "feature_columns": feature_cols,
        "class_distribution": {
            "pathogenic": int((y == 1).sum()),
            "benign": int((y == 0).sum()),
        },
        "missingness_rate": missingness_rate,
        "seed": seed,
    }
    metadata = {
        "model_version": model_version,
        "n_samples": len(X),
        "n_features": len(feature_cols),
        "n_genes": len(set(genes)),
        "genes": sorted(set(genes)),
    }
    save_training_manifest(output_dir, metrics, metadata)

    logger.info("Final model saved to %s (version: %s)", output_dir, model_version)
    return metrics


# =============================================================================
# STATUS COMMAND
# =============================================================================

def show_status(
    manifest_path: Optional[Path] = None,
    show_failed: bool = False,
    by_gene: bool = False,
) -> None:
    """Print progress summary of feature computation.

    Args:
        manifest_path: Path to variants_manifest.csv. Uses default if None.
        show_failed: If True, list failed variants with reasons.
        by_gene: If True, show per-gene breakdown.
    """
    if manifest_path is None:
        manifest_path = TRAINING_CACHE_DIR / "variants_manifest.csv"

    if not manifest_path.exists():
        logger.info("No manifest found at %s. Run 'select' first.", manifest_path)
        return

    manifest = pd.read_csv(manifest_path)
    total = len(manifest)
    computed_mask = manifest["computed"].fillna(False).astype(str).str.lower() == "true"
    computed = computed_mask.sum()
    failed_mask = manifest["failed_reason"].notna() & ~computed_mask
    failed = failed_mask.sum()
    pending = total - computed - failed

    logger.info("=" * 60)
    logger.info("M5 Training Pipeline Status")
    logger.info("=" * 60)
    logger.info("Total variants:   %d", total)
    logger.info("Computed:         %d (%.1f%%)", computed, 100 * computed / max(total, 1))
    logger.info("Pending:          %d (%.1f%%)", pending, 100 * pending / max(total, 1))
    logger.info("Failed:           %d (%.1f%%)", failed, 100 * failed / max(total, 1))
    logger.info(
        "Class balance:    pathogenic=%d, benign=%d",
        (manifest["label"] == 1).sum(),
        (manifest["label"] == 0).sum(),
    )

    if by_gene:
        logger.info("-" * 60)
        logger.info("Per-gene breakdown:")
        for gene in sorted(manifest["gene"].unique()):
            gene_df = manifest[manifest["gene"] == gene]
            gene_computed = (gene_df["computed"].fillna(False).astype(str).str.lower() == "true").sum()
            gene_total = len(gene_df)
            gene_failed = gene_df["failed_reason"].notna().sum() - gene_computed
            gene_failed = max(0, gene_failed)
            logger.info(
                "  %s: %d/%d computed, %d failed",
                gene.ljust(10), gene_computed, gene_total, gene_failed,
            )

    if show_failed and failed > 0:
        logger.info("-" * 60)
        logger.info("Failed variants:")
        failed_df = manifest[failed_mask]
        for _, row in failed_df.iterrows():
            logger.info(
                "  %s %s: %s",
                row["gene"], row["hgvs"], row.get("failed_reason", "unknown"),
            )


# =============================================================================
# CLI MAIN
# =============================================================================

def main() -> None:
    """Training pipeline CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Varis M5 Training Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    subparsers = parser.add_subparsers(dest="command", help="Training phase to run")

    # select subcommand
    sel = subparsers.add_parser(
        "select", help="Phase A: Select training variants from ClinVar",
    )
    sel.add_argument(
        "--genes", type=str, default=None,
        help="Comma-separated gene list (default: 20 cancer genes)",
    )
    sel.add_argument(
        "--per-gene", type=int, default=MAX_PER_CLASS_PER_GENE,
        help="Max variants per class per gene (default: 15)",
    )
    sel.add_argument(
        "--min-stars", type=int, default=1,
        help="Minimum ClinVar review stars (default: 1)",
    )
    sel.add_argument("--seed", type=int, default=42)

    # compute subcommand
    comp = subparsers.add_parser(
        "compute", help="Phase B: Compute features via M1-M4 (resumable)",
    )
    comp.add_argument(
        "--sleep", type=float, default=API_SLEEP_SECONDS,
        help=f"Seconds between API calls (default: {API_SLEEP_SECONDS})",
    )
    comp.add_argument(
        "--max-variants", type=int, default=None,
        help="Max variants to process (default: all)",
    )
    comp.add_argument(
        "--skip-m4", action="store_true",
        help="Skip conservation analysis (faster)",
    )

    # train-only subcommand
    tr = subparsers.add_parser(
        "train-only", help="Phase C: Train from cached features only",
    )
    tr.add_argument(
        "--model-version", type=str, default=None,
        help="Model version tag (default: auto-generated)",
    )
    tr.add_argument("--seed", type=int, default=42)
    tr.add_argument("--cv-folds", type=int, default=5)

    # train subcommand (all phases)
    full = subparsers.add_parser(
        "train", help="Run all phases: select + compute + train",
    )
    full.add_argument(
        "--max-variants", type=int, default=None,
        help="Max variants to compute (default: all)",
    )
    full.add_argument(
        "--skip-m4", action="store_true",
        help="Skip conservation analysis (faster)",
    )
    full.add_argument("--seed", type=int, default=42)

    # status subcommand
    st = subparsers.add_parser("status", help="Show computation progress")
    st.add_argument("--show-failed", action="store_true")
    st.add_argument("--by-gene", action="store_true")

    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable debug logging",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s %(message)s",
        )

    command = args.command

    # Default (no subcommand) = "train" (all phases)
    if command is None:
        command = "train"
        # Re-parse with defaults for the full train command
        args.max_variants = None
        args.skip_m4 = False
        args.seed = 42

    if command == "select":
        genes = args.genes.split(",") if args.genes else None
        select_training_variants(
            target_genes=genes,
            max_per_class=args.per_gene,
            min_stars=args.min_stars,
            seed=args.seed,
        )

    elif command == "compute":
        compute_features(
            sleep_seconds=args.sleep,
            max_variants=args.max_variants,
            skip_m4=args.skip_m4,
        )

    elif command == "train-only":
        records, labels, genes = load_cached_records()
        if not records:
            logger.error("No cached records found. Run 'compute' first.")
            return
        train_and_evaluate(
            records, labels, genes,
            n_cv_folds=args.cv_folds,
            model_version=args.model_version,
            seed=args.seed,
        )

    elif command == "train":
        # Phase A
        logger.info("=" * 60)
        logger.info("PHASE A: Selecting training variants from ClinVar")
        logger.info("=" * 60)
        select_training_variants(seed=args.seed)

        # Phase B
        logger.info("=" * 60)
        logger.info("PHASE B: Computing features via M1-M4")
        logger.info("=" * 60)
        total, cached, computed, failed = compute_features(
            max_variants=args.max_variants,
            skip_m4=args.skip_m4,
        )

        if cached + computed == 0:
            logger.error("No variants computed successfully. Cannot train.")
            return

        # Phase C
        logger.info("=" * 60)
        logger.info("PHASE C: Training ensemble model")
        logger.info("=" * 60)
        records, labels, genes = load_cached_records()
        if not records:
            logger.error("No cached records available for training.")
            return
        train_and_evaluate(
            records, labels, genes,
            seed=args.seed,
        )

    elif command == "status":
        show_status(
            show_failed=args.show_failed,
            by_gene=args.by_gene,
        )

    else:
        parser.print_help()
