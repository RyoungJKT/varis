"""Benchmarks — Compare Varis ensemble against existing tools.

Rigorous evaluation against AlphaMissense, PolyPhen-2, SIFT, and CADD
on the same held-out test set. Reports ROC-AUC, PR-AUC, precision, recall, F1.
"""
import logging
logger = logging.getLogger(__name__)

def run_benchmarks(test_features, test_labels, ensemble) -> dict:
    """Run full benchmark comparison.

    Returns dict with metrics for each tool/model and comparison plots.
    """
    pass

def compare_with_alphamissense(varis_scores, am_scores, labels) -> dict:
    """Head-to-head comparison with AlphaMissense on same test set."""
    pass
