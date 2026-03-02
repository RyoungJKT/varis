"""Command-line interface for Varis.

Usage:
    python -m varis "BRCA1" "p.Arg1699Trp"
    python -m varis "BRCA1" "p.Arg1699Trp" --output results/
    python -m varis --validate  (run validation variants)
"""

import argparse
import logging
import json
import sys
from pathlib import Path


def main():
    """CLI entry point for Varis."""
    parser = argparse.ArgumentParser(
        description="Varis — Structural Evidence for Every Variant",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example: python -m varis BRCA1 p.Arg1699Trp",
    )
    parser.add_argument("gene", nargs="?", help="Gene symbol (e.g., BRCA1)")
    parser.add_argument("hgvs", nargs="?", help="HGVS protein notation (e.g., p.Arg1699Trp)")
    parser.add_argument("--output", "-o", default=".", help="Output directory for JSON results")
    parser.add_argument("--validate", action="store_true", help="Run validation variants")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable debug logging")

    args = parser.parse_args()

    # Configure logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    from varis.pipeline import investigate
    from varis.config import VALIDATION_VARIANTS

    if args.validate:
        _run_validation(investigate, VALIDATION_VARIANTS, args.output)
    elif args.gene and args.hgvs:
        _run_single(investigate, args.gene, args.hgvs, args.output)
    else:
        parser.print_help()
        sys.exit(1)


def _run_single(investigate_fn, gene: str, hgvs: str, output_dir: str):
    """Investigate a single variant and save results."""
    record = investigate_fn(gene, hgvs)
    output_path = Path(output_dir) / f"{record.variant_id}.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    record.save(str(output_path))
    print(f"\nInvestigation complete: {output_path}")
    print(f"Features available: {record.count_available_features()}/15")
    print(f"Ensemble score: {record.score_ensemble}")
    print(f"Classification: {record.classification}")


def _run_validation(investigate_fn, variants: list, output_dir: str):
    """Run all validation variants and report results."""
    print(f"\nRunning {len(variants)} validation variants...\n")
    for v in variants:
        record = investigate_fn(v["gene"], v["hgvs"])
        status = "✓" if record.classification == v["expected"] else "✗"
        print(f"  {status} {v['gene']} {v['hgvs']}: {record.classification} (expected: {v['expected']})")


if __name__ == "__main__":
    main()
