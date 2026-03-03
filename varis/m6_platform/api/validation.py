"""Input validation for the VarisDB API.

Two layers:
  1. Security: length caps, charset filtering (prevents injection)
  2. Domain: HGVS format recognition (provides friendly errors)
"""

import re
import logging

logger = logging.getLogger(__name__)

MAX_GENE_LENGTH = 50
MAX_HGVS_LENGTH = 100

# Security: only allow alphanumeric, dots, underscores, hyphens, parentheses, greater/less than
SAFE_GENE_PATTERN = re.compile(r"^[A-Za-z0-9_\-]+$")
SAFE_HGVS_PATTERN = re.compile(r"^[A-Za-z0-9_\.\>\*\(\)\+\-]+$")

# Domain: HGVS protein notation patterns
# Matches: p.Arg1699Trp, p.R1699W, p.Arg1699*, p.(Arg1699Trp)
HGVS_PROTEIN_PATTERN = re.compile(
    r"^p\."                              # p. prefix
    r"\(?"                               # optional opening paren
    r"("
    r"[A-Z][a-z]{2}\d+[A-Z][a-z]{2}"    # 3-letter: Arg1699Trp
    r"|[A-Z]\d+[A-Z\*]"                 # 1-letter: R1699W or R1699*
    r"|[A-Z][a-z]{2}\d+\*"              # 3-letter nonsense: Arg1699*
    r"|[A-Z][a-z]{2}\d+del"             # deletion: Arg1699del
    r"|[A-Z][a-z]{2}\d+_[A-Z][a-z]{2}\d+del"  # range deletion
    r"|[A-Z][a-z]{2}\d+dup"             # duplication
    r")"
    r"\)?$"                              # optional closing paren
)


def validate_variant_input(gene: str, hgvs: str) -> dict:
    """Validate gene and HGVS input for security and domain correctness.

    Args:
        gene: Gene symbol string.
        hgvs: HGVS protein notation string.

    Returns:
        Dict with 'valid' (bool) and optionally 'error' (str).
    """
    # Layer 1: Security — length
    if len(gene) > MAX_GENE_LENGTH:
        return {"valid": False, "error": f"Gene name exceeds maximum length of {MAX_GENE_LENGTH} characters."}
    if len(hgvs) > MAX_HGVS_LENGTH:
        return {"valid": False, "error": f"HGVS notation exceeds maximum length of {MAX_HGVS_LENGTH} characters."}

    # Layer 1: Security — charset
    if not SAFE_GENE_PATTERN.match(gene):
        return {"valid": False, "error": "Gene name contains invalid characters. Only letters, digits, hyphens, and underscores are allowed."}
    if not SAFE_HGVS_PATTERN.match(hgvs):
        return {"valid": False, "error": "HGVS notation contains invalid characters."}

    # Layer 2: Domain — gene not empty
    if len(gene.strip()) == 0:
        return {"valid": False, "error": "Gene name cannot be empty."}

    # Layer 2: Domain — HGVS format
    if not HGVS_PROTEIN_PATTERN.match(hgvs):
        return {
            "valid": False,
            "error": (
                "Expected HGVS protein notation like 'p.Arg1699Trp' or 'p.R1699W'. "
                "See https://varnomen.hgvs.org/recommendations/protein/ for format details."
            ),
        }

    return {"valid": True}
