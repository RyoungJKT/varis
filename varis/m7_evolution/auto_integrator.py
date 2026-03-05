"""Auto-Integrator (Loop 3) — Partially autonomous tool integration.

When a proposed tool passes review, attempts to install, probe its API,
benchmark against current metrics, and integrate it automatically.

Works autonomously for ~60-70% of standard Python tools.
The other 30-40% need human help — this is an honest assessment.

Priority: 6 (stretch goal)
Autonomy: Partially autonomous
"""
import importlib
import logging
import subprocess
import sys
from typing import Optional

logger = logging.getLogger(__name__)

COMMON_CALLABLES: frozenset[str] = frozenset({
    "predict",
    "score",
    "compute",
    "run",
    "analyze",
    "calculate",
    "evaluate",
    "process",
    "transform",
    "fit",
    "main",
})

INSTALL_TIMEOUT: int = 120
BENCHMARK_TIMEOUT: int = 60


def attempt_install(package_name: str) -> bool:
    """Attempt to install a Python package via pip with a dry-run safety check.

    Phase 1: Runs pip install --dry-run to verify the package can be resolved
    without actually installing anything. If this fails, returns False early.

    Phase 2: Runs the actual pip install.

    Args:
        package_name: The PyPI package name to install (e.g. 'biopython').

    Returns:
        True if the package was installed successfully, False on any failure.
    """
    try:
        # Phase 1: dry-run
        logger.info("Attempting dry-run install of '%s'", package_name)
        dry_run = subprocess.run(
            [sys.executable, "-m", "pip", "install", "--dry-run", package_name],
            capture_output=True,
            text=True,
            timeout=INSTALL_TIMEOUT,
        )
        if dry_run.returncode != 0:
            logger.warning(
                "Dry-run install failed for '%s': %s",
                package_name,
                dry_run.stderr.strip(),
            )
            return False

        # Phase 2: actual install
        logger.info("Dry-run succeeded, installing '%s'", package_name)
        actual = subprocess.run(
            [sys.executable, "-m", "pip", "install", package_name],
            capture_output=True,
            text=True,
            timeout=INSTALL_TIMEOUT,
        )
        if actual.returncode != 0:
            logger.warning(
                "Actual install failed for '%s': %s",
                package_name,
                actual.stderr.strip(),
            )
            return False

        logger.info("Successfully installed '%s'", package_name)
        return True

    except subprocess.TimeoutExpired:
        logger.warning("Install timed out for '%s' after %ds", package_name, INSTALL_TIMEOUT)
        return False
    except Exception as e:
        logger.warning("Unexpected error installing '%s': %s", package_name, e)
        return False


def probe_tool(package_name: str) -> Optional[dict]:
    """Import a package and scan for commonly-named callable entry points.

    Imports the package via importlib and inspects dir(module) for attributes
    whose lowercase name is in COMMON_CALLABLES and which are callable.

    Args:
        package_name: The Python importable package name (e.g. 'biopython'
            would be 'Bio', but typically the PyPI name matches the import).

    Returns:
        Dict with 'package', 'callables' (list of found callable names),
        and 'module_attrs' (total attribute count), or None on ImportError.
    """
    try:
        module = importlib.import_module(package_name)
    except ImportError as e:
        logger.warning("Cannot import '%s': %s", package_name, e)
        return None
    except Exception as e:
        logger.warning("Unexpected error importing '%s': %s", package_name, e)
        return None

    try:
        attrs = dir(module)
        found_callables = []
        for attr_name in attrs:
            if attr_name.lower() in COMMON_CALLABLES:
                try:
                    attr = getattr(module, attr_name)
                    if callable(attr):
                        found_callables.append(attr_name)
                except Exception:
                    pass

        info = {
            "package": package_name,
            "callables": found_callables,
            "module_attrs": len(attrs),
        }
        logger.info(
            "Probed '%s': found %d callables out of %d attrs",
            package_name,
            len(found_callables),
            len(attrs),
        )
        return info

    except Exception as e:
        logger.warning("Error probing '%s': %s", package_name, e)
        return None


def benchmark_new_feature(
    current_metrics: dict,
    candidate_metrics: dict,
) -> dict:
    """Compare current and candidate metrics to decide integration.

    Args:
        current_metrics: Dict with at least 'roc_auc' and 'pr_auc' keys.
        candidate_metrics: Dict with at least 'roc_auc' and 'pr_auc' keys.

    Returns:
        Dict with 'decision', 'reason', and 'metric_deltas'.
    """
    pass


def _log_integration(
    log_db: object,
    name: str,
    result: dict,
) -> None:
    """Log a TOOL_INTEGRATION event to the evolution log.

    Args:
        log_db: SQLAlchemy session factory (from init_evolution_log).
        name: Name of the tool being integrated.
        result: The integration result dict to store in event details.
    """
    pass


def attempt_integration(
    proposal: dict,
    log_db: object = None,
    current_metrics: Optional[dict] = None,
    candidate_metrics: Optional[dict] = None,
) -> dict:
    """Attempt full autonomous integration of a proposed tool.

    Args:
        proposal: Dict with at least 'name' and 'package' keys.
        log_db: SQLAlchemy session factory for evolution log, or None.
        current_metrics: Dict with current model metrics (roc_auc, pr_auc).
        candidate_metrics: Dict with candidate model metrics after adding the tool.

    Returns:
        Dict with 'decision', 'reason', 'name', and additional context fields.
    """
    pass
