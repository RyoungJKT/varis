"""Allow running M7 as: python -m varis.m7_evolution <command>

Subcommands:
    retrain   - Run the auto-retrain loop (select variants, train, evaluate, deploy/reject)
    rollback  - Rollback to the previous model version
    status    - Show current model version and recent deploy history
    log       - Display evolution log events
    scout     - Run tool discovery scan (Loop 2)
    integrate - Attempt integration of a tool (Loop 3)

Examples:
    python -m varis.m7_evolution retrain
    python -m varis.m7_evolution rollback --reason "Regression in BRCA1 predictions"
    python -m varis.m7_evolution status
    python -m varis.m7_evolution log --limit 20 --type DEPLOY
"""
import argparse
import json
import logging
import os
import sys
from pathlib import Path

from varis.config import DATA_DIR, MODELS_DIR

logger = logging.getLogger(__name__)

# Default database URL — same pattern as config.py
DEFAULT_DATABASE_URL = os.getenv("DATABASE_URL", f"sqlite:///{DATA_DIR}/varis.db")


def _setup_logging(verbose: bool) -> None:
    """Configure logging level based on verbosity flag.

    Args:
        verbose: If True, set DEBUG level. Otherwise INFO.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _get_log_db(database_url: str | None = None):
    """Initialize the evolution log database and return session factory.

    Args:
        database_url: SQLAlchemy database URL. Defaults to DEFAULT_DATABASE_URL.

    Returns:
        A sessionmaker factory for the evolution log.
    """
    from varis.m7_evolution.evolution_log import init_evolution_log

    db_url = database_url or DEFAULT_DATABASE_URL
    return init_evolution_log(db_url)


def cmd_retrain(args: argparse.Namespace) -> None:
    """Execute the retrain subcommand.

    Args:
        args: Parsed CLI arguments.
    """
    from varis.m7_evolution.auto_retrain import run_retrain_loop

    log_db = _get_log_db()
    archive_root = Path(args.archive_root) if args.archive_root else MODELS_DIR

    logger.info("Starting retrain loop (archive_root=%s)", archive_root)
    result = run_retrain_loop(archive_root=archive_root, log_db=log_db)

    if result["completed"]:
        logger.info(
            "Retrain completed: decision=%s, reason=%s, duration=%.1fs",
            result["decision"],
            result["reason"],
            result["duration"],
        )
    else:
        logger.warning("Retrain did not complete: %s", result["reason"])

    # Print summary to stdout for scripting
    print(json.dumps(result, indent=2, default=str))


def cmd_rollback(args: argparse.Namespace) -> None:
    """Execute the rollback subcommand.

    Args:
        args: Parsed CLI arguments.
    """
    from varis.m7_evolution.model_archive import rollback

    archive_root = Path(args.archive_root) if args.archive_root else MODELS_DIR
    reason = args.reason

    if not reason:
        print("Error: --reason is required for rollback", file=sys.stderr)
        sys.exit(1)

    logger.info("Rolling back model (reason: %s)", reason)
    success = rollback(reason, archive_root=archive_root)

    if success:
        print(f"Rollback successful: {reason}")
    else:
        print("Rollback failed: no previous version available", file=sys.stderr)
        sys.exit(1)


def cmd_status(args: argparse.Namespace) -> None:
    """Execute the status subcommand.

    Args:
        args: Parsed CLI arguments.
    """
    from varis.m7_evolution.model_archive import get_current_version, list_versions

    archive_root = Path(args.archive_root) if args.archive_root else MODELS_DIR

    current = get_current_version(archive_root)
    versions = list_versions(archive_root)

    print(f"Current production version: {current or '(none)'}")
    print(f"Total archived versions: {len(versions)}")
    print()

    if versions:
        print("Version history:")
        print(f"  {'Version':<14} {'Status':<14} {'Deployed At':<28} {'ROC-AUC'}")
        print(f"  {'-'*14} {'-'*14} {'-'*28} {'-'*8}")
        for v in reversed(versions):
            version = v.get("version", "?")
            status = v.get("status", "?")
            deployed_at = v.get("deployed_at", "")
            metrics = v.get("metrics", {})
            roc_auc = metrics.get("roc_auc")
            roc_str = f"{roc_auc:.4f}" if roc_auc is not None else "N/A"
            print(f"  {version:<14} {status:<14} {deployed_at:<28} {roc_str}")
    else:
        print("No model versions found.")


def cmd_log(args: argparse.Namespace) -> None:
    """Execute the log subcommand.

    Args:
        args: Parsed CLI arguments.
    """
    from varis.m7_evolution.evolution_log import get_log

    log_db = _get_log_db()
    events = get_log(log_db, limit=args.limit, event_type=args.type)

    if not events:
        print("No evolution log events found.")
        return

    print(f"Evolution Log ({len(events)} events):")
    print(f"  {'ID':<6} {'Timestamp':<28} {'Type':<20} {'Version':<14} Details")
    print(f"  {'-'*6} {'-'*28} {'-'*20} {'-'*14} {'-'*30}")

    for event in events:
        eid = event.get("id", "?")
        ts = event.get("timestamp", "")
        etype = event.get("event_type", "?")
        version = event.get("model_version", "")
        details = event.get("details")
        details_str = ""
        if details and isinstance(details, dict):
            # Show a compact summary of details
            parts = []
            for k, v in list(details.items())[:3]:
                parts.append(f"{k}={v}")
            details_str = ", ".join(parts)
            if len(details) > 3:
                details_str += "..."
        print(f"  {eid:<6} {ts:<28} {etype:<20} {version:<14} {details_str}")


def cmd_scout(args: argparse.Namespace) -> None:
    """Execute the scout subcommand.

    Args:
        args: Parsed CLI arguments.
    """
    from varis.m7_evolution.tool_scout import run_scout_loop

    log_db = _get_log_db()
    logger.info("Starting tool scout scan")
    result = run_scout_loop(log_db=log_db)
    print(json.dumps(result, indent=2))


def cmd_integrate(args: argparse.Namespace) -> None:
    """Execute the integrate subcommand.

    Args:
        args: Parsed CLI arguments.
    """
    from varis.m7_evolution.auto_integrator import attempt_integration

    log_db = _get_log_db()
    proposal = {"name": args.tool, "package": args.tool, "source": "manual"}

    logger.info("Attempting integration of %s", args.tool)
    result = attempt_integration(proposal, log_db=log_db)
    print(json.dumps(result, indent=2, default=str))


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser.

    Returns:
        Configured ArgumentParser with all subcommands.
    """
    parser = argparse.ArgumentParser(
        prog="python -m varis.m7_evolution",
        description="M7: Self-Evolution — Model management and evolution log",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable debug logging",
    )
    parser.add_argument(
        "--archive-root",
        type=str,
        default=None,
        help=f"Model archive root directory (default: {MODELS_DIR})",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # retrain
    subparsers.add_parser("retrain", help="Run the auto-retrain loop")

    # rollback
    rollback_parser = subparsers.add_parser("rollback", help="Rollback to previous model version")
    rollback_parser.add_argument(
        "--reason",
        type=str,
        required=True,
        help="Reason for the rollback (required)",
    )

    # status
    subparsers.add_parser("status", help="Show current model version and history")

    # log
    log_parser = subparsers.add_parser("log", help="Display evolution log events")
    log_parser.add_argument(
        "--limit",
        type=int,
        default=20,
        help="Maximum number of events to show (default: 20)",
    )
    log_parser.add_argument(
        "--type",
        type=str,
        default=None,
        help="Filter by event type (e.g. DEPLOY, REJECT, RETRAIN_START)",
    )

    # scout
    subparsers.add_parser("scout", help="Run tool discovery scan")

    # integrate
    integrate_parser = subparsers.add_parser("integrate", help="Attempt integration of a tool")
    integrate_parser.add_argument(
        "--tool",
        type=str,
        required=True,
        help="Package name to attempt installing and benchmarking",
    )

    return parser


def main() -> None:
    """Entry point for the M7 CLI."""
    parser = build_parser()
    args = parser.parse_args()

    _setup_logging(args.verbose)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    commands = {
        "retrain": cmd_retrain,
        "rollback": cmd_rollback,
        "status": cmd_status,
        "log": cmd_log,
        "scout": cmd_scout,
        "integrate": cmd_integrate,
    }

    handler = commands.get(args.command)
    if handler is None:
        parser.print_help()
        sys.exit(1)

    try:
        handler(args)
    except Exception as e:
        logger.error("Command '%s' failed: %s", args.command, e, exc_info=args.verbose)
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
