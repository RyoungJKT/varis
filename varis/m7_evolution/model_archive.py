"""Model Archive — Version storage with atomic deployment and rollback.

Manages model version lifecycle: archive, deploy, rollback, cleanup.
Every model version is stored with SHA256 hashes for integrity verification,
and deployment uses symlink-based atomic switching.

Directory layout under archive_root:
    archive/
        v2026.01/
            version_metadata.json
            model.cbm
            model.json
            ...
        v2026.02/
            ...
    current -> archive/v2026.02   (symlink to production version)

Priority: Part of M7 governance — build BEFORE enabling auto-retrain.
"""

import hashlib
import json
import logging
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from varis.config import DATA_DIR

logger = logging.getLogger(__name__)

# File extensions considered model files for hashing
MODEL_FILE_EXTENSIONS = {".cbm", ".json", ".txt", ".pkl"}

# Default archive root
_DEFAULT_ARCHIVE_ROOT = DATA_DIR / "models"


def _get_archive_root(archive_root: Path | None) -> Path:
    """Resolve the archive root directory.

    Args:
        archive_root: Explicit root, or None to use default.

    Returns:
        The resolved archive root path.
    """
    return archive_root if archive_root is not None else _DEFAULT_ARCHIVE_ROOT


def _compute_file_hashes(directory: Path) -> dict[str, str]:
    """Compute SHA256 hashes for all model files in a directory.

    Args:
        directory: Path to the directory containing model files.

    Returns:
        Dict mapping filename to "sha256:{hexdigest}" for each model file.
    """
    hashes: dict[str, str] = {}
    if not directory.is_dir():
        return hashes

    for filepath in sorted(directory.iterdir()):
        if filepath.is_file() and filepath.suffix in MODEL_FILE_EXTENSIONS:
            sha256 = hashlib.sha256()
            with open(filepath, "rb") as f:
                for chunk in iter(lambda: f.read(8192), b""):
                    sha256.update(chunk)
            hashes[filepath.name] = f"sha256:{sha256.hexdigest()}"

    return hashes


def _get_reproducibility_info() -> dict[str, Any]:
    """Gather reproducibility information about the current environment.

    Returns:
        Dict with Python version, ML library versions, feature schema version,
        and git commit hash if available.
    """
    info: dict[str, Any] = {
        "python_version": sys.version,
        "platform": platform.platform(),
    }

    # ML library versions — optional imports
    for lib_name in ("catboost", "xgboost", "lightgbm", "shap"):
        try:
            mod = __import__(lib_name)
            info[f"{lib_name}_version"] = getattr(mod, "__version__", "unknown")
        except ImportError:
            info[f"{lib_name}_version"] = "not_installed"

    # Feature schema version from VariantRecord
    try:
        from varis.models.variant_record import RECORD_SCHEMA_VERSION
        info["feature_schema_version"] = RECORD_SCHEMA_VERSION
    except ImportError:
        info["feature_schema_version"] = "unknown"

    # Git commit hash
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            info["git_commit"] = result.stdout.strip()
    except Exception:
        pass

    return info


def _read_metadata(version_dir: Path) -> dict[str, Any] | None:
    """Read version_metadata.json from a version directory.

    Args:
        version_dir: Path to the version directory.

    Returns:
        Parsed metadata dict, or None if not found or invalid.
    """
    meta_path = version_dir / "version_metadata.json"
    if not meta_path.exists():
        return None
    try:
        with open(meta_path, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.warning("Failed to read metadata from %s: %s", meta_path, e)
        return None


def _write_metadata(version_dir: Path, metadata: dict[str, Any]) -> None:
    """Write version_metadata.json to a version directory.

    Args:
        version_dir: Path to the version directory.
        metadata: The metadata dict to write.
    """
    meta_path = version_dir / "version_metadata.json"
    with open(meta_path, "w") as f:
        json.dump(metadata, f, indent=2, default=str)


def _validate_hashes(version_dir: Path, stored_hashes: dict[str, str]) -> bool:
    """Validate that files in a version directory match stored hashes.

    Args:
        version_dir: Path to the version directory.
        stored_hashes: Dict of filename -> "sha256:hexdigest" from metadata.

    Returns:
        True if all hashes match, False if any mismatch or missing file.
    """
    current_hashes = _compute_file_hashes(version_dir)

    for filename, expected_hash in stored_hashes.items():
        actual_hash = current_hashes.get(filename)
        if actual_hash != expected_hash:
            logger.warning(
                "Hash mismatch for %s: expected %s, got %s",
                filename, expected_hash, actual_hash,
            )
            return False

    return True


def _find_production_version(archive_dir: Path) -> str | None:
    """Scan archive directory for the version with status='production'.

    Args:
        archive_dir: Path to the archive directory.

    Returns:
        Version string of the production version, or None.
    """
    if not archive_dir.is_dir():
        return None

    for version_dir in sorted(archive_dir.iterdir()):
        if not version_dir.is_dir():
            continue
        metadata = _read_metadata(version_dir)
        if metadata and metadata.get("status") == "production":
            return metadata.get("version")

    return None


def archive_version(
    source_dir: Path,
    version: str,
    metrics: dict[str, Any],
    archive_root: Path | None = None,
) -> bool:
    """Archive a model version from source_dir to the version archive.

    Copies all model files, computes SHA256 hashes, and writes metadata.

    Args:
        source_dir: Directory containing model files to archive.
        version: Version string (e.g. "v2026.03").
        metrics: Dict of evaluation metrics at archive time.
        archive_root: Root directory for the archive. Defaults to DATA_DIR/models.

    Returns:
        True on success, False on failure.
    """
    try:
        root = _get_archive_root(archive_root)
        archive_dir = root / "archive" / version

        if archive_dir.exists():
            logger.warning("Version %s already exists in archive, overwriting", version)
            shutil.rmtree(archive_dir)

        archive_dir.mkdir(parents=True, exist_ok=True)

        # Copy all files from source_dir
        source_path = Path(source_dir)
        if not source_path.is_dir():
            logger.warning("Source directory does not exist: %s", source_dir)
            return False

        for item in source_path.iterdir():
            if item.is_file():
                shutil.copy2(item, archive_dir / item.name)

        # Compute hashes for model files
        file_hashes = _compute_file_hashes(archive_dir)

        # Write metadata
        metadata = {
            "version": version,
            "status": "archived",
            "archived_at": datetime.now(timezone.utc).isoformat(),
            "metrics": metrics,
            "file_hashes": file_hashes,
            "reproducibility": _get_reproducibility_info(),
        }
        _write_metadata(archive_dir, metadata)

        logger.info("Archived model version %s to %s", version, archive_dir)
        return True

    except Exception as e:
        logger.warning("Failed to archive version %s: %s", version, e)
        return False


def deploy_version(
    version: str,
    archive_root: Path | None = None,
) -> bool:
    """Deploy a version by creating a 'current' symlink to it.

    Validates file hashes, updates metadata status, and atomically
    switches the current symlink. On Windows (no symlinks), copies files instead.

    Args:
        version: Version string to deploy (must be in the archive).
        archive_root: Root directory for the archive. Defaults to DATA_DIR/models.

    Returns:
        True on success, False on failure.
    """
    try:
        root = _get_archive_root(archive_root)
        archive_dir = root / "archive"
        version_dir = archive_dir / version
        current_link = root / "current"

        if not version_dir.is_dir():
            logger.warning("Version %s not found in archive", version)
            return False

        # Read and validate metadata
        metadata = _read_metadata(version_dir)
        if metadata is None:
            logger.warning("No metadata found for version %s", version)
            return False

        # Validate file hashes
        stored_hashes = metadata.get("file_hashes", {})
        if stored_hashes and not _validate_hashes(version_dir, stored_hashes):
            logger.warning("Hash validation failed for version %s, aborting deploy", version)
            return False

        # Find and demote current production version
        previous_production = _find_production_version(archive_dir)
        if previous_production and previous_production != version:
            prev_dir = archive_dir / previous_production
            prev_metadata = _read_metadata(prev_dir)
            if prev_metadata:
                prev_metadata["status"] = "archived"
                _write_metadata(prev_dir, prev_metadata)
                logger.info("Demoted previous production version %s to archived", previous_production)

        # Update target version metadata
        metadata["status"] = "production"
        metadata["deployed_at"] = datetime.now(timezone.utc).isoformat()
        _write_metadata(version_dir, metadata)

        # Create or update symlink (or copy on Windows)
        if sys.platform == "win32":
            # Windows fallback: copy files
            if current_link.exists():
                shutil.rmtree(current_link)
            shutil.copytree(version_dir, current_link)
        else:
            # Atomic symlink swap on Unix
            tmp_link = root / "current.tmp"
            if tmp_link.is_symlink() or tmp_link.exists():
                tmp_link.unlink()
            tmp_link.symlink_to(version_dir)
            tmp_link.rename(current_link)

        logger.info("Deployed model version %s as production", version)
        return True

    except Exception as e:
        logger.warning("Failed to deploy version %s: %s", version, e)
        return False


def rollback(
    reason: str,
    archive_root: Path | None = None,
) -> bool:
    """Rollback from current production to the most recent archived version.

    Finds the current production version, marks it as rolled_back with a reason,
    and switches to the most recent non-rejected archived version.

    Args:
        reason: Human-readable reason for the rollback.
        archive_root: Root directory for the archive. Defaults to DATA_DIR/models.

    Returns:
        True if rollback succeeded, False if no previous version available.
    """
    try:
        root = _get_archive_root(archive_root)
        archive_dir = root / "archive"

        if not archive_dir.is_dir():
            logger.warning("No archive directory found at %s", archive_dir)
            return False

        # Find current production version
        current_version = _find_production_version(archive_dir)
        if current_version is None:
            logger.warning("No production version found to rollback from")
            return False

        # Find most recent archived (non-rejected) version
        candidates: list[tuple[str, dict[str, Any]]] = []
        for version_dir in sorted(archive_dir.iterdir()):
            if not version_dir.is_dir():
                continue
            metadata = _read_metadata(version_dir)
            if metadata is None:
                continue
            ver = metadata.get("version", "")
            status = metadata.get("status", "")
            if ver != current_version and status in ("archived", "rolled_back"):
                candidates.append((ver, metadata))

        if not candidates:
            logger.warning("No previous version available for rollback")
            return False

        # Pick the most recent candidate by version string (sorted descending)
        candidates.sort(key=lambda x: x[0], reverse=True)
        target_version, target_metadata = candidates[0]
        target_dir = archive_dir / target_version

        # Validate target hashes
        stored_hashes = target_metadata.get("file_hashes", {})
        if stored_hashes and not _validate_hashes(target_dir, stored_hashes):
            logger.warning("Hash validation failed for rollback target %s", target_version)
            return False

        # Mark current production as rolled_back
        current_dir = archive_dir / current_version
        current_metadata = _read_metadata(current_dir)
        if current_metadata:
            current_metadata["status"] = "rolled_back"
            current_metadata["rollback_reason"] = reason
            current_metadata["rolled_back_at"] = datetime.now(timezone.utc).isoformat()
            _write_metadata(current_dir, current_metadata)

        # Deploy the target version
        target_metadata["status"] = "production"
        target_metadata["deployed_at"] = datetime.now(timezone.utc).isoformat()
        _write_metadata(target_dir, target_metadata)

        # Swap symlink
        current_link = root / "current"
        if sys.platform == "win32":
            if current_link.exists():
                shutil.rmtree(current_link)
            shutil.copytree(target_dir, current_link)
        else:
            tmp_link = root / "current.tmp"
            if tmp_link.is_symlink() or tmp_link.exists():
                tmp_link.unlink()
            tmp_link.symlink_to(target_dir)
            tmp_link.rename(current_link)

        logger.info(
            "Rolled back from %s to %s (reason: %s)",
            current_version, target_version, reason,
        )
        return True

    except Exception as e:
        logger.warning("Rollback failed: %s", e)
        return False


def get_current_version(archive_root: Path | None = None) -> str | None:
    """Get the currently deployed (production) model version.

    Reads the 'current' symlink target or scans metadata for status='production'.

    Args:
        archive_root: Root directory for the archive. Defaults to DATA_DIR/models.

    Returns:
        Version string of the production version, or None if not set.
    """
    try:
        root = _get_archive_root(archive_root)
        current_link = root / "current"

        # Try reading symlink target first
        if current_link.is_symlink():
            target = current_link.resolve()
            metadata = _read_metadata(target)
            if metadata and metadata.get("status") == "production":
                return metadata.get("version")

        # Try reading as directory (Windows copy case)
        if current_link.is_dir():
            metadata = _read_metadata(current_link)
            if metadata and metadata.get("status") == "production":
                return metadata.get("version")

        # Fallback: scan archive for production status
        archive_dir = root / "archive"
        return _find_production_version(archive_dir)

    except Exception as e:
        logger.warning("Failed to get current version: %s", e)
        return None


def list_versions(archive_root: Path | None = None) -> list[dict[str, Any]]:
    """List all archived model versions with their metadata.

    Args:
        archive_root: Root directory for the archive. Defaults to DATA_DIR/models.

    Returns:
        List of metadata dicts, sorted by version string.
    """
    try:
        root = _get_archive_root(archive_root)
        archive_dir = root / "archive"

        if not archive_dir.is_dir():
            return []

        versions: list[dict[str, Any]] = []
        for version_dir in sorted(archive_dir.iterdir()):
            if not version_dir.is_dir():
                continue
            metadata = _read_metadata(version_dir)
            if metadata is not None:
                versions.append(metadata)

        return versions

    except Exception as e:
        logger.warning("Failed to list versions: %s", e)
        return []


def mark_rejected(
    version: str,
    reason: str,
    archive_root: Path | None = None,
) -> None:
    """Mark a version as rejected, preventing it from being deployed or rolled back to.

    Args:
        version: Version string to reject.
        reason: Human-readable reason for rejection.
        archive_root: Root directory for the archive. Defaults to DATA_DIR/models.
    """
    try:
        root = _get_archive_root(archive_root)
        version_dir = root / "archive" / version

        if not version_dir.is_dir():
            logger.warning("Version %s not found in archive", version)
            return

        metadata = _read_metadata(version_dir)
        if metadata is None:
            logger.warning("No metadata found for version %s", version)
            return

        metadata["status"] = "rejected"
        metadata["rejection_reason"] = reason
        metadata["rejected_at"] = datetime.now(timezone.utc).isoformat()
        _write_metadata(version_dir, metadata)

        logger.info("Marked version %s as rejected: %s", version, reason)

    except Exception as e:
        logger.warning("Failed to mark version %s as rejected: %s", version, e)


def cleanup_old_versions(
    archive_root: Path | None = None,
    max_candidates: int = 5,
    max_months: int = 12,
) -> None:
    """Remove old versions beyond retention limits.

    Removes candidate (archived, non-production) versions beyond max_candidates
    (keeping the newest), and removes archived versions older than max_months
    (but never removes production or rolled_back versions).

    Args:
        archive_root: Root directory for the archive. Defaults to DATA_DIR/models.
        max_candidates: Maximum number of archived candidate versions to keep.
        max_months: Maximum age in months for archived versions.
    """
    try:
        root = _get_archive_root(archive_root)
        archive_dir = root / "archive"

        if not archive_dir.is_dir():
            return

        # Collect all versions with metadata
        all_versions: list[tuple[str, Path, dict[str, Any]]] = []
        for version_dir in sorted(archive_dir.iterdir()):
            if not version_dir.is_dir():
                continue
            metadata = _read_metadata(version_dir)
            if metadata is not None:
                all_versions.append((metadata.get("version", ""), version_dir, metadata))

        # Separate protected versions (production, rolled_back) from candidates
        protected_statuses = {"production", "rolled_back"}
        candidates: list[tuple[str, Path, dict[str, Any]]] = []
        for ver, vdir, meta in all_versions:
            if meta.get("status") not in protected_statuses:
                candidates.append((ver, vdir, meta))

        # Sort candidates by version string (newest first)
        candidates.sort(key=lambda x: x[0], reverse=True)

        # Remove excess candidates beyond max_candidates
        if len(candidates) > max_candidates:
            to_remove = candidates[max_candidates:]
            for ver, vdir, _meta in to_remove:
                logger.info("Removing old version %s (exceeded max_candidates=%d)", ver, max_candidates)
                shutil.rmtree(vdir)

        # Remove archived versions older than max_months
        now = datetime.now(timezone.utc)
        remaining_candidates = candidates[:max_candidates]
        for ver, vdir, meta in remaining_candidates:
            archived_at_str = meta.get("archived_at")
            if archived_at_str:
                try:
                    archived_at = datetime.fromisoformat(archived_at_str)
                    age_days = (now - archived_at).days
                    age_months = age_days / 30.44  # Average days per month
                    if age_months > max_months:
                        logger.info(
                            "Removing version %s (age %.1f months > max %d months)",
                            ver, age_months, max_months,
                        )
                        shutil.rmtree(vdir)
                except (ValueError, TypeError):
                    pass

    except Exception as e:
        logger.warning("Cleanup failed: %s", e)
