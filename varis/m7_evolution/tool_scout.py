"""Tool Scout (Loop 2) — Semi-autonomous tool discovery system.

Monitors PyPI, GitHub, and bioRxiv for new tools and papers that could
improve Varis's pipeline. Uses rule-based scoring to evaluate relevance
and logs discoveries to the evolution log.

Priority: 5 (stretch goal)
Frequency: Daily scan, weekly proposals
Autonomy: Semi-autonomous (AI proposes, system benchmarks or human reviews)
"""
import logging
from datetime import datetime, timedelta, timezone
from typing import Any, Optional

import httpx

from varis.m7_evolution.evolution_log import (
    EVENT_TOOL_DISCOVERY,
    get_log,
    log_event,
)

logger = logging.getLogger(__name__)

# Keywords relevant to Varis's structural variant analysis pipeline
PIPELINE_KEYWORDS: frozenset[str] = frozenset({
    "ddg",
    "delta-g",
    "stability",
    "foldx",
    "rosetta",
    "evoef",
    "conservation",
    "sasa",
    "solvent",
    "dssp",
    "secondary structure",
    "variant effect",
    "missense",
    "pathogenicity",
    "protein structure",
    "alphafold",
    "plddt",
    "contacts",
    "packing",
    "protein stability",
    "structural biology",
    "variant prediction",
    "bioinformatics",
})

# Minimum score for a candidate to be logged as a discovery
SCORE_THRESHOLD: int = 3

# Sources to monitor with search keywords
SCOUT_SOURCES: dict[str, dict[str, Any]] = {
    "pypi": {
        "url": "https://pypi.org/",
        "keywords": ["bioinformatics", "protein structure", "protein stability"],
    },
    "github": {
        "url": "https://github.com/",
        "keywords": ["variant effect", "protein stability", "variant prediction"],
    },
    "biorxiv": {
        "url": "https://www.biorxiv.org/",
        "keywords": ["variant effect prediction", "protein stability"],
    },
}


def score_candidate(candidate: dict) -> int:
    """Score a candidate tool/paper based on pipeline relevance.

    Rule-based scoring:
      - +1 per PIPELINE_KEYWORD found in lowercase(name + description + keywords)
      - +2 if updated field is within the last 90 days
      - +1 if stars > 100

    Args:
        candidate: Dict with keys: name, description, keywords, stars, updated.

    Returns:
        Integer relevance score.
    """
    score = 0

    # Build searchable text from name, description, and keywords
    name = (candidate.get("name") or "").lower()
    description = (candidate.get("description") or "").lower()
    kw_list = candidate.get("keywords") or []
    if isinstance(kw_list, list):
        kw_text = " ".join(str(k).lower() for k in kw_list)
    else:
        kw_text = str(kw_list).lower()

    combined_text = f"{name} {description} {kw_text}"

    # +1 per matching pipeline keyword
    for keyword in PIPELINE_KEYWORDS:
        if keyword in combined_text:
            score += 1

    # +2 if updated within last 90 days
    updated_str = candidate.get("updated")
    if updated_str:
        try:
            # Handle ISO format with or without timezone
            updated_clean = updated_str.replace("Z", "+00:00")
            # Try parsing date-only format first
            try:
                updated_dt = datetime.strptime(
                    updated_clean[:10], "%Y-%m-%d"
                ).replace(tzinfo=timezone.utc)
            except ValueError:
                updated_dt = datetime.fromisoformat(updated_clean)
                if updated_dt.tzinfo is None:
                    updated_dt = updated_dt.replace(tzinfo=timezone.utc)

            cutoff = datetime.now(timezone.utc) - timedelta(days=90)
            if updated_dt >= cutoff:
                score += 2
        except (ValueError, TypeError) as e:
            logger.debug("Could not parse updated date '%s': %s", updated_str, e)

    # +1 if stars > 100
    stars = candidate.get("stars") or 0
    try:
        if int(stars) > 100:
            score += 1
    except (ValueError, TypeError):
        pass

    return score


def deduplicate(
    candidates: list[dict],
    log_db: Optional[Any] = None,
) -> list[dict]:
    """Filter out candidates already discovered in the evolution log.

    Checks candidate names (lowercased) against existing TOOL_DISCOVERY
    events in the evolution log to avoid re-logging the same tool.

    Args:
        candidates: List of candidate dicts with at least a 'name' field.
        log_db: SQLAlchemy session factory for evolution log. If None,
            no deduplication against the log is performed.

    Returns:
        Filtered list of candidates not previously discovered.
    """
    if not log_db:
        return candidates

    try:
        existing_events = get_log(log_db, limit=1000, event_type=EVENT_TOOL_DISCOVERY)
    except Exception as e:
        logger.warning("Failed to retrieve existing discoveries for dedup: %s", e)
        return candidates

    # Extract names from existing discovery events
    existing_names: set[str] = set()
    for event in existing_events:
        details = event.get("details")
        if isinstance(details, dict):
            event_name = details.get("name", "")
            if event_name:
                existing_names.add(event_name.lower())

    # Filter out already-discovered candidates
    novel = []
    for candidate in candidates:
        cand_name = (candidate.get("name") or "").lower()
        if cand_name and cand_name not in existing_names:
            novel.append(candidate)
        elif cand_name in existing_names:
            logger.debug("Dedup: skipping already-discovered '%s'", cand_name)

    return novel


def scan_pypi(query: str) -> list[dict]:
    """Search PyPI for packages matching the query.

    Uses the PyPI search/JSON API to find relevant packages.

    Args:
        query: Search query string.

    Returns:
        List of candidate dicts with keys: name, source, url, description,
        keywords, stars, updated. Returns empty list on failure.
    """
    candidates: list[dict] = []
    try:
        # PyPI doesn't have a search API anymore; use the warehouse search endpoint
        url = "https://pypi.org/search/"
        params = {"q": query}
        response = httpx.get(url, params=params, timeout=30.0, follow_redirects=True)

        if response.status_code != 200:
            logger.warning("PyPI search returned status %d for query '%s'", response.status_code, query)
            return candidates

        # PyPI search returns HTML; parse package names from the results
        # For robustness, also try the JSON API for individual packages
        # Extract package names from search results HTML
        import re
        package_links = re.findall(
            r'/project/([^/]+)/',
            response.text,
        )
        # Deduplicate while preserving order
        seen: set[str] = set()
        unique_packages: list[str] = []
        for pkg in package_links:
            if pkg.lower() not in seen:
                seen.add(pkg.lower())
                unique_packages.append(pkg)

        # Fetch details for top results via JSON API
        for pkg_name in unique_packages[:10]:
            try:
                detail_resp = httpx.get(
                    f"https://pypi.org/pypi/{pkg_name}/json",
                    timeout=15.0,
                    follow_redirects=True,
                )
                if detail_resp.status_code == 200:
                    data = detail_resp.json()
                    info = data.get("info", {})
                    # Extract keywords
                    kw_str = info.get("keywords") or ""
                    kw_list = [k.strip() for k in kw_str.split(",") if k.strip()]

                    candidates.append({
                        "name": info.get("name", pkg_name),
                        "source": "pypi",
                        "url": info.get("project_url", f"https://pypi.org/project/{pkg_name}/"),
                        "description": info.get("summary", ""),
                        "keywords": kw_list,
                        "stars": 0,
                        "updated": (info.get("version") and _get_pypi_upload_date(data)) or "",
                    })
            except Exception as e:
                logger.debug("Failed to fetch PyPI details for '%s': %s", pkg_name, e)

    except Exception as e:
        logger.warning("PyPI scan failed for query '%s': %s", query, e)

    return candidates


def _get_pypi_upload_date(data: dict) -> str:
    """Extract the most recent upload date from PyPI package data.

    Args:
        data: Full PyPI JSON API response.

    Returns:
        ISO date string (YYYY-MM-DD) or empty string if unavailable.
    """
    try:
        releases = data.get("releases", {})
        info = data.get("info", {})
        version = info.get("version", "")
        if version and version in releases:
            version_files = releases[version]
            if version_files:
                upload_time = version_files[0].get("upload_time_iso_8601", "")
                if upload_time:
                    return upload_time[:10]
        # Fallback: check urls
        urls = data.get("urls", [])
        if urls:
            upload_time = urls[0].get("upload_time_iso_8601", "")
            if upload_time:
                return upload_time[:10]
    except Exception as e:
        logger.debug("Failed to extract PyPI upload date: %s", e)
    return ""


def scan_github(query: str) -> list[dict]:
    """Search GitHub repositories for tools matching the query.

    Uses the GitHub search repositories API.

    Args:
        query: Search query string.

    Returns:
        List of candidate dicts with keys: name, source, url, description,
        keywords, stars, updated. Returns empty list on failure.
    """
    candidates: list[dict] = []
    try:
        url = "https://api.github.com/search/repositories"
        params = {
            "q": query,
            "sort": "stars",
            "order": "desc",
            "per_page": 15,
        }
        headers = {
            "Accept": "application/vnd.github.v3+json",
        }
        response = httpx.get(
            url, params=params, headers=headers, timeout=30.0, follow_redirects=True,
        )

        if response.status_code != 200:
            logger.warning(
                "GitHub search returned status %d for query '%s'",
                response.status_code, query,
            )
            return candidates

        data = response.json()
        items = data.get("items", [])

        for item in items:
            candidates.append({
                "name": item.get("name", ""),
                "source": "github",
                "url": item.get("html_url", ""),
                "description": item.get("description") or "",
                "keywords": item.get("topics") or [],
                "stars": item.get("stargazers_count", 0),
                "updated": (item.get("updated_at") or "")[:10],
            })

    except Exception as e:
        logger.warning("GitHub scan failed for query '%s': %s", query, e)

    return candidates


def scan_biorxiv(query: str) -> list[dict]:
    """Search bioRxiv for recent preprints matching the query.

    Uses the bioRxiv content API to find papers published in the
    last 90 days, filtering by query keywords in title/abstract.

    Args:
        query: Search query string.

    Returns:
        List of candidate dicts with keys: name, source, url, description,
        keywords, stars, updated. Returns empty list on failure.
    """
    candidates: list[dict] = []
    try:
        # bioRxiv API: content detail endpoint for recent papers
        end_date = datetime.now(timezone.utc)
        start_date = end_date - timedelta(days=90)
        start_str = start_date.strftime("%Y-%m-%d")
        end_str = end_date.strftime("%Y-%m-%d")

        # bioRxiv API uses date ranges; search bioinformatics server
        url = (
            f"https://api.biorxiv.org/details/biorxiv/{start_str}/{end_str}/0/25"
        )
        response = httpx.get(url, timeout=30.0, follow_redirects=True)

        if response.status_code != 200:
            logger.warning(
                "bioRxiv API returned status %d for query '%s'",
                response.status_code, query,
            )
            return candidates

        data = response.json()
        papers = data.get("collection", [])

        # Filter by query keywords in title or abstract
        query_terms = query.lower().split()

        for paper in papers:
            title = (paper.get("title") or "").lower()
            abstract = (paper.get("abstract") or "").lower()
            combined = f"{title} {abstract}"

            # Check if any query term matches
            if any(term in combined for term in query_terms):
                doi = paper.get("doi", "")
                candidates.append({
                    "name": paper.get("title", "untitled"),
                    "source": "biorxiv",
                    "url": f"https://doi.org/{doi}" if doi else "",
                    "description": (paper.get("abstract") or "")[:500],
                    "keywords": paper.get("category", "").split(";") if paper.get("category") else [],
                    "stars": 0,
                    "updated": paper.get("date", ""),
                })

    except Exception as e:
        logger.warning("bioRxiv scan failed for query '%s': %s", query, e)

    return candidates


def scan_sources() -> list[dict]:
    """Scan all configured sources for potentially relevant new tools/papers.

    Calls scan_pypi, scan_github, and scan_biorxiv with their configured
    keywords and aggregates the results.

    Returns:
        Aggregated list of candidate dicts from all sources.
    """
    all_candidates: list[dict] = []

    # Scan PyPI
    pypi_config = SCOUT_SOURCES.get("pypi", {})
    for keyword in pypi_config.get("keywords", []):
        try:
            results = scan_pypi(keyword)
            all_candidates.extend(results)
        except Exception as e:
            logger.warning("PyPI scan for '%s' failed: %s", keyword, e)

    # Scan GitHub
    github_config = SCOUT_SOURCES.get("github", {})
    for keyword in github_config.get("keywords", []):
        try:
            results = scan_github(keyword)
            all_candidates.extend(results)
        except Exception as e:
            logger.warning("GitHub scan for '%s' failed: %s", keyword, e)

    # Scan bioRxiv
    biorxiv_config = SCOUT_SOURCES.get("biorxiv", {})
    for keyword in biorxiv_config.get("keywords", []):
        try:
            results = scan_biorxiv(keyword)
            all_candidates.extend(results)
        except Exception as e:
            logger.warning("bioRxiv scan for '%s' failed: %s", keyword, e)

    logger.info("Scanned all sources, found %d total candidates", len(all_candidates))
    return all_candidates


def run_scout_loop(log_db: Optional[Any] = None) -> dict:
    """Execute the full tool discovery loop.

    Orchestrates: scan_sources -> score each -> filter by threshold ->
    deduplicate -> log TOOL_DISCOVERY events.

    Args:
        log_db: SQLAlchemy session factory for evolution log. If None,
            discoveries are not logged persistently.

    Returns:
        Dict with keys: scanned (int), logged (int), candidates (list of dicts).
    """
    result: dict[str, Any] = {"scanned": 0, "logged": 0, "candidates": []}

    try:
        # Step 1: Scan all sources
        candidates = scan_sources()
        result["scanned"] = len(candidates)

        # Step 2: Score each candidate
        scored: list[tuple[dict, int]] = []
        for candidate in candidates:
            try:
                s = score_candidate(candidate)
                scored.append((candidate, s))
            except Exception as e:
                logger.warning("Failed to score candidate '%s': %s", candidate.get("name"), e)

        # Step 3: Filter by threshold
        above_threshold = [
            (cand, s) for cand, s in scored if s >= SCORE_THRESHOLD
        ]
        logger.info(
            "%d of %d candidates scored above threshold (%d)",
            len(above_threshold), len(scored), SCORE_THRESHOLD,
        )

        # Step 4: Deduplicate against existing discoveries
        threshold_candidates = [cand for cand, _ in above_threshold]
        novel_candidates = deduplicate(threshold_candidates, log_db=log_db)

        # Step 5: Log TOOL_DISCOVERY events for novel candidates
        logged_count = 0
        for candidate in novel_candidates:
            if log_db:
                try:
                    # Find the score for this candidate
                    cand_score = 0
                    for cand, s in above_threshold:
                        if cand is candidate:
                            cand_score = s
                            break

                    log_event(
                        log_db,
                        EVENT_TOOL_DISCOVERY,
                        details={
                            "name": candidate.get("name", ""),
                            "source": candidate.get("source", ""),
                            "url": candidate.get("url", ""),
                            "description": candidate.get("description", ""),
                            "score": cand_score,
                        },
                    )
                    logged_count += 1
                    logger.info(
                        "Logged TOOL_DISCOVERY for '%s' (score=%d)",
                        candidate.get("name"), cand_score,
                    )
                except Exception as e:
                    logger.warning(
                        "Failed to log discovery for '%s': %s",
                        candidate.get("name"), e,
                    )
            else:
                logged_count += 1

        result["logged"] = logged_count
        result["candidates"] = novel_candidates

    except Exception as e:
        logger.error("Scout loop failed: %s", e)

    return result
