"""Clustal Omega Client — Multiple sequence alignment via EBI REST API.

Submits orthologous protein sequences to the EBI Clustal Omega service,
polls for completion, and parses the aligned FASTA output. The aligned
sequences preserve gap characters introduced by the alignment algorithm.

API docs: https://www.ebi.ac.uk/Tools/services/rest/clustalo

Priority: 2 (primary alignment tool)
Fallback: MAFFT, or skip alignment and use ConSurf pre-computed scores.
"""

import logging
import time
from typing import Optional

import httpx

from varis.models.variant_record import VariantRecord

logger = logging.getLogger(__name__)

# =============================================================================
# CONSTANTS
# =============================================================================

_BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
_SUBMIT_URL = f"{_BASE_URL}/run"
_STATUS_URL = f"{_BASE_URL}/status/{{job_id}}"
_RESULT_URL = f"{_BASE_URL}/result/{{job_id}}/aln-fasta"

_EMAIL = "varis@russellgenetics.com"
_TIMEOUT_SECONDS = 60
_DEFAULT_MAX_POLLS = 30
_DEFAULT_POLL_INTERVAL = 5  # seconds


# =============================================================================
# PUBLIC API
# =============================================================================

def run_alignment(
    variant_record: VariantRecord,
    orthologs: Optional[dict],
    client: Optional[httpx.Client] = None,
    max_polls: int = _DEFAULT_MAX_POLLS,
    poll_interval: float = _DEFAULT_POLL_INTERVAL,
) -> tuple[VariantRecord, Optional[dict]]:
    """Align orthologous sequences using EBI Clustal Omega REST API.

    Submits a FASTA-formatted set of sequences, polls until the job
    completes (or times out), and parses the aligned FASTA result back
    into a dict matching the input orthologs structure.

    Args:
        variant_record: The shared VariantRecord (used for logging context).
        orthologs: Dict with keys:
            - "sequences": dict of seq_id -> amino acid sequence
            - "query_id": str identifying the query protein
            - "taxonomy": dict of seq_id -> NCBI taxon ID (int)
            If None, alignment is skipped.
        client: Optional httpx.Client for dependency injection in tests.
            If None, a new client is created with default timeout.
        max_polls: Maximum number of status polls before giving up.
        poll_interval: Seconds to sleep between status polls.

    Returns:
        A tuple of (variant_record, alignment_dict) where alignment_dict
        has the same structure as the input orthologs dict but with
        gap-containing aligned sequences. Returns (variant_record, None)
        if orthologs is None, submission fails, polling times out, or
        any error occurs.
    """
    if orthologs is None:
        logger.warning(
            "No orthologs provided for %s — skipping alignment",
            variant_record.variant_id,
        )
        return variant_record, None

    sequences = orthologs.get("sequences")
    if not sequences or len(sequences) < 2:
        logger.warning(
            "Fewer than 2 sequences for %s — skipping alignment",
            variant_record.variant_id,
        )
        return variant_record, None

    owns_client = client is None
    if owns_client:
        client = httpx.Client(
            timeout=_TIMEOUT_SECONDS,
            follow_redirects=True,
        )

    try:
        # Step 1: Submit job
        job_id = _submit_job(sequences, client)
        if job_id is None:
            return variant_record, None

        # Step 2: Poll for completion
        finished = _poll_status(job_id, client, max_polls, poll_interval)
        if not finished:
            logger.warning(
                "Clustal Omega job %s did not finish within %d polls",
                job_id, max_polls,
            )
            return variant_record, None

        # Step 3: Retrieve aligned FASTA
        aligned_fasta = _fetch_result(job_id, client)
        if aligned_fasta is None:
            return variant_record, None

        # Step 4: Parse aligned FASTA into dict
        aligned_sequences = _parse_aligned_fasta(aligned_fasta)
        if not aligned_sequences:
            logger.warning(
                "Failed to parse aligned FASTA from Clustal Omega job %s",
                job_id,
            )
            return variant_record, None

        # Build result preserving taxonomy and query_id from input
        alignment = {
            "sequences": aligned_sequences,
            "query_id": orthologs.get("query_id", "query"),
            "taxonomy": orthologs.get("taxonomy", {}),
        }

        logger.info(
            "Clustal Omega alignment complete for %s: %d sequences aligned",
            variant_record.variant_id,
            len(aligned_sequences),
        )

        return variant_record, alignment

    except Exception as e:
        logger.warning(
            "Clustal Omega alignment failed for %s: %s",
            variant_record.variant_id,
            e,
        )
        return variant_record, None

    finally:
        if owns_client:
            client.close()


# =============================================================================
# PRIVATE HELPERS
# =============================================================================

def _sequences_to_fasta(sequences: dict[str, str]) -> str:
    """Convert a dict of seq_id -> sequence to FASTA format string.

    Args:
        sequences: Dictionary mapping sequence identifiers to amino acid
            sequences.

    Returns:
        A FASTA-formatted string with one entry per sequence.
    """
    lines = []
    for seq_id, seq in sequences.items():
        lines.append(f">{seq_id}")
        lines.append(seq)
    return "\n".join(lines) + "\n"


def _submit_job(
    sequences: dict[str, str],
    client: httpx.Client,
) -> Optional[str]:
    """Submit sequences to EBI Clustal Omega for alignment.

    Args:
        sequences: Dict of seq_id -> amino acid sequence.
        client: httpx.Client to use for the request.

    Returns:
        The job ID string, or None on failure.
    """
    fasta_str = _sequences_to_fasta(sequences)

    try:
        response = client.post(
            _SUBMIT_URL,
            data={
                "email": _EMAIL,
                "sequence": fasta_str,
                "outfmt": "fa",
            },
        )

        if response.status_code != 200:
            logger.warning(
                "Clustal Omega submit returned HTTP %d: %s",
                response.status_code,
                response.text[:200],
            )
            return None

        job_id = response.text.strip()
        if not job_id:
            logger.warning("Clustal Omega submit returned empty job ID")
            return None

        logger.info("Clustal Omega job submitted: %s", job_id)
        return job_id

    except httpx.TimeoutException:
        logger.warning("Clustal Omega submit request timed out")
        return None
    except httpx.HTTPError as e:
        logger.warning("HTTP error submitting Clustal Omega job: %s", e)
        return None


def _poll_status(
    job_id: str,
    client: httpx.Client,
    max_polls: int,
    poll_interval: float,
) -> bool:
    """Poll EBI for job completion status.

    Args:
        job_id: The Clustal Omega job ID.
        client: httpx.Client to use for status requests.
        max_polls: Maximum number of polls before giving up.
        poll_interval: Seconds to sleep between polls.

    Returns:
        True if job finished successfully, False if timed out or errored.
    """
    url = _STATUS_URL.format(job_id=job_id)

    for poll_num in range(max_polls):
        try:
            response = client.get(url)

            if response.status_code != 200:
                logger.warning(
                    "Clustal Omega status check returned HTTP %d for job %s",
                    response.status_code,
                    job_id,
                )
                return False

            status = response.text.strip()

            if status == "FINISHED":
                logger.info(
                    "Clustal Omega job %s finished after %d polls",
                    job_id, poll_num + 1,
                )
                return True

            if status in ("FAILURE", "ERROR", "NOT_FOUND"):
                logger.warning(
                    "Clustal Omega job %s returned status: %s",
                    job_id, status,
                )
                return False

            # Still running — wait before next poll
            if poll_interval > 0:
                time.sleep(poll_interval)

        except httpx.TimeoutException:
            logger.warning(
                "Clustal Omega status poll %d timed out for job %s",
                poll_num + 1, job_id,
            )
            # Continue polling — transient timeout
            if poll_interval > 0:
                time.sleep(poll_interval)
        except httpx.HTTPError as e:
            logger.warning(
                "HTTP error polling Clustal Omega job %s: %s",
                job_id, e,
            )
            return False

    return False


def _fetch_result(
    job_id: str,
    client: httpx.Client,
) -> Optional[str]:
    """Fetch the aligned FASTA result from a completed Clustal Omega job.

    Args:
        job_id: The Clustal Omega job ID.
        client: httpx.Client to use for the request.

    Returns:
        The aligned FASTA text, or None on failure.
    """
    url = _RESULT_URL.format(job_id=job_id)

    try:
        response = client.get(url)

        if response.status_code != 200:
            logger.warning(
                "Clustal Omega result fetch returned HTTP %d for job %s",
                response.status_code,
                job_id,
            )
            return None

        result_text = response.text
        if not result_text or not result_text.strip():
            logger.warning(
                "Clustal Omega returned empty result for job %s",
                job_id,
            )
            return None

        return result_text

    except httpx.TimeoutException:
        logger.warning(
            "Clustal Omega result fetch timed out for job %s",
            job_id,
        )
        return None
    except httpx.HTTPError as e:
        logger.warning(
            "HTTP error fetching Clustal Omega result for job %s: %s",
            job_id, e,
        )
        return None


def _parse_aligned_fasta(fasta_text: str) -> dict[str, str]:
    """Parse aligned FASTA text into a dict of seq_id -> aligned sequence.

    Aligned sequences contain gap characters ('-') introduced by the
    alignment algorithm. Each sequence ID is taken from the FASTA header
    (first whitespace-delimited token after '>').

    Args:
        fasta_text: The aligned FASTA text from Clustal Omega.

    Returns:
        Dict mapping sequence IDs to aligned sequences (with gaps).
        Empty dict if parsing fails.
    """
    sequences: dict[str, str] = {}
    current_id: Optional[str] = None
    current_lines: list[str] = []

    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            # Save previous entry
            if current_id is not None:
                sequences[current_id] = "".join(current_lines)

            # Parse new header — use first token as ID
            header = line[1:].strip()
            current_id = header.split()[0] if header.split() else header
            current_lines = []
        else:
            current_lines.append(line)

    # Save last entry
    if current_id is not None:
        sequences[current_id] = "".join(current_lines)

    return sequences
