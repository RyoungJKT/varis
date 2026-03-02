"""AlphaFold DB Client — Downloads predicted 3D protein structures.

Downloads PDB/mmCIF files from the AlphaFold Protein Structure Database.
This is the primary structure source. ESMFold (in M2) is the fallback.

Populates: structure_source, pdb_path (stored in data/structures/).
"""

import logging
from pathlib import Path
from varis.models.variant_record import VariantRecord
from varis.config import ALPHAFOLD_DB_URL, STRUCTURES_DIR

logger = logging.getLogger(__name__)


def fetch_alphafold_structure(variant_record: VariantRecord) -> VariantRecord:
    """Download AlphaFold predicted structure for the protein.

    Args:
        variant_record: Must have uniprot_id set.

    Returns:
        VariantRecord with pdb_path set to downloaded file, or None if unavailable.
    """
    pass


def _download_pdb(uniprot_id: str, output_dir: Path) -> Path | None:
    """Download PDB file from AlphaFold DB.

    Args:
        uniprot_id: UniProt accession ID.
        output_dir: Directory to save the PDB file.

    Returns:
        Path to downloaded PDB file, or None on failure.
    """
    pass
