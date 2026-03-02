"""Structure Utilities — pLDDT extraction and coordinate helpers."""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def extract_plddt(variant_record: VariantRecord) -> VariantRecord:
    """Extract pLDDT at mutation site and mean across protein from B-factor column."""
    pass

def get_residue_coordinates(pdb_path: str, position: int) -> dict | None:
    """Extract CA atom {x, y, z} coordinates for a residue."""
    pass
