"""AlphaFold Retriever — Downloads predicted structures from AlphaFold DB."""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def retrieve_alphafold(variant_record: VariantRecord) -> VariantRecord:
    """Download AlphaFold structure and set pdb_path. Needs uniprot_id from M1."""
    pass
