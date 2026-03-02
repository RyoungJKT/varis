"""ESMFold Predictor — Fallback structure prediction via Meta's ESMFold API."""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def predict_esmfold(variant_record: VariantRecord) -> VariantRecord:
    """Predict structure using ESMFold. Needs protein_sequence from M1."""
    pass
