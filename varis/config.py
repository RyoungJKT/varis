"""Varis configuration — environment variables, constants, and paths.

All API keys and configuration are loaded from environment variables.
See .env.example for the required variables.
"""

import os
from pathlib import Path

# =========================================================================
# PATHS
# =========================================================================
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
VALIDATION_DIR = DATA_DIR / "validation"
STRUCTURES_DIR = DATA_DIR / "structures"
MODELS_DIR = DATA_DIR / "models"
LOGS_DIR = DATA_DIR / "logs"

# =========================================================================
# API ENDPOINTS
# =========================================================================
CLINVAR_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
UNIPROT_BASE_URL = "https://rest.uniprot.org"
ALPHAFOLD_DB_URL = "https://alphafold.ebi.ac.uk/api"
ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb"
ALPHAMISSENSE_URL = "https://zenodo.org/records/8208688"
NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
INTERPRO_API_URL = "https://www.ebi.ac.uk/interpro/api"
CONSURF_API_URL = "https://consurf.tau.ac.il/results"
CLUSTAL_OMEGA_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"

# =========================================================================
# API KEYS (from environment variables)
# =========================================================================
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")  # Optional but increases rate limit
ALPHAFOLD_API_KEY = os.getenv("ALPHAFOLD_API_KEY", "")  # Not currently required

# =========================================================================
# PIPELINE SETTINGS
# =========================================================================
PIPELINE_VERSION = "v1.0"
DDG_DAMAGING_THRESHOLD = 2.0  # kcal/mol — above this is destabilizing
PLDDT_CONFIDENCE_THRESHOLD = 70.0  # Below this, structural predictions are uncertain
CONSERVATION_HIGH_THRESHOLD = 0.9  # Above this, position is highly conserved
GNOMAD_RARE_THRESHOLD = 0.0001  # Below this, variant is rare

# =========================================================================
# ML MODEL SETTINGS
# =========================================================================
ENSEMBLE_MODELS = ["catboost", "xgboost", "lightgbm"]
CATBOOST_PARAMS = {
    "iterations": 1000,
    "learning_rate": 0.05,
    "depth": 6,
    "eval_metric": "AUC",
    "random_seed": 42,
    "verbose": 0,
}
XGBOOST_PARAMS = {
    "n_estimators": 1000,
    "learning_rate": 0.05,
    "max_depth": 6,
    "eval_metric": "auc",
    "random_state": 42,
    "verbosity": 0,
}
LIGHTGBM_PARAMS = {
    "n_estimators": 1000,
    "learning_rate": 0.05,
    "max_depth": 6,
    "metric": "auc",
    "random_state": 42,
    "verbose": -1,
}

# =========================================================================
# VALIDATION VARIANTS — Well-studied variants for testing
# =========================================================================
VALIDATION_VARIANTS = [
    {"gene": "BRCA1", "hgvs": "p.Arg1699Trp", "expected": "likely_pathogenic"},
    {"gene": "BRCA2", "hgvs": "p.Asp2723His", "expected": "likely_pathogenic"},
    {"gene": "TP53", "hgvs": "p.Arg175His", "expected": "likely_pathogenic"},
    {"gene": "CFTR", "hgvs": "p.Gly551Asp", "expected": "likely_pathogenic"},
    {"gene": "BRCA1", "hgvs": "p.Lys1183Arg", "expected": "likely_benign"},
]

# =========================================================================
# AMINO ACID PROPERTIES — For charge change calculation
# =========================================================================
AA_CHARGE = {
    "Arg": "+", "Lys": "+", "His": "+",  # Positive
    "Asp": "-", "Glu": "-",              # Negative
    "Ala": "0", "Cys": "0", "Phe": "0", "Gly": "0", "Ile": "0",
    "Leu": "0", "Met": "0", "Asn": "0", "Pro": "0", "Gln": "0",
    "Ser": "0", "Thr": "0", "Val": "0", "Trp": "0", "Tyr": "0",
}

AA_THREE_TO_ONE = {
    "Ala": "A", "Cys": "C", "Asp": "D", "Glu": "E", "Phe": "F",
    "Gly": "G", "His": "H", "Ile": "I", "Lys": "K", "Leu": "L",
    "Met": "M", "Asn": "N", "Pro": "P", "Gln": "Q", "Arg": "R",
    "Ser": "S", "Thr": "T", "Val": "V", "Trp": "W", "Tyr": "Y",
}

AA_ONE_TO_THREE = {v: k for k, v in AA_THREE_TO_ONE.items()}

# =========================================================================
# M6 PLATFORM SETTINGS — VarisDB web platform
# =========================================================================
DATABASE_URL = os.getenv("DATABASE_URL", "sqlite:///data/varis.db")
CORS_ORIGINS = os.getenv("CORS_ORIGINS", "http://localhost:5173").split(",")
MAX_INVESTIGATION_WORKERS = int(os.getenv("MAX_INVESTIGATION_WORKERS", "2"))
JOB_TIMEOUT_SECONDS = int(os.getenv("JOB_TIMEOUT_SECONDS", "600"))
RATE_LIMIT_PER_MINUTE = int(os.getenv("RATE_LIMIT_PER_MINUTE", "10"))
