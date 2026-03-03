"""Shared test fixtures for the Varis test suite."""
import pytest
from varis.models.variant_record import VariantRecord, create_variant_record

@pytest.fixture
def empty_record():
    """A fresh VariantRecord with just gene and HGVS set."""
    return create_variant_record("BRCA1", "p.Arg1699Trp")

@pytest.fixture
def m1_completed_record():
    """VariantRecord with M1 fields populated."""
    record = create_variant_record("BRCA1", "p.Arg1699Trp")
    record.residue_position = 1699
    record.ref_amino_acid = "Arg"
    record.alt_amino_acid = "Trp"
    record.ref_aa_single = "R"
    record.alt_aa_single = "W"
    record.charge_change = "+ve → neutral"
    record.clinvar_id = "VCV000055361"
    record.clinvar_classification = "Uncertain significance"
    record.uniprot_id = "P38398"
    record.protein_name = "Breast cancer type 1 susceptibility protein"
    record.protein_length = 1863
    record.gnomad_frequency = 0.00003
    record.alphamissense_score = 0.934
    record.alphamissense_class = "likely_pathogenic"
    record.mark_module_completed("M1")
    return record

@pytest.fixture
def fully_populated_record():
    """VariantRecord with all modules completed — for M5 and report testing."""
    record = create_variant_record("BRCA1", "p.Arg1699Trp")
    record.residue_position = 1699
    record.ref_amino_acid = "Arg"
    record.alt_amino_acid = "Trp"
    record.ref_aa_single = "R"
    record.alt_aa_single = "W"
    record.charge_change = "+ve → neutral"
    record.clinvar_id = "VCV000055361"
    record.uniprot_id = "P38398"
    record.gnomad_frequency = 0.00003
    record.alphamissense_score = 0.934
    record.structure_source = "alphafold"
    record.pdb_path = "data/structures/AF-P38398-F1-model_v4.pdb"
    record.pdb_hash = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    record.mutation_site_present = True
    record.mutation_site_plddt = 92.4
    record.plddt_mean = 85.1
    record.plddt_available = True
    record.mutation_site_confidence_bucket = "very_high"
    record.numbering_scheme = "uniprot"
    record.structure_quality_summary = "High confidence AlphaFold structure"
    record.preparation_steps = ["download", "validate_residue"]
    record.ddg_foldx = 3.7
    record.ddg_pyrosetta = 4.1
    record.ddg_mean = 3.9
    record.solvent_accessibility_relative = 0.03
    record.burial_category = "core"
    record.secondary_structure = "H"
    record.secondary_structure_name = "helix"
    record.contacts_wt = 12
    record.hbonds_wt = 3
    record.packing_density = 0.72
    record.domain_name = "BRCT"
    record.domain_id = "PF00533"
    record.domain_criticality = "critical"
    record.domain_start = 1646
    record.domain_end = 1736
    record.conservation_score = 1.0
    record.conservation_method = "clustal_omega"
    record.num_orthologs = 45
    record.conserved_across_mammals = True
    record.msa_num_sequences = 46
    record.msa_gap_fraction_at_site = 0.02
    record.msa_column_index = 1699
    record.score_ensemble = 0.91
    record.score_catboost = 0.93
    record.score_xgboost = 0.90
    record.score_lightgbm = 0.89
    record.confidence_lower = 0.85
    record.confidence_upper = 0.96
    record.classification = "likely_pathogenic"
    record.model_agreement = "high"
    record.features_used = 15
    record.ensemble_version = "v1.2"
    record.acmg_codes = ["PM1", "PP3", "PS3-proxy"]
    record.modules_completed = ["M1", "M2", "M3", "M4", "M5"]
    return record
