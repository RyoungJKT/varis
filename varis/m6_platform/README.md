# M6: VarisDB Platform

## Overview

VarisDB is the public-facing web platform for Varis investigations. It serves two purposes:
(1) a searchable database of all investigated variants, and (2) an interactive investigation
viewer that presents structural evidence to clinicians and researchers.

See CLAUDE.md for M6 UI Architecture Rules. See MASTER_PROMPT.md Phase 5 for the full spec.

## Architecture

Frontend (React + Tailwind):
- Left Panel: Mol* 3D Viewer (pdbe-molstar) with mapping-aware residue highlighting
- Right Panel: Reliability Strip + SHAP Waterfall (animated, Recharts) + ACMG Evidence Panel
- Bottom: Data quality metadata, schema version, source database links

Backend (FastAPI):
- GET /api/v1/investigate/{variant_id} -> full investigation JSON payload
- GET /api/v1/variants?q={query} -> search
- PostgreSQL storage of VariantRecords

## API Contract: /api/v1/investigate/{variant_id}

Returns JSON with four sections:
- structure: PDB source, chain, residue_index, pLDDT, coordinate_mapping_confidence
- features: list of { name, value, units, evidence_tag, available }
- prediction: score, classification, confidence bounds, model_agreement, individual_scores
- explanation: list of { feature, value, shap } sorted by |shap| from ACTUAL trained models

CRITICAL: SHAP values computed server-side by real models. Frontend is visualization only.

## Frontend Components

### 1. Mol* 3D Viewer (MolstarViewer.jsx)
- Uses pdbe-molstar (same viewer as PDB and AlphaFold DB)
- confidence "exact"/"high": highlight residue, label "p.Gly123Arg (Chain A: 123)"
- confidence "low": yellow banner "Residue mapping uncertain", dashed highlight
- confidence "failed": red banner, NO residue highlight at all

### 2. Reliability Strip (ReliabilityStrip.jsx)
- Structure confidence: green (pLDDT > 90), yellow (70-90), red (< 70)
- Feature availability: icons per block with null_reason on hover
- Ensemble agreement: green check or orange warning

### 3. SHAP Waterfall (ShapWaterfall.jsx)
- Recharts stacked bar with invisible baselines (upgrade to D3 only if needed)
- Animated: features appear one by one, highest |shap| first
- Each bar: feature name + actual value + SHAP contribution
- Final score at bottom with confidence interval band

### 4. ACMG Evidence Panel (AcmgPanel.jsx)
- Header: "Suggested Evidence Codes" (never "Assigned" or "Determined")
- Each code with supporting measurement and confidence indicator

## Build Order
1. Layer 1: FastAPI + PostgreSQL + investigation endpoint (MVP)
2. Layer 2: React search + variant detail page (text only)
3. Layer 3: SHAP waterfall (Recharts) + Reliability Strip (demo-ready)
4. Layer 4: Mol* 3D viewer with mapping-aware highlighting (impressive)
5. Layer 5: PDF reports (WeasyPrint)
6. Layer 6: Elasticsearch for production search

## Key Rules
- SHAP from trained models only. No mock/placeholder data in charts.
- Never highlight a residue when coordinate_mapping_confidence is "failed".
- Always show Reliability Strip to prevent overclaiming.
- ACMG codes labeled "Suggested evidence" not "Assigned evidence".
- Investigation endpoint reads from the same VariantRecord schema used everywhere.
