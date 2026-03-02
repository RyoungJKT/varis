"""Report Generator — Creates PDF clinical investigation reports.

Generates downloadable PDF reports with structural evidence, SHAP explanations,
3D visualization snapshots, and ACMG evidence summary.

Primary: WeasyPrint. Fallback: ReportLab. Worst case: Markdown export.
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def generate_pdf_report(variant_record: VariantRecord, output_path: str) -> str:
    """Generate a PDF clinical investigation report.

    Returns path to generated PDF file.
    """
    pass

def generate_html_report(variant_record: VariantRecord) -> str:
    """Generate HTML report (intermediate step for PDF, also standalone)."""
    pass
