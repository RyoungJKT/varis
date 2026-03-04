import { useState, useEffect } from "react";
import { useParams, Link } from "react-router-dom";
import { getInvestigation, downloadReport, getClinvarSubmission } from "../api/client";
import PredictionBadge from "../components/PredictionBadge";
import ProvenanceFooter from "../components/ProvenanceFooter";
import ReliabilityStrip from "../components/ReliabilityStrip";
import ShapWaterfall from "../components/ShapWaterfall";
import EvidencePanel from "../components/EvidencePanel";
import MolstarViewer from "../components/MolstarViewer";
import ErrorBanner from "../components/ErrorBanner";

export default function InvestigationPage() {
  const { variantId } = useParams();
  const [data, setData] = useState(null);
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(true);
  const [clinvarData, setClinvarData] = useState(null);
  const [clinvarLoading, setClinvarLoading] = useState(false);
  const [showClinvar, setShowClinvar] = useState(false);

  useEffect(() => {
    setLoading(true);
    setError(null);
    getInvestigation(variantId)
      .then(setData)
      .catch((err) => setError(err.message))
      .finally(() => setLoading(false));
  }, [variantId]);

  async function handleDownloadReport() {
    try {
      const resp = await downloadReport(variantId);
      const blob = await resp.blob();
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `${variantId}_report.html`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (err) {
      setError(err.message);
    }
  }

  async function handleClinvarPreview() {
    if (clinvarData) {
      setShowClinvar(!showClinvar);
      return;
    }
    setClinvarLoading(true);
    try {
      const result = await getClinvarSubmission(variantId);
      setClinvarData(result);
      setShowClinvar(true);
    } catch (err) {
      setError(err.message);
    } finally {
      setClinvarLoading(false);
    }
  }

  if (loading) {
    return (
      <div className="max-w-7xl mx-auto px-4 py-16 text-center">
        <p className="text-gray-500">Loading investigation...</p>
      </div>
    );
  }

  if (error && !data) {
    return (
      <div className="max-w-7xl mx-auto px-4 py-16 text-center">
        <div className="bg-red-50 border border-red-200 rounded-lg p-6 inline-block">
          <h2 className="text-lg font-semibold text-red-700 mb-2">Not Found</h2>
          <p className="text-red-600 text-sm mb-4">{error}</p>
          <Link to="/" className="text-blue-600 hover:underline text-sm">Back to search</Link>
        </div>
      </div>
    );
  }

  if (!data) return null;

  const evidenceTags = data.evidence_tags || [];

  return (
    <div className="max-w-7xl mx-auto px-4 py-6">
      {/* Error banner for non-fatal errors */}
      {error && data && (
        <div className="mb-4">
          <ErrorBanner error={error} onDismiss={() => setError(null)} />
        </div>
      )}

      {/* Header */}
      <div className="mb-6">
        <div className="flex items-baseline gap-3">
          <h2 className="text-2xl font-bold text-gray-900">{data.gene}</h2>
          <span className="text-lg text-gray-600">{data.hgvs}</span>
          {data.clinvar_id && (
            <a
              href={`https://www.ncbi.nlm.nih.gov/clinvar/variation/${data.clinvar_id.replace("VCV", "")}`}
              target="_blank"
              rel="noopener noreferrer"
              className="text-sm text-blue-500 hover:underline"
            >
              {data.clinvar_id}
            </a>
          )}
        </div>
        {/* Action buttons */}
        <div className="flex gap-2 mt-3">
          <button
            onClick={handleDownloadReport}
            className="px-3 py-1.5 text-xs font-medium bg-white border border-gray-300 rounded-md hover:bg-gray-50 text-gray-700"
          >
            Download Report
          </button>
          <button
            onClick={handleClinvarPreview}
            disabled={clinvarLoading}
            className="px-3 py-1.5 text-xs font-medium bg-white border border-gray-300 rounded-md hover:bg-gray-50 text-gray-700 disabled:opacity-50"
          >
            {clinvarLoading ? "Loading..." : "ClinVar Preview"}
          </button>
        </div>
      </div>

      {/* ClinVar submission preview panel */}
      {showClinvar && clinvarData && (
        <div className="mb-6 bg-white rounded-lg shadow border">
          <div className="flex items-center justify-between px-4 py-3 border-b bg-gray-50 rounded-t-lg">
            <h3 className="text-sm font-semibold text-gray-700">ClinVar Submission Preview</h3>
            <button
              onClick={() => setShowClinvar(false)}
              className="text-gray-400 hover:text-gray-600 text-lg leading-none"
              aria-label="Close ClinVar preview"
            >
              &times;
            </button>
          </div>
          <div className="p-4">
            {clinvarData.eligible ? (
              <div>
                <span className="inline-block px-2 py-1 text-xs font-medium bg-green-100 text-green-700 rounded-full mb-3">
                  Eligible for submission
                </span>
                <pre className="text-xs bg-gray-50 rounded p-3 overflow-x-auto max-h-64 overflow-y-auto">
                  {JSON.stringify(clinvarData.submission, null, 2)}
                </pre>
              </div>
            ) : (
              <div className="text-sm text-gray-500">
                <span className="inline-block px-2 py-1 text-xs font-medium bg-gray-100 text-gray-600 rounded-full mb-2">
                  Not eligible
                </span>
                <p>This variant does not meet ClinVar submission criteria (requires M5 scoring, non-uncertain classification, and high model agreement).</p>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Two-column layout */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Left: 3D viewer */}
        <div className="bg-white rounded-lg shadow overflow-hidden h-96">
          <MolstarViewer structure={data.structure} />
        </div>

        {/* Right: Panels */}
        <div className="space-y-4">
          <PredictionBadge prediction={data.prediction} />
          <ReliabilityStrip structure={data.structure} features={data.features} />
          <ShapWaterfall
            explanation={data.explanation}
            score={data.prediction?.score}
          />
          <EvidencePanel evidenceTags={evidenceTags} />
        </div>
      </div>

      {/* Provenance Footer */}
      <div className="mt-6">
        <ProvenanceFooter provenance={data.provenance} />
      </div>
    </div>
  );
}
