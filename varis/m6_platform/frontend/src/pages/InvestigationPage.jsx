import { useState, useEffect } from "react";
import { useParams, Link } from "react-router-dom";
import { getInvestigation } from "../api/client";
import PredictionBadge from "../components/PredictionBadge";
import ProvenanceFooter from "../components/ProvenanceFooter";
import ReliabilityStrip from "../components/ReliabilityStrip";
import ShapWaterfall from "../components/ShapWaterfall";
import EvidencePanel from "../components/EvidencePanel";
import MolstarViewer from "../components/MolstarViewer";

export default function InvestigationPage() {
  const { variantId } = useParams();
  const [data, setData] = useState(null);
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    setLoading(true);
    getInvestigation(variantId)
      .then(setData)
      .catch((err) => setError(err.message))
      .finally(() => setLoading(false));
  }, [variantId]);

  if (loading) {
    return (
      <div className="max-w-7xl mx-auto px-4 py-16 text-center">
        <p className="text-gray-500">Loading investigation...</p>
      </div>
    );
  }

  if (error) {
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

  // Extract evidence tags from features or provenance
  const evidenceTags = data.features
    ?.filter((f) => f.evidence_tag)
    .map((f) => f.evidence_tag)
    .filter((v, i, a) => a.indexOf(v) === i) || [];

  return (
    <div className="max-w-7xl mx-auto px-4 py-6">
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
      </div>

      {/* Two-column layout */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Left: 3D viewer */}
        <div className="bg-white rounded-lg shadow overflow-hidden">
          <div className="h-96">
            <MolstarViewer structure={data.structure} />
          </div>
        </div>

        {/* Right: Panels */}
        <div className="space-y-4">
          {/* Prediction */}
          <PredictionBadge prediction={data.prediction} />

          {/* Reliability Strip */}
          <ReliabilityStrip structure={data.structure} features={data.features} />

          {/* SHAP Waterfall */}
          <ShapWaterfall
            explanation={data.explanation}
            score={data.prediction?.score}
          />

          {/* Evidence Panel */}
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
