const EVIDENCE_DESCRIPTIONS = {
  computational_support: {
    label: "Computational Support",
    acmg: "PP3-like",
    description: "Multiple computational tools agree on pathogenicity prediction.",
  },
  rarity_evidence: {
    label: "Population Rarity",
    acmg: "PM2-like",
    description: "Variant is absent or extremely rare in population databases (gnomAD).",
  },
  energetics_support: {
    label: "Energetic Destabilization",
    acmg: "Structural",
    description: "Predicted change in protein stability (\u0394\u0394G) exceeds damaging threshold.",
  },
  domain_context: {
    label: "Critical Domain",
    acmg: "PM1-like",
    description: "Variant located in a functionally critical protein domain.",
  },
};

export default function EvidencePanel({ evidenceTags }) {
  if (!evidenceTags || evidenceTags.length === 0) {
    return null;
  }

  return (
    <div className="bg-white rounded-lg shadow p-4">
      <h3 className="text-sm font-semibold text-gray-700 mb-3">
        Suggested Computational Evidence
      </h3>
      <div className="space-y-2">
        {evidenceTags.map((tag) => {
          const info = EVIDENCE_DESCRIPTIONS[tag] || {
            label: tag.replace(/_/g, " "),
            acmg: "",
            description: "",
          };
          return (
            <div key={tag} className="flex items-start gap-2 bg-gray-50 rounded px-3 py-2">
              <span className="text-xs font-mono bg-blue-100 text-blue-700 px-1.5 py-0.5 rounded flex-shrink-0">
                {info.acmg}
              </span>
              <div>
                <div className="text-sm font-medium text-gray-800">{info.label}</div>
                {info.description && (
                  <div className="text-xs text-gray-500 mt-0.5">{info.description}</div>
                )}
              </div>
            </div>
          );
        })}
      </div>
      <p className="text-xs text-gray-400 mt-3 italic">
        These are computational suggestions, not clinical adjudications.
      </p>
    </div>
  );
}
