export default function ReliabilityStrip({ structure, features }) {
  const plddt = structure?.plddt_at_residue;
  const plddtBucket = structure?.confidence_bucket;

  let plddtColor = "bg-gray-200 text-gray-600";
  let plddtLabel = "N/A";
  if (plddt != null) {
    if (plddt >= 90) {
      plddtColor = "bg-green-100 text-green-700";
      plddtLabel = `pLDDT ${Math.round(plddt)}`;
    } else if (plddt >= 70) {
      plddtColor = "bg-yellow-100 text-yellow-700";
      plddtLabel = `pLDDT ${Math.round(plddt)}`;
    } else {
      plddtColor = "bg-red-100 text-red-700";
      plddtLabel = `pLDDT ${Math.round(plddt)}`;
    }
  }

  const featureIndicators = (features || []).map((f) => ({
    name: f.name.replace(/_/g, " "),
    available: f.available,
    reason: f.null_reason,
  }));

  return (
    <div className="bg-white rounded-lg shadow px-4 py-3">
      <div className="flex items-center gap-3 flex-wrap">
        {/* pLDDT badge */}
        <span className={`text-xs font-medium px-2.5 py-1 rounded-full ${plddtColor}`}>
          {plddtLabel}
        </span>

        {/* Coordinate mapping confidence */}
        {structure?.coordinate_mapping_confidence && (
          <span className={`text-xs px-2.5 py-1 rounded-full ${
            structure.coordinate_mapping_confidence === "failed"
              ? "bg-red-100 text-red-700"
              : structure.coordinate_mapping_confidence === "low"
              ? "bg-yellow-100 text-yellow-700"
              : "bg-green-100 text-green-700"
          }`}>
            mapping: {structure.coordinate_mapping_confidence}
          </span>
        )}

        <div className="h-4 w-px bg-gray-200" />

        {/* Feature availability indicators */}
        {featureIndicators.map((f) => (
          <span
            key={f.name}
            title={f.available ? `${f.name}: available` : `${f.name}: ${f.reason || "unavailable"}`}
            className={`text-xs cursor-default ${f.available ? "text-green-600" : "text-red-400"}`}
          >
            {f.available ? "\u2713" : "\u2717"}{" "}
            <span className="text-gray-500">{f.name}</span>
          </span>
        ))}
      </div>
    </div>
  );
}
