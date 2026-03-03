export default function ProvenanceFooter({ provenance }) {
  if (!provenance) return null;

  return (
    <div className="bg-gray-50 border-t border-gray-200 rounded-b-lg px-6 py-4 text-xs text-gray-500">
      <div className="flex flex-wrap gap-x-6 gap-y-1">
        {provenance.data_sources?.length > 0 && (
          <div>
            <span className="font-medium text-gray-600">Sources:</span>{" "}
            {provenance.data_sources.join(", ")}
          </div>
        )}
        {provenance.pipeline_version && (
          <div>
            <span className="font-medium text-gray-600">Pipeline:</span>{" "}
            {provenance.pipeline_version}
          </div>
        )}
        {provenance.modules_completed?.length > 0 && (
          <div>
            <span className="font-medium text-gray-600">Modules:</span>{" "}
            {provenance.modules_completed.join(", ")}
          </div>
        )}
        {provenance.investigation_timestamp && (
          <div>
            <span className="font-medium text-gray-600">Investigated:</span>{" "}
            {new Date(provenance.investigation_timestamp).toLocaleDateString()}
          </div>
        )}
        {provenance.processing_time_seconds != null && (
          <div>
            <span className="font-medium text-gray-600">Time:</span>{" "}
            {provenance.processing_time_seconds}s
          </div>
        )}
      </div>
    </div>
  );
}
