export default function PredictionBadge({ prediction }) {
  if (!prediction || prediction.score == null) {
    return (
      <div className="bg-white rounded-lg shadow p-6 text-center">
        <p className="text-gray-400 text-sm">No prediction available yet.</p>
      </div>
    );
  }

  const score = prediction.score;
  const pct = Math.round(score * 100);
  const classification = (prediction.classification || "unknown").replace("_", " ");

  let colorClass = "text-yellow-600";
  let bgClass = "bg-yellow-50 border-yellow-200";
  let ringColor = "#ca8a04";
  if (score >= 0.7) {
    colorClass = "text-red-600";
    bgClass = "bg-red-50 border-red-200";
    ringColor = "#dc2626";
  } else if (score < 0.3) {
    colorClass = "text-green-600";
    bgClass = "bg-green-50 border-green-200";
    ringColor = "#16a34a";
  }

  const circumference = 2 * Math.PI * 36;
  const dashoffset = circumference * (1 - score);

  return (
    <div className={`rounded-lg shadow p-6 border ${bgClass}`}>
      <div className="flex items-center gap-6">
        <div className="relative w-24 h-24 flex-shrink-0">
          <svg viewBox="0 0 80 80" className="w-24 h-24">
            <circle cx="40" cy="40" r="36" fill="none" stroke="#e5e7eb" strokeWidth="6" />
            <circle
              cx="40" cy="40" r="36" fill="none"
              stroke={ringColor} strokeWidth="6"
              strokeLinecap="round"
              strokeDasharray={circumference}
              strokeDashoffset={dashoffset}
              transform="rotate(-90 40 40)"
            />
          </svg>
          <div className="absolute inset-0 flex items-center justify-center">
            <span className={`text-xl font-bold ${colorClass}`}>{pct}%</span>
          </div>
        </div>

        <div>
          <div className={`text-lg font-semibold capitalize ${colorClass}`}>
            {classification}
          </div>
          {prediction.confidence_lower != null && prediction.confidence_upper != null && (
            <div className="text-xs text-gray-500 mt-1">
              CI: {Math.round(prediction.confidence_lower * 100)}% &ndash; {Math.round(prediction.confidence_upper * 100)}%
            </div>
          )}
          {prediction.model_agreement && (
            <div className="mt-2">
              <span className={`text-xs px-2 py-0.5 rounded-full ${
                prediction.model_agreement === "high"
                  ? "bg-green-100 text-green-700"
                  : prediction.model_agreement === "moderate"
                  ? "bg-yellow-100 text-yellow-700"
                  : "bg-red-100 text-red-700"
              }`}>
                {prediction.model_agreement} agreement
              </span>
            </div>
          )}
          {prediction.features_used != null && (
            <div className="text-xs text-gray-400 mt-1">
              {prediction.features_used} features used
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
