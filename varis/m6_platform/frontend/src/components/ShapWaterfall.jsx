import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  Tooltip,
  Cell,
  ResponsiveContainer,
  ReferenceLine,
} from "recharts";

function CustomTooltip({ active, payload }) {
  if (!active || !payload?.length) return null;
  const d = payload[0].payload;
  return (
    <div className="bg-white shadow-lg rounded px-3 py-2 text-xs border">
      <div className="font-medium text-gray-800">{d.feature}</div>
      {d.value != null && (
        <div className="text-gray-500">Value: {d.value}</div>
      )}
      <div className={d.shap >= 0 ? "text-red-600" : "text-blue-600"}>
        SHAP: {d.shap > 0 ? "+" : ""}{d.shap.toFixed(3)}
      </div>
    </div>
  );
}

export default function ShapWaterfall({ explanation, score }) {
  if (!explanation || explanation.length === 0) {
    return null;
  }

  const chartData = explanation.map((item) => ({
    feature: item.feature.replace(/_/g, " "),
    shap: item.shap,
    value: item.value,
  }));

  return (
    <div className="bg-white rounded-lg shadow p-4">
      <div className="flex items-baseline justify-between mb-3">
        <h3 className="text-sm font-semibold text-gray-700">
          Feature Contributions (SHAP)
        </h3>
        {score != null && (
          <span className="text-xs text-gray-500">
            Final score: {Math.round(score * 100)}%
          </span>
        )}
      </div>
      <ResponsiveContainer width="100%" height={Math.max(200, explanation.length * 32)}>
        <BarChart
          data={chartData}
          layout="vertical"
          margin={{ top: 5, right: 30, left: 100, bottom: 5 }}
        >
          <XAxis type="number" tick={{ fontSize: 11 }} />
          <YAxis
            type="category"
            dataKey="feature"
            tick={{ fontSize: 11 }}
            width={95}
          />
          <Tooltip content={<CustomTooltip />} />
          <ReferenceLine x={0} stroke="#9ca3af" />
          <Bar
            dataKey="shap"
            animationDuration={800}
            animationBegin={100}
          >
            {chartData.map((entry, i) => (
              <Cell
                key={i}
                fill={entry.shap >= 0 ? "#ef4444" : "#3b82f6"}
              />
            ))}
          </Bar>
        </BarChart>
      </ResponsiveContainer>
    </div>
  );
}
