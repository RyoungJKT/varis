import { useState, useEffect } from "react";
import { useNavigate, useSearchParams } from "react-router-dom";
import { listVariants, getStats } from "../api/client";

export default function VariantListPage() {
  const navigate = useNavigate();
  const [searchParams, setSearchParams] = useSearchParams();
  const filterParam = searchParams.get("filter");

  const [variants, setVariants] = useState([]);
  const [total, setTotal] = useState(0);
  const [page, setPage] = useState(1);
  const [stats, setStats] = useState(null);
  const [loading, setLoading] = useState(true);
  const ITEMS_PER_PAGE = 20;

  useEffect(() => {
    getStats().then(setStats).catch(() => {});
  }, []);

  useEffect(() => {
    setLoading(true);
    listVariants(page, ITEMS_PER_PAGE)
      .then((data) => {
        let results = data.results || [];
        if (filterParam === "pathogenic") {
          results = results.filter((v) => v.classification === "likely_pathogenic");
        } else if (filterParam === "benign") {
          results = results.filter((v) => v.classification === "likely_benign");
        }
        setVariants(results);
        setTotal(data.total || 0);
      })
      .catch(() => {
        setVariants([]);
        setTotal(0);
      })
      .finally(() => setLoading(false));
  }, [page, filterParam]);

  const filterLabel = filterParam === "pathogenic"
    ? "Pathogenic Variants"
    : filterParam === "benign"
    ? "Benign Variants"
    : "All Investigated Variants";

  return (
    <div className="max-w-5xl mx-auto px-4 py-8">
      <button
        onClick={() => navigate("/")}
        className="flex items-center gap-1 text-sm text-gray-500 hover:text-gray-800 mb-6"
      >
        <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" />
        </svg>
        Back to search
      </button>

      <div className="flex items-center justify-between mb-6">
        <h2 className="text-2xl font-bold text-gray-900">
          {filterLabel}
          {total > 0 && !filterParam && (
            <span className="text-base font-normal text-gray-500 ml-2">({total})</span>
          )}
        </h2>

        {filterParam && (
          <button
            onClick={() => setSearchParams({})}
            className="text-sm text-blue-600 hover:text-blue-800"
          >
            Show all
          </button>
        )}
      </div>

      {stats && (
        <div className="grid grid-cols-4 gap-4 mb-8">
          <button
            onClick={() => setSearchParams({})}
            className={`bg-white rounded-lg shadow p-4 text-center hover:ring-2 hover:ring-blue-300 transition ${
              !filterParam ? "ring-2 ring-blue-500" : ""
            }`}
          >
            <div className="text-2xl font-bold text-gray-900">{stats.total_variants}</div>
            <div className="text-sm text-gray-500">Variants</div>
          </button>
          <div className="bg-white rounded-lg shadow p-4 text-center">
            <div className="text-2xl font-bold text-gray-900">{stats.genes_covered}</div>
            <div className="text-sm text-gray-500">Genes</div>
          </div>
          <button
            onClick={() => setSearchParams({ filter: "pathogenic" })}
            className={`bg-white rounded-lg shadow p-4 text-center hover:ring-2 hover:ring-red-300 transition ${
              filterParam === "pathogenic" ? "ring-2 ring-red-500" : ""
            }`}
          >
            <div className="text-2xl font-bold text-red-600">{stats.total_pathogenic}</div>
            <div className="text-sm text-gray-500">Pathogenic</div>
          </button>
          <button
            onClick={() => setSearchParams({ filter: "benign" })}
            className={`bg-white rounded-lg shadow p-4 text-center hover:ring-2 hover:ring-green-300 transition ${
              filterParam === "benign" ? "ring-2 ring-green-500" : ""
            }`}
          >
            <div className="text-2xl font-bold text-green-600">{stats.total_benign}</div>
            <div className="text-sm text-gray-500">Benign</div>
          </button>
        </div>
      )}

      {loading ? (
        <p className="text-gray-400 text-sm">Loading variants...</p>
      ) : variants.length > 0 ? (
        <>
          <div className="bg-white rounded-lg shadow overflow-hidden">
            <table className="w-full text-sm">
              <thead className="bg-gray-50 border-b border-gray-200">
                <tr>
                  <th className="text-left px-4 py-3 font-medium text-gray-600">Gene</th>
                  <th className="text-left px-4 py-3 font-medium text-gray-600">Variant</th>
                  <th className="text-left px-4 py-3 font-medium text-gray-600">Classification</th>
                  <th className="text-right px-4 py-3 font-medium text-gray-600">Score</th>
                  <th className="text-right px-4 py-3 font-medium text-gray-600">Date</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-100">
                {variants.map((v) => (
                  <tr
                    key={v.variant_id}
                    onClick={() => navigate(`/variant/${v.variant_id}`)}
                    className="hover:bg-gray-50 cursor-pointer"
                  >
                    <td className="px-4 py-3 font-medium text-gray-900">{v.gene}</td>
                    <td className="px-4 py-3 text-gray-600">{v.hgvs_protein}</td>
                    <td className="px-4 py-3">
                      {v.classification ? (
                        <span className={`text-xs px-2 py-1 rounded-full ${
                          v.classification === "likely_pathogenic"
                            ? "bg-red-100 text-red-700"
                            : v.classification === "likely_benign"
                            ? "bg-green-100 text-green-700"
                            : "bg-yellow-100 text-yellow-700"
                        }`}>
                          {v.classification.replace("_", " ")}
                        </span>
                      ) : (
                        <span className="text-gray-400">--</span>
                      )}
                    </td>
                    <td className="px-4 py-3 text-right text-gray-600">
                      {v.score != null ? v.score.toFixed(3) : "--"}
                    </td>
                    <td className="px-4 py-3 text-right text-gray-400">
                      {v.investigation_timestamp
                        ? new Date(v.investigation_timestamp).toLocaleDateString()
                        : "--"}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          {!filterParam && total > ITEMS_PER_PAGE && (
            <div className="flex items-center justify-between mt-4">
              <button
                onClick={() => setPage((p) => Math.max(1, p - 1))}
                disabled={page <= 1}
                className="px-3 py-1.5 text-sm border border-gray-300 rounded-md hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
              >
                Previous
              </button>
              <span className="text-sm text-gray-500">
                Page {page} of {Math.ceil(total / ITEMS_PER_PAGE)}
              </span>
              <button
                onClick={() => setPage((p) => p + 1)}
                disabled={page >= Math.ceil(total / ITEMS_PER_PAGE)}
                className="px-3 py-1.5 text-sm border border-gray-300 rounded-md hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
              >
                Next
              </button>
            </div>
          )}
        </>
      ) : (
        <div className="bg-white rounded-lg shadow p-8 text-center">
          <p className="text-gray-500">
            {filterParam ? `No ${filterParam} variants found.` : "No variants investigated yet."}
          </p>
        </div>
      )}
    </div>
  );
}
