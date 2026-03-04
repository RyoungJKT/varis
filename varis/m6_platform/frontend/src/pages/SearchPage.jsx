import { useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { submitInvestigation, searchVariants, listVariants, getStats } from "../api/client";

export default function SearchPage() {
  const navigate = useNavigate();
  const [gene, setGene] = useState("");
  const [hgvs, setHgvs] = useState("");
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(false);
  const [searchQuery, setSearchQuery] = useState("");
  const [searchResults, setSearchResults] = useState([]);
  const [stats, setStats] = useState(null);
  const [allVariants, setAllVariants] = useState([]);
  const [totalVariants, setTotalVariants] = useState(0);
  const [page, setPage] = useState(1);
  const [loadingList, setLoadingList] = useState(true);
  const ITEMS_PER_PAGE = 20;

  useEffect(() => {
    getStats().then(setStats).catch(() => {});
  }, []);

  useEffect(() => {
    setLoadingList(true);
    listVariants(page, ITEMS_PER_PAGE)
      .then((data) => {
        setAllVariants(data.results || []);
        setTotalVariants(data.total || 0);
      })
      .catch(() => {
        setAllVariants([]);
        setTotalVariants(0);
      })
      .finally(() => setLoadingList(false));
  }, [page]);

  useEffect(() => {
    if (!searchQuery.trim()) {
      setSearchResults([]);
      return;
    }
    const timer = setTimeout(() => {
      searchVariants(searchQuery)
        .then((data) => setSearchResults(data.results || []))
        .catch(() => setSearchResults([]));
    }, 300);
    return () => clearTimeout(timer);
  }, [searchQuery]);

  async function handleSubmit(e) {
    e.preventDefault();
    setError(null);
    setLoading(true);
    try {
      const data = await submitInvestigation(gene.trim(), hgvs.trim());
      if (data.status === "exists") {
        navigate(`/variant/${data.variant_id}`);
      } else if (data.job_id) {
        navigate(`/jobs/${data.job_id}`);
      }
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  }

  return (
    <div className="max-w-4xl mx-auto px-4 py-12">
      <div className="text-center mb-10">
        <h2 className="text-3xl font-bold text-gray-900 mb-2">
          Investigate a Variant
        </h2>
        <p className="text-gray-500">
          Enter a gene and HGVS protein notation to run a structural investigation.
        </p>
      </div>

      <form onSubmit={handleSubmit} className="bg-white shadow rounded-lg p-6 mb-8">
        <div className="flex gap-4 mb-4">
          <div className="flex-1">
            <label className="block text-sm font-medium text-gray-700 mb-1">Gene</label>
            <input
              type="text"
              value={gene}
              onChange={(e) => setGene(e.target.value)}
              placeholder="BRCA1"
              className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm focus:ring-blue-500 focus:border-blue-500"
              required
            />
          </div>
          <div className="flex-2">
            <label className="block text-sm font-medium text-gray-700 mb-1">HGVS Protein</label>
            <input
              type="text"
              value={hgvs}
              onChange={(e) => setHgvs(e.target.value)}
              placeholder="p.Arg1699Trp"
              className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm focus:ring-blue-500 focus:border-blue-500"
              required
            />
          </div>
        </div>
        {error && (
          <div className="text-red-600 text-sm mb-3 bg-red-50 border border-red-200 rounded px-3 py-2">
            {error}
          </div>
        )}
        <button
          type="submit"
          disabled={loading || !gene.trim() || !hgvs.trim()}
          className="w-full bg-blue-600 text-white py-2 px-4 rounded-md hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed text-sm font-medium"
        >
          {loading ? "Submitting..." : "Investigate"}
        </button>
      </form>

      <div className="flex gap-6 mb-8">
        <div className="flex-1">
          <h3 className="text-lg font-semibold text-gray-900 mb-3">Search Existing</h3>
          <input
            type="text"
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            placeholder="Search by gene, variant, or ClinVar ID..."
            className="w-full border border-gray-300 rounded-md px-3 py-2 text-sm mb-3"
          />
          {searchResults.length > 0 ? (
            <ul className="divide-y divide-gray-200 bg-white rounded-lg shadow">
              {searchResults.map((v) => (
                <li key={v.variant_id}>
                  <button
                    onClick={() => navigate(`/variant/${v.variant_id}`)}
                    className="w-full text-left px-4 py-3 hover:bg-gray-50 flex justify-between items-center"
                  >
                    <div>
                      <span className="font-medium text-gray-900">{v.gene}</span>{" "}
                      <span className="text-gray-600">{v.hgvs_protein}</span>
                    </div>
                    {v.classification && (
                      <span className={`text-xs px-2 py-1 rounded-full ${
                        v.classification === "likely_pathogenic"
                          ? "bg-red-100 text-red-700"
                          : v.classification === "likely_benign"
                          ? "bg-green-100 text-green-700"
                          : "bg-yellow-100 text-yellow-700"
                      }`}>
                        {v.classification.replace("_", " ")}
                      </span>
                    )}
                  </button>
                </li>
              ))}
            </ul>
          ) : searchQuery.trim() ? (
            <p className="text-gray-400 text-sm">No variants found.</p>
          ) : null}
        </div>

        {stats && (
          <div className="w-64">
            <h3 className="text-lg font-semibold text-gray-900 mb-3">Database</h3>
            <div
              onClick={() => navigate("/variants")}
              className="bg-white rounded-lg shadow p-4 space-y-2 text-sm cursor-pointer hover:ring-2 hover:ring-blue-300 transition"
            >
              <div className="flex justify-between">
                <span className="text-gray-500">Variants</span>
                <span className="font-medium">{stats.total_variants}</span>
              </div>
              <div className="flex justify-between">
                <span className="text-gray-500">Genes</span>
                <span className="font-medium">{stats.genes_covered}</span>
              </div>
              <div className="flex justify-between">
                <span className="text-gray-500">Pathogenic</span>
                <span className="font-medium text-red-600">{stats.total_pathogenic}</span>
              </div>
              <div className="flex justify-between">
                <span className="text-gray-500">Benign</span>
                <span className="font-medium text-green-600">{stats.total_benign}</span>
              </div>
              <div className="text-center pt-2 border-t border-gray-100">
                <span className="text-xs text-blue-600">View all variants →</span>
              </div>
            </div>
          </div>
        )}
      </div>

      <div>
        <div className="flex items-center justify-between mb-4">
          <h3 className="text-lg font-semibold text-gray-900">
            All Investigated Variants
            {totalVariants > 0 && (
              <span className="text-sm font-normal text-gray-500 ml-2">
                ({totalVariants} total)
              </span>
            )}
          </h3>
        </div>

        {loadingList ? (
          <p className="text-gray-400 text-sm">Loading variants...</p>
        ) : allVariants.length > 0 ? (
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
                  {allVariants.map((v) => (
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

            {totalVariants > ITEMS_PER_PAGE && (
              <div className="flex items-center justify-between mt-4">
                <button
                  onClick={() => setPage((p) => Math.max(1, p - 1))}
                  disabled={page <= 1}
                  className="px-3 py-1.5 text-sm border border-gray-300 rounded-md hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
                >
                  Previous
                </button>
                <span className="text-sm text-gray-500">
                  Page {page} of {Math.ceil(totalVariants / ITEMS_PER_PAGE)}
                </span>
                <button
                  onClick={() => setPage((p) => p + 1)}
                  disabled={page >= Math.ceil(totalVariants / ITEMS_PER_PAGE)}
                  className="px-3 py-1.5 text-sm border border-gray-300 rounded-md hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed"
                >
                  Next
                </button>
              </div>
            )}
          </>
        ) : (
          <div className="bg-white rounded-lg shadow p-8 text-center">
            <p className="text-gray-500">No variants investigated yet.</p>
            <p className="text-gray-400 text-sm mt-1">
              Use the form above to investigate your first variant.
            </p>
          </div>
        )}
      </div>
    </div>
  );
}
