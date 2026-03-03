import { useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { submitInvestigation, searchVariants, getStats } from "../api/client";

export default function SearchPage() {
  const navigate = useNavigate();
  const [gene, setGene] = useState("");
  const [hgvs, setHgvs] = useState("");
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(false);
  const [searchQuery, setSearchQuery] = useState("");
  const [searchResults, setSearchResults] = useState([]);
  const [stats, setStats] = useState(null);

  useEffect(() => {
    getStats().then(setStats).catch(() => {});
  }, []);

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

      <div className="flex gap-6">
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
            <div className="bg-white rounded-lg shadow p-4 space-y-2 text-sm">
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
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
