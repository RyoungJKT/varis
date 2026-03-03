import { BrowserRouter, Routes, Route, Link } from "react-router-dom";
import SearchPage from "./pages/SearchPage";
import InvestigationPage from "./pages/InvestigationPage";
import JobStatusPage from "./pages/JobStatusPage";

export default function App() {
  return (
    <BrowserRouter>
      <div className="min-h-screen bg-gray-50">
        <header className="bg-white shadow-sm border-b">
          <div className="max-w-7xl mx-auto px-4 py-3 flex items-center justify-between">
            <Link to="/" className="text-xl font-bold text-gray-900 no-underline">
              Varis<span className="text-blue-600">DB</span>
            </Link>
            <span className="text-sm text-gray-400">Structural Variant Investigation</span>
          </div>
        </header>
        <Routes>
          <Route path="/" element={<SearchPage />} />
          <Route path="/variant/:variantId" element={<InvestigationPage />} />
          <Route path="/jobs/:jobId" element={<JobStatusPage />} />
        </Routes>
      </div>
    </BrowserRouter>
  );
}
