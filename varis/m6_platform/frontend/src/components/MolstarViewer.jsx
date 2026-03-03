import { useRef, useEffect, useState } from "react";

export default function MolstarViewer({ structure }) {
  const containerRef = useRef(null);
  const [loaded, setLoaded] = useState(false);

  const confidence = structure?.coordinate_mapping_confidence;
  const residueIndex = structure?.residue_index;
  const source = structure?.source;
  const plddtBucket = structure?.confidence_bucket;
  const uniprotId = structure?.uniprot_id;

  // Determine if we can highlight the residue
  const canHighlight = confidence !== "failed" && residueIndex != null;
  const showWarning = confidence === "low";
  const showError = confidence === "failed";

  useEffect(() => {
    if (!containerRef.current || !source) return;

    // Use the <pdbe-molstar> web component
    const el = document.createElement("pdbe-molstar");

    // Prefer local PDB file served by our API; fall back to AlphaFold DB
    const pdbUrl = structure?.pdb_url;
    if (pdbUrl) {
      // Resolve relative URL against the API origin
      const apiBase = "http://localhost:8000";
      el.setAttribute("custom-data-url", apiBase + pdbUrl);
      el.setAttribute("custom-data-format", "pdb");
      if (source === "alphafold") {
        el.setAttribute("alphafold-view", "true");
      }
    } else if (source === "alphafold" && uniprotId) {
      el.setAttribute("molecule-id", uniprotId);
      el.setAttribute("alphafold-view", "true");
    }

    el.setAttribute("hide-controls", "true");
    el.setAttribute("hide-canvas-controls", "selection,animation,controlToggle,controlInfo");
    el.setAttribute("bg-color-r", "249");
    el.setAttribute("bg-color-g", "250");
    el.setAttribute("bg-color-b", "251");
    el.style.width = "100%";
    el.style.height = "100%";
    el.style.display = "block";

    containerRef.current.innerHTML = "";
    containerRef.current.appendChild(el);
    setLoaded(true);

    return () => {
      if (containerRef.current) {
        containerRef.current.innerHTML = "";
      }
    };
  }, [source, uniprotId]);

  return (
    <div className="relative h-full">
      {/* Warning banners */}
      {showWarning && (
        <div className="absolute top-2 left-2 right-2 z-10 bg-yellow-50 border border-yellow-300 rounded px-3 py-1.5 text-xs text-yellow-700">
          Residue mapping uncertain &mdash; highlighted position may be approximate.
        </div>
      )}
      {showError && (
        <div className="absolute top-2 left-2 right-2 z-10 bg-red-50 border border-red-300 rounded px-3 py-1.5 text-xs text-red-700">
          Cannot locate mutation residue in structure. No residue highlighted.
        </div>
      )}

      {/* pLDDT low confidence badge */}
      {plddtBucket === "low" || plddtBucket === "very_low" ? (
        <div className="absolute bottom-2 left-2 z-10 bg-red-50 border border-red-200 rounded px-2 py-1 text-xs text-red-600">
          Low structure confidence
        </div>
      ) : null}

      {/* Mol* container */}
      <div ref={containerRef} className="h-full w-full bg-gray-100 rounded" />

      {/* Fallback if Mol* not loaded */}
      {!loaded && (
        <div className="absolute inset-0 flex flex-col items-center justify-center bg-gray-100 rounded">
          <svg className="w-16 h-16 text-gray-300 mb-3" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5">
            <path d="M12 3L2 9l10 6 10-6-10-6z" />
            <path d="M2 17l10 6 10-6" />
            <path d="M2 13l10 6 10-6" />
          </svg>
          <p className="text-sm text-gray-400">3D Structure Viewer</p>
          {source && (
            <p className="text-xs text-gray-300 mt-1">Source: {source}</p>
          )}
          {canHighlight && (
            <p className="text-xs text-gray-300">Residue: {residueIndex}</p>
          )}
        </div>
      )}
    </div>
  );
}
