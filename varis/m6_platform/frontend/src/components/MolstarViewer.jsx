import { useRef, useEffect, useState } from "react";

export default function MolstarViewer({ structure }) {
  const containerRef = useRef(null);
  const viewerRef = useRef(null);
  const [loaded, setLoaded] = useState(false);
  const [error, setError] = useState(null);

  const confidence = structure?.coordinate_mapping_confidence;
  const residueIndex = structure?.residue_index;
  const source = structure?.source;
  const plddtBucket = structure?.confidence_bucket;
  const uniprotId = structure?.uniprot_id;
  const pdbUrl = structure?.pdb_url;

  // Determine if we can highlight the residue
  const canHighlight = confidence !== "failed" && residueIndex != null;
  const showWarning = confidence === "low";
  const showError = confidence === "failed";

  useEffect(() => {
    if (!containerRef.current || !source) return;

    // Check that the plugin script loaded
    if (!window.PDBeMolstarPlugin) {
      console.error("PDBeMolstarPlugin not found — CDN script may not have loaded");
      setError("3D viewer failed to load");
      return;
    }

    const viewer = new window.PDBeMolstarPlugin();
    viewerRef.current = viewer;

    const options = {
      hideControls: true,
      hideCanvasControls: ["selection", "animation", "controlToggle", "controlInfo"],
      bgColor: { r: 249, g: 250, b: 251 },
      subscribeEvents: false,
      landscape: true,
    };

    if (pdbUrl) {
      // Use relative URL — goes through Vite proxy in dev, same-origin in prod
      options.customData = {
        url: pdbUrl,
        format: "pdb",
        binary: false,
      };
      if (source === "alphafold") {
        options.alphafoldView = true;
      }
    } else if (source === "alphafold" && uniprotId) {
      // Fallback: fetch directly from AlphaFold DB
      options.customData = {
        url: `https://alphafold.ebi.ac.uk/files/AF-${uniprotId}-F1-model_v4.cif`,
        format: "cif",
        binary: false,
      };
      options.alphafoldView = true;
    }

    viewer
      .render(containerRef.current, options)
      .then(() => {
        setLoaded(true);
      })
      .catch((err) => {
        console.error("Mol* render failed:", err);
        setError("Structure could not be displayed");
      });

    return () => {
      try {
        viewer.clear();
      } catch {
        // ignore cleanup errors
      }
      viewerRef.current = null;
    };
  }, [source, uniprotId, pdbUrl]);

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

      {/* Fallback / error states */}
      {!loaded && (
        <div className="absolute inset-0 flex flex-col items-center justify-center bg-gray-100 rounded">
          {error ? (
            <p className="text-sm text-red-400">{error}</p>
          ) : (
            <>
              <svg className="w-16 h-16 text-gray-300 mb-3" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5">
                <path d="M12 3L2 9l10 6 10-6-10-6z" />
                <path d="M2 17l10 6 10-6" />
                <path d="M2 13l10 6 10-6" />
              </svg>
              <p className="text-sm text-gray-400">Loading 3D structure...</p>
              {source && (
                <p className="text-xs text-gray-300 mt-1">Source: {source}</p>
              )}
            </>
          )}
        </div>
      )}
    </div>
  );
}
