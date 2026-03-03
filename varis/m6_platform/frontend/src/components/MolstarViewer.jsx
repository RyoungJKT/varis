import { useRef, useEffect, useState } from "react";

export default function MolstarViewer({ structure }) {
  const containerRef = useRef(null);
  const [loaded, setLoaded] = useState(false);
  const viewerRef = useRef(null);

  const confidence = structure?.coordinate_mapping_confidence;
  const residueIndex = structure?.residue_index;
  const source = structure?.source;
  const plddtBucket = structure?.confidence_bucket;

  // Determine if we can highlight the residue
  const canHighlight = confidence !== "failed" && residueIndex != null;
  const showWarning = confidence === "low";
  const showError = confidence === "failed";

  useEffect(() => {
    if (!containerRef.current || !source) return;

    // Only initialize pdbe-molstar if it's available (loaded from CDN)
    if (typeof window.PDBeMolstarPlugin === "undefined") {
      setLoaded(false);
      return;
    }

    const viewerInstance = new window.PDBeMolstarPlugin();
    viewerRef.current = viewerInstance;

    const options = {
      hideControls: true,
      bgColor: { r: 249, g: 250, b: 251 },
      alphafoldView: source === "alphafold",
    };

    // For AlphaFold structures, use UniProt ID
    if (source === "alphafold" && structure.ref_aa) {
      // Use AlphaFold DB with a dummy structure for now
      // In production, this would load from the actual PDB path
      options.customData = {
        url: structure.pdb_url || "",
        format: "pdb",
      };
    }

    viewerInstance.render(containerRef.current, options);
    setLoaded(true);

    return () => {
      if (viewerRef.current) {
        try {
          viewerRef.current.clear();
        } catch {
          // Viewer already disposed
        }
      }
    };
  }, [source, structure]);

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
