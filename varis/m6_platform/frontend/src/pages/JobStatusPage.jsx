import { useState, useEffect, useRef } from "react";
import { useParams, useNavigate } from "react-router-dom";
import { getJobStatus } from "../api/client";
import ErrorBanner from "../components/ErrorBanner";

const PIPELINE_STEPS = [
  "M1: Data Ingestion",
  "M2: Structure Engine",
  "M3: Structural Analysis",
  "M4: Conservation",
  "M5: ML Scoring",
];

function getStepStatus(currentStep, stepName, jobStatus) {
  if (jobStatus === "succeeded") return "complete";
  if (jobStatus === "failed") {
    const currentIdx = PIPELINE_STEPS.indexOf(currentStep);
    const stepIdx = PIPELINE_STEPS.indexOf(stepName);
    if (stepIdx < currentIdx) return "complete";
    if (stepIdx === currentIdx) return "failed";
    return "pending";
  }
  if (!currentStep) return "pending";
  const currentIdx = PIPELINE_STEPS.indexOf(currentStep);
  const stepIdx = PIPELINE_STEPS.indexOf(stepName);
  if (stepIdx < currentIdx) return "complete";
  if (stepIdx === currentIdx) return "running";
  return "pending";
}

export default function JobStatusPage() {
  const { jobId } = useParams();
  const navigate = useNavigate();
  const [job, setJob] = useState(null);
  const [error, setError] = useState(null);
  const intervalRef = useRef(null);

  useEffect(() => {
    function poll() {
      getJobStatus(jobId)
        .then((data) => {
          setJob(data);
          if (data.status === "succeeded" && data.variant_id) {
            clearInterval(intervalRef.current);
            setTimeout(() => navigate(`/variant/${data.variant_id}`), 1000);
          } else if (data.status === "failed") {
            clearInterval(intervalRef.current);
          }
        })
        .catch((err) => {
          setError(err.message);
          clearInterval(intervalRef.current);
        });
    }

    poll();
    intervalRef.current = setInterval(poll, 3000);
    return () => clearInterval(intervalRef.current);
  }, [jobId, navigate]);

  if (error) {
    return (
      <div className="max-w-lg mx-auto px-4 py-16">
        <ErrorBanner error={error} onDismiss={() => setError(null)} />
      </div>
    );
  }

  if (!job) {
    return (
      <div className="max-w-lg mx-auto px-4 py-16 text-center">
        <p className="text-gray-500">Loading job status...</p>
      </div>
    );
  }

  return (
    <div className="max-w-lg mx-auto px-4 py-16">
      <div className="bg-white shadow rounded-lg p-6">
        <h2 className="text-xl font-bold text-gray-900 mb-1">Investigation in Progress</h2>
        <p className="text-sm text-gray-500 mb-6">{job.variant_id}</p>

        <div className="space-y-4">
          {PIPELINE_STEPS.map((step) => {
            const status = getStepStatus(job.current_step, step, job.status);
            return (
              <div key={step} className="flex items-center gap-3">
                <div className="flex-shrink-0 w-6 h-6 flex items-center justify-center">
                  {status === "complete" && (
                    <span className="text-green-500 text-lg">&#10003;</span>
                  )}
                  {status === "running" && (
                    <span className="text-blue-500 text-lg animate-pulse">&#9679;</span>
                  )}
                  {status === "pending" && (
                    <span className="text-gray-300 text-lg">&#9675;</span>
                  )}
                  {status === "failed" && (
                    <span className="text-red-500 text-lg">&#10007;</span>
                  )}
                </div>
                <span className={`text-sm ${
                  status === "running" ? "text-blue-700 font-medium" :
                  status === "complete" ? "text-gray-700" :
                  status === "failed" ? "text-red-600" :
                  "text-gray-400"
                }`}>
                  {step}
                </span>
              </div>
            );
          })}
        </div>

        {job.status === "succeeded" && (
          <div className="mt-6 bg-green-50 border border-green-200 rounded p-3 text-sm text-green-700">
            Investigation complete. Redirecting...
          </div>
        )}

        {job.status === "failed" && job.error_message && (
          <div className="mt-6 bg-red-50 border border-red-200 rounded p-3 text-sm text-red-700">
            {job.error_message}
          </div>
        )}
      </div>
    </div>
  );
}
