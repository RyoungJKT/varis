import { useState } from "react";

const ERROR_MESSAGES = {
  429: "Too many requests. Please wait a moment and try again.",
  500: "Something went wrong on the server. Please try again later.",
  503: "The service is temporarily unavailable. Please try again later.",
};

export default function ErrorBanner({ error, onDismiss }) {
  const [visible, setVisible] = useState(true);

  if (!error || !visible) return null;

  // Extract status code from error message if present (e.g. "...failed: 429")
  const statusMatch = error.match(/(\d{3})$/);
  const statusCode = statusMatch ? parseInt(statusMatch[1]) : null;
  const message = ERROR_MESSAGES[statusCode] || error;

  const bgColor = statusCode === 429
    ? "bg-amber-50 border-amber-200"
    : "bg-red-50 border-red-200";
  const textColor = statusCode === 429 ? "text-amber-700" : "text-red-700";

  function handleDismiss() {
    setVisible(false);
    onDismiss?.();
  }

  return (
    <div className={`${bgColor} border rounded-lg p-4 flex items-start gap-3`}>
      <div className="flex-1">
        <p className={`text-sm font-medium ${textColor}`}>{message}</p>
      </div>
      <button
        onClick={handleDismiss}
        className="text-gray-400 hover:text-gray-600 text-lg leading-none"
        aria-label="Dismiss error"
      >
        &times;
      </button>
    </div>
  );
}
