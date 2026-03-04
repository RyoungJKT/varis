const API_BASE = "/api/v1";

export async function getInvestigation(variantId) {
  const resp = await fetch(`${API_BASE}/investigations/${encodeURIComponent(variantId)}`);
  if (!resp.ok) throw new Error(`Investigation not found: ${resp.status}`);
  return resp.json();
}

export async function submitInvestigation(gene, hgvs) {
  const resp = await fetch(`${API_BASE}/investigations`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ gene, hgvs }),
  });
  if (!resp.ok) {
    const err = await resp.json().catch(() => ({}));
    throw new Error(err.detail || `Submission failed: ${resp.status}`);
  }
  return resp.json();
}

export async function getJobStatus(jobId) {
  const resp = await fetch(`${API_BASE}/jobs/${encodeURIComponent(jobId)}`);
  if (!resp.ok) throw new Error(`Job not found: ${resp.status}`);
  return resp.json();
}

export async function searchVariants(query) {
  const resp = await fetch(`${API_BASE}/variants?q=${encodeURIComponent(query)}`);
  if (!resp.ok) throw new Error(`Search failed: ${resp.status}`);
  return resp.json();
}

export async function listVariants(page = 1, limit = 20) {
  const resp = await fetch(`${API_BASE}/variants?page=${page}&limit=${limit}`);
  if (!resp.ok) throw new Error(`List failed: ${resp.status}`);
  return resp.json();
}

export async function getStats() {
  const resp = await fetch(`${API_BASE}/variants/stats`);
  if (!resp.ok) throw new Error(`Stats failed: ${resp.status}`);
  return resp.json();
}
