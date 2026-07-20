// ------------------------------------------------------------------
// NGS 三級分析 web tool
// Reads tertiary-output TSV/JSON from this server (no GitHub access).
// ------------------------------------------------------------------

const API_BASE = "/api";

// ---------- State ---------------------------------------------------

const state = {
  index:        null,      // [{LIS_ID, Name, MRN, Test, Category}]
  options:      null,      // { category_options: [...] }
  data:         null,      // variants JSON payload
  reports:      null,      // { status, edits, panels, category, updated_at }
  currentLIS:   null,
  dirty:        false,
  diagnosticAnalysisVisible: true,
};
const SAMPLE_TEST_TYPES = ["WES", "WGS", "TITAN-WGS"];
const sampleTestFilters = new Set(SAMPLE_TEST_TYPES);

function normalizeSampleTestType(value = "", sampleId = "") {
  if (/^\d{2}T/i.test(String(sampleId || "").trim())) return "TITAN-WGS";
  return String(value || "").trim().toUpperCase();
}

function isWgsTestType(value) {
  return ["WGS", "TITAN-WGS"].includes(String(value || "").trim().toUpperCase());
}

function currentSampleTestType() {
  const meta = state.data?.meta || {};
  return normalizeSampleTestType(
    meta.Test || meta.test_type || "",
    meta.LIS_ID || meta.lis_id || state.currentLIS || "",
  );
}

function applyDiagnosticAnalysisVisibility() {
  const isTitan = currentSampleTestType() === "TITAN-WGS";
  if (!isTitan) state.diagnosticAnalysisVisible = true;
  const visible = !isTitan || state.diagnosticAnalysisVisible;

  document.body.classList.toggle(
    "titan-diagnostic-analysis-hidden",
    isTitan && !visible,
  );

  const row = document.getElementById("diagnostic-analysis-toggle-row");
  row?.classList.toggle("hidden", !isTitan);
  const btn = document.getElementById("btn-toggle-diagnostic-analysis");
  if (btn) {
    btn.setAttribute("aria-expanded", visible ? "true" : "false");
    btn.textContent = visible ? "▾ 隱藏診斷分析" : "▸ 顯示診斷分析";
  }
}

function setupDiagnosticAnalysisToggle() {
  document.getElementById("btn-toggle-diagnostic-analysis")?.addEventListener("click", () => {
    if (currentSampleTestType() !== "TITAN-WGS") return;
    state.diagnosticAnalysisVisible = !state.diagnosticAnalysisVisible;
    applyDiagnosticAnalysisVisibility();
  });
}

function applyTitanSecondaryFindingsDefault() {
  if (currentSampleTestType() !== "TITAN-WGS") return;
  const btn = document.querySelector(".secondary-findings-toggle");
  const body = btn?.nextElementSibling;
  btn?.setAttribute("aria-expanded", "true");
  body?.classList.add("open");
}

// Tracks blocks the user has manually toggled (by host id).
// If a block id is in this set, we respect its wasOpen dataset;
// otherwise we use defaultOpen from the section def.
const toggledBlocks = new Set();

// lucide-react `Copy` icon, vertically mirrored so the foreground square
// sits top-right and the back outline at bottom-left (the orientation the
// claude.ai code-block button uses).
const COPY_ICON_SVG =
  '<svg viewBox="0 0 24 24" width="12" height="12" fill="none" '
  + 'stroke="currentColor" stroke-width="2" stroke-linecap="round" '
  + 'stroke-linejoin="round" aria-hidden="true">'
  + '<rect width="14" height="14" x="8" y="2" rx="2" ry="2"/>'
  + '<path d="M4 8c-1.1 0-2 .9-2 2v10c0 1.1.9 2 2 2h10c1.1 0 2-.9 2-2"/>'
  + '</svg>';

// ---------- Backend fetch ------------------------------------------

// All requests carry the session cookie; same-origin so credentials
// flow automatically, but spelling it out keeps the intent obvious.
async function apiFetch(path, init = {}) {
  const resp = await fetch(`${API_BASE}${path}`, {
    cache: "no-store",
    credentials: "same-origin",
    headers: { "Accept": "application/json", ...(init.headers || {}) },
    ...init,
  });
  if (resp.status === 401) {
    showLoginModal();
    throw new Error("not authenticated");
  }
  if (resp.status === 404) return null;
  if (!resp.ok) throw new Error(`${resp.status} ${resp.statusText} on ${path}`);
  return await resp.json();
}

async function apiPut(path, body) {
  const resp = await fetch(`${API_BASE}${path}`, {
    method: "PUT",
    credentials: "same-origin",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (resp.status === 401) { showLoginModal(); throw new Error("not authenticated"); }
  if (!resp.ok) throw new Error(`${resp.status} ${resp.statusText} on ${path}`);
  return await resp.json();
}

// ---------- Sample loading -----------------------------------------

async function loadIndex() {
  // Backend returns one entry per sample with keys
  // {sample_id, lis_id, name, mrn, test_type, category, ...}.
  // The combobox / matchSamples below expect the legacy upper-case keys
  // (LIS_ID / Name / MRN / Test / Category), so map both shapes here.
  const data = await apiFetch("/samples");
  const list = Array.isArray(data) ? data : [];
  state.index = list.map(r => {
    const lisId = r.lis_id || r.sample_id || r.LIS_ID || "";
    return {
      LIS_ID:    lisId,
      Name:      r.name || r.Name || "",
      MRN:       r.mrn || r.MRN || "",
      Test:      normalizeSampleTestType(r.test_type || r.Test || "", lisId),
      Category:  r.category || r.Category || "",
      Tag:       (r.tags || []).join(",") || r.Tag || "",
      sample_id: r.sample_id || r.lis_id || r.LIS_ID || "",
    };
  });
  const opts = await apiFetch("/options").catch(() => null);
  state.options = opts && typeof opts === "object"
    ? opts
    : { category_options: [], tag_suggestions: [] };
  return list;
}

function matchSamples(q) {
  if (!state.index) return [];
  const s = (q || "").trim().toLowerCase();
  return state.index.filter(r => {
    if (!sampleTestFilters.has((r.Test || "").toUpperCase())) return false;
    if (!s) return true;
    return ["LIS_ID", "Name", "MRN"].some(k => {
      const v = (r[k] || "").toString().toLowerCase();
      return v && v.includes(s);
    });
  }).slice(0, 50);
}

function resolveLIS(query) {
  if (!state.index) return null;
  const s = (query || "").trim();
  if (!s) return null;
  const visible = state.index.filter(r => sampleTestFilters.has((r.Test || "").toUpperCase()));
  // Exact match on LIS_ID first, then MRN, then single-hit Name
  const hitLis = visible.find(r => (r.LIS_ID || "").trim() === s);
  if (hitLis) return hitLis;
  const hitMrn = visible.find(r => (r.MRN || "").trim() === s);
  if (hitMrn) return hitMrn;
  const hitNames = visible.filter(r => (r.Name || "").trim() === s);
  if (hitNames.length === 1) return hitNames[0];
  if (hitNames.length > 1) throw new Error(`找到 ${hitNames.length} 筆同名樣本，請改用 LIS_ID 或病歷號`);
  // Fallback: single partial match across fields
  const partial = matchSamples(s);
  if (partial.length === 1) return partial[0];
  return null;
}

let _sampleLoadingDepth = 0;
function showSampleLoading() {
  _sampleLoadingDepth += 1;
  document.getElementById("sample-loading-modal")?.classList.remove("hidden");
}
function hideSampleLoading() {
  _sampleLoadingDepth = Math.max(0, _sampleLoadingDepth - 1);
  if (!_sampleLoadingDepth) {
    document.getElementById("sample-loading-modal")?.classList.add("hidden");
  }
}

function resetSampleScopedUiState() {
  activeTierTab = null;
  activeCnvSvTab = null;
  activeMitoTab = null;
  activeStrTab = null;
  document.querySelectorAll(".gene-search-input").forEach(inp => {
    inp.value = "";
  });
  const modalInput = document.getElementById("gene-search-modal-input");
  if (modalInput) {
    modalInput.value = "";
    modalInput.style.display = "";
    modalInput.dataset.kind = "snv";
  }
  const modalTitle = document.getElementById("gene-search-title");
  if (modalTitle) modalTitle.textContent = "基因變異搜尋";
  const modalResults = document.getElementById("gene-search-results");
  if (modalResults) modalResults.innerHTML = "";
  const gnomadFilter = document.getElementById("gene-search-filter-gnomad-af");
  if (gnomadFilter) gnomadFilter.checked = true;
  const omimFilter = document.getElementById("gene-search-filter-omim");
  if (omimFilter) omimFilter.checked = true;
  document.getElementById("gene-search-modal")?.classList.remove("hide-omim");
  hideModal("gene-search-modal");
  _geneSearchToken += 1;
}

async function _loadSample(LIS_ID) {
  // The combobox carries the row's `sample_id` (which equals LIS_ID for
  // legacy samples). Look it up so we use the directory name the backend
  // wants on the URL path.
  const previousTestType = currentSampleTestType();
  const sameTitanSample = state.currentLIS === LIS_ID && previousTestType === "TITAN-WGS";
  const row = (state.index || []).find(r => r.LIS_ID === LIS_ID);
  const sid = row?.sample_id || LIS_ID;
  const data = await apiFetch(`/samples/${encodeURIComponent(sid)}`);
  if (!data) throw new Error(`找不到 sample ${sid}`);
  const reports = await apiFetch(`/samples/${encodeURIComponent(sid)}/report`) || {
    status: {}, edits: {}, panels: {}, clinical_description: "",
    genetic_counseling: "", comment: "",
    category: null, sry_confirmed: false, yield: 0, updated_at: null,
  };
  if (!reports.secondary_findings || typeof reports.secondary_findings !== "object") {
    reports.secondary_findings = {};
  }
  if (reports.clinical_description == null) reports.clinical_description = "";
  if (reports.genetic_counseling   == null) reports.genetic_counseling   = data.genetic_counseling || "";
  if (reports.comment == null)               reports.comment = "";
  if (!Array.isArray(reports.tags))          reports.tags = [];
  if (!Array.isArray(reports.manual_variants)) reports.manual_variants = [];
  if (!Array.isArray(reports.cnv_sv_merges)) reports.cnv_sv_merges = [];
  state.data       = data;
  state.snvSearchVariants = {};
  state.reports    = reports;
  state.currentLIS = LIS_ID;
  const loadedTestType = currentSampleTestType();
  state.diagnosticAnalysisVisible = loadedTestType === "TITAN-WGS"
    ? (sameTitanSample && state.diagnosticAnalysisVisible)
    : true;
  state.dirty      = false;
  _saveError       = "";
  _lastSavedAt     = null;
  clearTimeout(_autoSaveTimer);
  resetSampleScopedUiState();
  // Reset manual-toggle tracking between samples so defaultOpen applies fresh.
  toggledBlocks.clear();
  applyTitanSecondaryFindingsDefault();

  // Staged loading: the core payload above carries empty CNV/SV + Mito + STR + PGx
  // side-channels (aux_pending) and secondary SNV panel categories
  // (secondary_pending). Pull them in the background so the
  // SNV/Indel view + report sections appear immediately; each card
  // re-renders itself when its data lands. A monotonic token drops a
  // stale response that arrives after the user switched samples.
  if (data.aux_pending) {
    const token = (state._auxLoadToken = (state._auxLoadToken || 0) + 1);
    state.cnvPending = true;
    state.svPending = true;
    state.mitoPending  = true;
    state.strPending = true;
    state.pgxPending = true;
    state.secondaryPending = true;
    apiFetch(`/samples/${encodeURIComponent(sid)}/secondary-snv`)
      .then(aux => {
        if (token !== state._auxLoadToken || !state.data) return;
        if (aux) {
          state.data.variants = { ...(state.data.variants || {}), ...(aux.variants || {}) };
          state.data.categories = { ...(state.data.categories || {}), ...(aux.categories || {}) };
          state.data.secondary_pending = false;
        }
        state.secondaryPending = false;
        try { renderReportSections(); } catch (_e) {}
        try { renderCandidateSections(); } catch (_e) {}
      })
      .catch(() => {
        if (token !== state._auxLoadToken) return;
        state.secondaryPending = false;
        if (state.data) state.data.secondary_pending = false;
        try { renderReportSections(); } catch (_e) {}
        try { renderCandidateSections(); } catch (_e) {}
      });
    apiFetch(`/samples/${encodeURIComponent(sid)}/cnv`)
      .then(aux => {
        if (token !== state._auxLoadToken || !state.data) return;
        if (aux) Object.assign(state.data, aux);
        state.cnvPending = false;
        try { renderCnvSvTabBar(); } catch (_e) {}
        // Causative / Other sections may carry CNV/SV ids from
        // state.reports.status that were rendered as "missing
        // variant" while the aux payload was still empty.
        try { renderReportSections(); } catch (_e) {}
      })
      .catch(() => {
        if (token !== state._auxLoadToken) return;
        state.cnvPending = false;
        try { renderCnvSvTabBar(); } catch (_e) {}
      });
    apiFetch(`/samples/${encodeURIComponent(sid)}/sv`)
      .then(aux => {
        if (token !== state._auxLoadToken || !state.data) return;
        if (aux) Object.assign(state.data, aux);
        state.svPending = false;
        try { renderCnvSvTabBar(); } catch (_e) {}
        try { renderReportSections(); } catch (_e) {}
      })
      .catch(() => {
        if (token !== state._auxLoadToken) return;
        state.svPending = false;
        try { renderCnvSvTabBar(); } catch (_e) {}
      });
    apiFetch(`/samples/${encodeURIComponent(sid)}/mito`)
      .then(aux => {
        if (token !== state._auxLoadToken || !state.data) return;
        if (aux) Object.assign(state.data, aux);
        state.mitoPending = false;
        try { renderMitoTabBar(); } catch (_e) {}
        // Same as above for Mito — refresh the report sections so
        // mito ids resolve once the payload lands.
        try { renderReportSections(); } catch (_e) {}
      })
      .catch(() => {
        if (token !== state._auxLoadToken) return;
        state.mitoPending = false;
        try { renderMitoTabBar(); } catch (_e) {}
      });
    apiFetch(`/samples/${encodeURIComponent(sid)}/str`)
      .then(aux => {
        if (token !== state._auxLoadToken || !state.data) return;
        if (aux) Object.assign(state.data, aux);
        state.strPending = false;
        try { renderStrTabBar(); } catch (_e) {}
      })
      .catch(() => {
        if (token !== state._auxLoadToken) return;
        state.strPending = false;
        try { renderStrTabBar(); } catch (_e) {}
      });
    apiFetch(`/samples/${encodeURIComponent(sid)}/pgx`)
      .then(aux => {
        if (token !== state._auxLoadToken || !state.data) return;
        if (aux) Object.assign(state.data, aux);
        state.pgxPending = false;
        try { renderPharmcatBlock("sec-pharmcat"); } catch (_e) {}
        try { renderPharmcatBlock("cat-pharmcat-c"); } catch (_e) {}
      })
      .catch(() => {
        if (token !== state._auxLoadToken) return;
        state.pgxPending = false;
        try { renderPharmcatBlock("sec-pharmcat"); } catch (_e) {}
        try { renderPharmcatBlock("cat-pharmcat-c"); } catch (_e) {}
      });
  } else {
    state.cnvPending = false;
    state.svPending = false;
    state.mitoPending  = false;
    state.strPending = false;
    state.pgxPending = false;
  }
}

async function loadSample(LIS_ID, opts = {}) {
  const showLoading = opts.showLoading !== false;
  if (showLoading) showSampleLoading();
  try {
    const result = await _loadSample(LIS_ID);
    const sampleInput = document.getElementById("q-lis");
    if (sampleInput) sampleInput.value = LIS_ID || "";
    updateWelcomeVisibility();
    return result;
  } finally {
    if (showLoading) hideSampleLoading();
  }
}

// ---------- Formatting helpers --------------------------------------

const CLINVAR_ABBREV = {
  "Pathogenic": "P",
  "Likely_pathogenic": "LP",
  "Pathogenic/Likely_pathogenic": "P/LP",
  "Uncertain_significance": "VUS",
  "Benign": "B",
  "Likely_benign": "LB",
  "Benign/Likely_benign": "B/LB",
  "Conflicting_classifications_of_pathogenicity": "Conflict",
};

function _formatClinvarPart(part) {
  const text = String(part || "").trim().replace(/^_/, "");
  if (!text) return "";
  const m = text.match(/^(.+?)\((\d+)\)$/);
  if (!m) return CLINVAR_ABBREV[text] || text;
  return (CLINVAR_ABBREV[m[1]] || m[1]) + "(" + m[2] + ")";
}

function _formatClinvarParts(text) {
  return String(text || "")
    .split(/[,&|]/)
    .map(_formatClinvarPart)
    .filter(Boolean)
    .join("|");
}

function formatClinvar(sig, conf, stars) {
  // Treat pipeline placeholders ('.', 'NA', '') as "no ClinVar data"
  // — otherwise the cell renders as `.(0★)` instead of `—`.
  const sigStr = (sig == null ? "" : String(sig)).trim();
  if (!sigStr || sigStr === "." || sigStr.toUpperCase() === "NA" || sigStr.toUpperCase() === "N/A") {
    return "—";
  }
  const starTxt = (stars != null && stars !== "") ? `(${stars}★)` : "";
  if (sigStr.startsWith("Conflicting") && conf) {
    return (_formatClinvarParts(conf) || (CLINVAR_ABBREV[sigStr] || sigStr)) + starTxt;
  }
  if (/[,&|]/.test(sigStr)) {
    const out = _formatClinvarParts(sigStr);
    if (out) return out + starTxt;
  }
  return (CLINVAR_ABBREV[sigStr] || sigStr) + starTxt;
}

// Map any ClinVar / ACMG classification string (canonical, abbreviated,
// with or without underscores) to one of the five color buckets. Returns
// null when the text doesn't match any known classification — in that
// case no background color is applied (e.g. conflicting calls, blanks,
// free-form user notes). P/LP → P (red), B/LB → B (dark green).
const SIG_CLASSES = ["sig-p", "sig-lp", "sig-vus", "sig-lb", "sig-b"];
function classifySignificance(text) {
  if (text == null) return null;
  const t = String(text).trim().toLowerCase().replace(/_/g, " ");
  switch (t) {
    case "pathogenic":                        case "p":     return "sig-p";
    case "pathogenic/likely pathogenic":
    case "likely pathogenic/pathogenic":      case "p/lp":  return "sig-p";
    case "likely pathogenic":                 case "lp":    return "sig-lp";
    case "uncertain significance":            case "vus":   return "sig-vus";
    case "likely benign":                     case "lb":    return "sig-lb";
    case "benign":                            case "b":     return "sig-b";
    case "benign/likely benign":              case "b/lb":  return "sig-b";
    default:                                                return null;
  }
}

// In-silico annotations use tool-specific evidence calibration.  Do not infer
// a generic PP3/BP4 strength from a score merely because it is 0–1: only the
// tools backed by Pejaver 2022, Bergquist 2025, or ClinGen SVI Splicing are
// labelled PP3/BP4 below.  Other tools are explicitly marked as model cutoffs.
const IN_SILICO_REFERENCES = {
  pejaver:   { url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC9748256/", pmid: "36413997" },
  bergquist: { url: "https://pubmed.ncbi.nlm.nih.gov/40084623/", pmid: "40084623" },
  splicing:  { url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC10357475/", pmid: "37352859" },
  pangolin:  { url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC9022248/", pmid: "35449021" },
  pknn:      { url: "https://doi.org/10.1101/2025.09.24.678417", pmid: "" },
  metarnn:   { url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC9548151/", pmid: "36209109" },
  dann:      { url: "https://pubmed.ncbi.nlm.nih.gov/25338716/", pmid: "25338716" },
  phactboost:{ url: "https://pubmed.ncbi.nlm.nih.gov/38934805/", pmid: "38934805" },
  loftool:   { url: "https://doi.org/10.1093/bioinformatics/btv602", pmid: "27563026" },
  logofunc:  { url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC10688473/", pmid: "38037155" },
  maxentscan:{ url: "https://pubmed.ncbi.nlm.nih.gov/15285897/", pmid: "15285897" },
  pdivas:    { url: "https://pmc.ncbi.nlm.nih.gov/articles/PMC10563346/", pmid: "37817060" },
};

const IN_SILICO_CALIBRATIONS = {
  pknn: {
    lines: [
      "PP3_Strong: ≥ 4", "PP3_Moderate: 2 to < 4",
      "PP3_Supporting: 1 to < 2", "Uncertain: > -1 to < 1",
      "BP4_Supporting: > -2 to ≤ -1", "BP4_Moderate: > -4 to ≤ -2",
      "BP4_Strong: ≤ -4",
    ], reference: IN_SILICO_REFERENCES.pknn,
  },
  alphamissense: {
    lines: [
      "PP3_Strong: ≥ 0.990", "PP3_3-point: 0.972–0.989",
      "PP3_Moderate: 0.906–0.971", "PP3_Supporting: 0.792–0.905",
      "Indeterminate: 0.170–0.791", "BP4_Supporting: 0.100–0.169",
      "BP4_Moderate: 0.071–0.099", "BP4_3-point: ≤ 0.070",
    ], reference: IN_SILICO_REFERENCES.bergquist,
  },
  pangolin: {
    lines: [
      "Predicted splice impact: |score| ≥ 0.20",
      "Indeterminate: |score| < 0.20",
      "這是原始模型 benchmark cutoff，並非通用 PP3/BP4 calibration",
    ], reference: IN_SILICO_REFERENCES.pangolin,
  },
  esm1b: {
    lines: [
      "PP3_Strong: ≤ -24.0", "PP3_3-point: -23.9 to -14.0",
      "PP3_Moderate: -13.9 to -12.2", "PP3_Supporting: -12.1 to -10.7",
      "Indeterminate: -10.6 to -6.4", "BP4_Supporting: -6.3 to -3.2",
      "BP4_Moderate: -3.1 to 8.7", "BP4_3-point: ≥ 8.8",
    ], reference: IN_SILICO_REFERENCES.bergquist,
  },
  varity: {
    lines: [
      "PP3_Strong: ≥ 0.965", "PP3_3-point: 0.915–0.964",
      "PP3_Moderate: 0.842–0.914", "PP3_Supporting: 0.675–0.841",
      "Indeterminate: 0.252–0.674", "BP4_Supporting: 0.117–0.251",
      "BP4_Moderate: 0.064–0.116", "BP4_3-point: 0.037–0.063",
      "BP4_Strong: ≤ 0.036",
    ], reference: IN_SILICO_REFERENCES.bergquist,
  },
  bayesdel: {
    lines: [
      "PP3_Strong: ≥ 0.50", "PP3_Moderate: 0.27 to < 0.50",
      "PP3_Supporting: 0.13 to < 0.27", "Indeterminate: -0.18 to < 0.13",
      "BP4_Supporting: > -0.36 to -0.18", "BP4_Moderate: ≤ -0.36",
    ], reference: IN_SILICO_REFERENCES.pejaver,
  },
  revel: {
    lines: [
      "PP3_Strong: ≥ 0.932", "PP3_Moderate: 0.773–0.931",
      "PP3_Supporting: 0.644–0.772", "Indeterminate: 0.291–0.643",
      "BP4_Supporting: 0.184–0.290", "BP4_Moderate: 0.017–0.183",
      "BP4_Strong: ≤ 0.016",
    ], reference: IN_SILICO_REFERENCES.pejaver,
  },
  spliceai: {
    lines: [
      "PP3: ≥ 0.20", "BP4: ≤ 0.10", "Indeterminate: > 0.10 and < 0.20",
      "須配合 variant context 與 ClinGen SVI splicing decision tree 使用",
    ], reference: IN_SILICO_REFERENCES.splicing,
  },
  metarnn: {
    lines: [
      "Model pathogenic class: ≥ 0.50", "Model benign class: < 0.50",
      "尚無可靠的通用 PP3/BP4 strength calibration；顏色只表示模型分類",
    ], reference: IN_SILICO_REFERENCES.metarnn,
  },
  dann: {
    lines: [
      "分數越高表示模型預測越 deleterious",
      "作者未提供可通用於臨床分類的 cutoff，亦無可靠 PP3/BP4 strength calibration",
    ], reference: IN_SILICO_REFERENCES.dann,
  },
  phactboost: {
    lines: [
      "Model pathogenic class: ≥ 0.50", "Model benign class: < 0.50",
      "尚無可靠的通用 PP3/BP4 strength calibration；顏色只表示模型分類",
    ], reference: IN_SILICO_REFERENCES.phactboost,
  },
  phylop: {
    lines: [
      "PP3_Moderate: ≥ 9.741", "PP3_Supporting: 7.367 to < 9.741",
      "Indeterminate: > 1.879 to < 7.367", "BP4_Supporting: > 0.021 to 1.879",
      "BP4_Moderate: ≤ 0.021", "文獻註明 PhyloP validation 為 marginal",
    ], reference: IN_SILICO_REFERENCES.pejaver,
  },
  gerp: {
    lines: [
      "Indeterminate / no PP3 cutoff: > 2.70",
      "BP4_Supporting: > -4.54 to 2.70", "BP4_Moderate: ≤ -4.54",
    ], reference: IN_SILICO_REFERENCES.pejaver,
  },
  sift: {
    lines: [
      "PP3_Moderate: 0", "PP3_Supporting: > 0 to 0.001",
      "Indeterminate: > 0.001 to < 0.080", "BP4_Supporting: 0.080 to < 0.327",
      "BP4_Moderate: ≥ 0.327",
    ], reference: IN_SILICO_REFERENCES.pejaver,
  },
  loftool: {
    lines: [
      "LoF-intolerant gene quartile: ≤ 0.25", "Other genes: > 0.25",
    ], reference: IN_SILICO_REFERENCES.loftool,
  },
};

// In-silico predictors in display order. The first three tools with
// actual values land on the primary row of the card; the 4th onwards
// go under ▾ More. Empty cells are skipped entirely — `IN_SILICO_TOOLS`
// is just the priority order, not a fixed slot list.
// Ordering is reviewer-defined; remaining populated tools live under More.
const IN_SILICO_TOOLS = [
  { key: "pknn",          label: "P-KNN LLR",     scoreField: "PKNN_LLR",            extraField: "PKNN_evidence" },
  { key: "alphamissense", label: "AlphaMissense", scoreField: "AlphaMissense_score", predField: "AlphaMissense_pred" },
  // Pangolin is signed: -0.87 = strong splice loss, +0.87 = strong
  // splice gain. Classify by |score| so both colour-code the same way;
  // we still display the signed number so reviewers see the direction.
  { key: "pangolin",      label: "Pangolin",      scoreField: "Pangolin_score" },
  { key: "revel",         label: "REVEL",         scoreField: "REVEL_score" },
  { key: "spliceai",      label: "SpliceAI",      scoreField: "SpliceAI_score" },
  // -- everything below this line lives in More when all three above exist --
  { key: "esm1b",         label: "ESM1b",         scoreField: "ESM1b_score",         predField: "ESM1b_pred" },
  { key: "varity",        label: "VARITY_R",      scoreField: "VARITY_R" },
  { key: "bayesdel",      label: "BayesDel",      scoreField: "BayesDel",            predField: "BayesDel_pred" },
  { key: "metarnn",       label: "MetaRNN",       scoreField: "MetaRNN_score" },
  { key: "dann",          label: "DANN",          scoreField: "DANN" },
  { key: "phactboost",    label: "PhactBoost",    scoreField: "PhactBoost" },
  { key: "phylop",        label: "PhyloP",        scoreField: "PhyloP" },
  { key: "gerp",          label: "GERP",          scoreField: "GERP" },
  { key: "sift",          label: "SIFT",          scoreField: "SIFT_score",          predField: "SIFT_pred" },
  { key: "loftool",       label: "LOFTOOL",       scoreField: "LOFTOOL" },
];
const IN_SILICO_PRIMARY_COUNT = 3;

function _hasNum(x) {
  return x != null && x !== "" && Number.isFinite(Number(x));
}

function _evidence(cls, label) { return { cls, label }; }

function _pknnEvidence(text) {
  const value = String(text || "").trim().replace(/[ -]+/g, "_");
  if (/^PP3_(Strong|Very_Strong)$/i.test(value)) return _evidence("sig-p", value);
  if (/^PP3_/i.test(value)) return _evidence("sig-lp", value);
  if (/^BP4_(Strong|Very_Strong|Moderate)$/i.test(value)) return _evidence("sig-b", value);
  if (/^BP4_/i.test(value)) return _evidence("sig-lb", value);
  return _evidence("sig-vus", value || "Uncertain");
}

function _toolEvidence(v, tool) {
  const x = Number(v[tool.scoreField]);
  if (!Number.isFinite(x)) return _evidence("", "Not available");
  switch (tool.key) {
    case "pknn": return _pknnEvidence(v.PKNN_evidence);
    case "alphamissense":
      if (x >= .990) return _evidence("sig-p", "PP3_Strong");
      if (x >= .972) return _evidence("sig-p", "PP3_3-point");
      if (x >= .906) return _evidence("sig-lp", "PP3_Moderate");
      if (x >= .792) return _evidence("sig-lp", "PP3_Supporting");
      if (x >= .170) return _evidence("sig-vus", "Indeterminate");
      if (x >= .100) return _evidence("sig-lb", "BP4_Supporting");
      if (x >= .071) return _evidence("sig-b", "BP4_Moderate");
      return _evidence("sig-b", "BP4_3-point");
    case "pangolin": return Math.abs(x) >= .20
      ? _evidence("sig-lp", "Predicted splice impact")
      : _evidence("sig-vus", "Indeterminate");
    case "esm1b":
      if (x <= -24.0) return _evidence("sig-p", "PP3_Strong");
      if (x <= -14.0) return _evidence("sig-p", "PP3_3-point");
      if (x <= -12.2) return _evidence("sig-lp", "PP3_Moderate");
      if (x <= -10.7) return _evidence("sig-lp", "PP3_Supporting");
      if (x <= -6.4) return _evidence("sig-vus", "Indeterminate");
      if (x <= -3.2) return _evidence("sig-lb", "BP4_Supporting");
      if (x <= 8.7) return _evidence("sig-b", "BP4_Moderate");
      return _evidence("sig-b", "BP4_3-point");
    case "varity":
      if (x >= .965) return _evidence("sig-p", "PP3_Strong");
      if (x >= .915) return _evidence("sig-p", "PP3_3-point");
      if (x >= .842) return _evidence("sig-lp", "PP3_Moderate");
      if (x >= .675) return _evidence("sig-lp", "PP3_Supporting");
      if (x >= .252) return _evidence("sig-vus", "Indeterminate");
      if (x >= .117) return _evidence("sig-lb", "BP4_Supporting");
      if (x >= .064) return _evidence("sig-b", "BP4_Moderate");
      if (x >= .037) return _evidence("sig-b", "BP4_3-point");
      return _evidence("sig-b", "BP4_Strong");
    case "bayesdel":
      if (x >= .50) return _evidence("sig-p", "PP3_Strong");
      if (x >= .27) return _evidence("sig-lp", "PP3_Moderate");
      if (x >= .13) return _evidence("sig-lp", "PP3_Supporting");
      if (x > -.18) return _evidence("sig-vus", "Indeterminate");
      if (x > -.36) return _evidence("sig-lb", "BP4_Supporting");
      return _evidence("sig-b", "BP4_Moderate");
    case "revel":
      if (x >= .932) return _evidence("sig-p", "PP3_Strong");
      if (x >= .773) return _evidence("sig-lp", "PP3_Moderate");
      if (x >= .644) return _evidence("sig-lp", "PP3_Supporting");
      if (x >= .291) return _evidence("sig-vus", "Indeterminate");
      if (x >= .184) return _evidence("sig-lb", "BP4_Supporting");
      if (x > .016) return _evidence("sig-b", "BP4_Moderate");
      return _evidence("sig-b", "BP4_Strong");
    case "spliceai":
      if (x >= .20) return _evidence("sig-lp", "PP3");
      if (x <= .10) return _evidence("sig-b", "BP4");
      return _evidence("sig-vus", "Indeterminate");
    case "metarnn": return x >= .50
      ? _evidence("sig-lp", "Model pathogenic") : _evidence("sig-lb", "Model benign");
    case "dann": return _evidence("", "No calibrated evidence");
    case "phactboost": return x >= .50
      ? _evidence("sig-lp", "Model pathogenic") : _evidence("sig-lb", "Model benign");
    case "phylop":
      if (x >= 9.741) return _evidence("sig-lp", "PP3_Moderate");
      if (x >= 7.367) return _evidence("sig-lp", "PP3_Supporting");
      if (x > 1.879) return _evidence("sig-vus", "Indeterminate");
      if (x > .021) return _evidence("sig-lb", "BP4_Supporting");
      return _evidence("sig-b", "BP4_Moderate");
    case "gerp":
      if (x > 2.70) return _evidence("sig-vus", "No calibrated PP3");
      if (x > -4.54) return _evidence("sig-lb", "BP4_Supporting");
      return _evidence("sig-b", "BP4_Moderate");
    case "sift":
      if (x === 0) return _evidence("sig-lp", "PP3_Moderate");
      if (x <= .001) return _evidence("sig-lp", "PP3_Supporting");
      if (x < .080) return _evidence("sig-vus", "Indeterminate");
      if (x < .327) return _evidence("sig-lb", "BP4_Supporting");
      return _evidence("sig-b", "BP4_Moderate");
    case "loftool": return x <= .25
      ? _evidence("sig-lp", "LoF-intolerant gene quartile")
      : _evidence("sig-vus", "Gene-level context only");
    default: return _evidence("sig-vus", "No calibrated evidence");
  }
}

function _annotationHint(title, lines, reference = null, options = {}) {
  const body = [`<strong>${escapeHtml(title)}</strong>`]
    .concat((lines || []).map(line => `<div>${escapeHtml(line)}</div>`));
  if (reference?.url) {
    const pmid = reference.pmid ? ` (PMID: ${escapeHtml(reference.pmid)})` : "";
    body.push(`<div class="tip-reference">Reference: <a href="${escapeAttr(reference.url)}" target="_blank" rel="noopener">${escapeHtml(reference.url)}</a>${pmid}</div>`);
  }
  const className = ["info-hint", options.className || ""].filter(Boolean).join(" ");
  const dataId = options.dataId != null ? ` data-id="${escapeAttr(options.dataId)}"` : "";
  return `<button type="button" class="${className}"${dataId} aria-label="${escapeAttr(title)} 註解" data-tip-html="${escapeAttr(body.join(""))}">ⓘ</button>`;
}

function _renderInSilicoCell(v, tool) {
  const score = v[tool.scoreField];
  const has = _hasNum(score);
  const pred = tool.predField ? (v[tool.predField] || "").trim() : "";
  const evidence = has ? _toolEvidence(v, tool) : _evidence("", "Not available");
  const calibration = IN_SILICO_CALIBRATIONS[tool.key] || {};
  const dynamicLines = [`本位點：${evidence.label}`].concat(calibration.lines || []);
  const hint = _annotationHint(tool.label, dynamicLines, calibration.reference);
  const valueTxt = has
    ? (pred ? `${fmtNum(score)} (${escapeHtml(pred)})` : fmtNum(score))
    : "—";
  return `<span class="k">${escapeHtml(tool.label)}${hint}</span>`
       + `<span class="v ${evidence.cls}">${valueTxt}</span>`;
}

// LoGoFunc emits strings like "GOF (0.123)*", "LOF (0.456)", or "Neutral (...)".
// A trailing star means probability > class-specific cutoff (deeper red);
// no star but class is GOF / LOF gives a lighter red. Neutral is a model-level
// negative call (light green), not calibrated ACMG BP4 evidence.
function classifyLoGoFunc(text) {
  if (text == null || text === "" || text === "—") return null;
  const s = String(text).trim();
  const m = s.match(/^(GOF|LOF|Neutral)/i);
  if (!m) return null;
  if (m[1].toUpperCase() === "NEUTRAL") return "sig-lb";
  return s.endsWith("*") ? "sig-p" : "sig-lp";
}

// MaxEntScan_diff has no reliable, general-purpose PP3/BP4 calibration.
// Keep it visibly yellow as contextual evidence instead of mapping an
// arbitrary raw-score difference onto ACMG evidence strength.
function classifyMaxEntScan(score) {
  if (score == null || score === "") return null;
  return Number.isFinite(Number(score)) ? "sig-vus" : null;
}

// PDIVAS paper sensitivity operating points (not PP3/BP4 calibration):
// ~0.501 gives 80% sensitivity; 0.082 gives 95% sensitivity.
// Source: Kurosawa R et al, BMC Genomics 2023.
function classifyPDIVAS(score) {
  if (score == null || score === "") return null;
  const x = Number(score);
  if (!Number.isFinite(x)) return null;
  if (x >= 0.501) return "sig-p";
  if (x >= 0.082) return "sig-lp";
  return "sig-vus";
}

// in_silico_prediction comes through as "<n_pathogenic> - <n_vus> - <n_benign>".
// Render each count as a coloured chip so the direction reads at a glance.
function fmtInSilico(text) {
  if (text == null || text === "") return "—";
  const m = String(text).match(/^\s*(\d+)\s*-\s*(\d+)\s*-\s*(\d+)\s*$/);
  if (!m) return escapeHtml(String(text));
  return `<span class="sig-p">${m[1]}</span> - `
       + `<span class="sig-vus">${m[2]}</span> - `
       + `<span class="sig-b">${m[3]}</span>`;
}

function formatClinvarDate(d) {  if (!d) return "";
  const s = String(d);
  const m = s.match(/^(\d{4})(\d{2})(\d{2})$/);
  if (m) return `${m[1]}-${m[2]}-${m[3]}`;
  return s;
}

function fmtNum(v, digits = 3) {
  if (v == null || v === "") return "—";
  const n = Number(v);
  if (!Number.isFinite(n)) return String(v);
  return n.toFixed(digits).replace(/\.?0+$/, "");
}

function fmtInt(v) {
  if (v == null || v === "") return "—";
  return Math.round(Number(v)).toString();
}

function fmtTxt(v) {
  if (v == null || v === "") return "—";
  return String(v);
}

// "21,18 (0.46)" — AD with VAF in parens. Either half falls back to a dash
// if the underlying field is missing, so a partially populated sample still
// renders cleanly.
// Sum AD ('4,6' / '4,6,0') → total DP. Returns null when AD is blank
// or contains non-numeric parts.
function adTotalDp(ad) {
  if (ad == null || ad === "") return null;
  const parts = String(ad).split(",").map(s => s.trim()).filter(Boolean);
  let sum = 0;
  for (const p of parts) {
    if (p === "." || p.toUpperCase() === "NA") continue;
    const n = Number(p);
    if (!Number.isFinite(n)) return null;
    sum += n;
  }
  return sum;
}

// `sig-lp` (the red ACMG-LP background) when total DP < 10 — clinical
// WGS convention; reviewer should IGV-confirm or Sanger-validate before
// reporting. Empty / unparseable AD → no class.
function lowDpClass(ad, cutoff = 10) {
  const dp = adTotalDp(ad);
  return (dp != null && dp < cutoff) ? "sig-lp" : "";
}

function fmtAdVaf(ad, vaf) {
  const adPart  = (ad == null || ad === "") ? "—" : String(ad);
  const vafPart = (vaf == null || vaf === "" || !Number.isFinite(Number(vaf)))
    ? "—"
    : fmtNum(vaf);
  return `${adPart} (${vafPart})`;
}

// ── IGV.js integration ─────────────────────────────────────────────
// Reviewer-side BAM viewer. Click the IGV button on any SNV / CNV /
// SV card → modal opens with igv.js pinned to that locus, BAM auto-
// detected from the sample's SID (sibling BAMs available as add-ons).
// Backend serves BAM/BAI via /api/igv/file with HTTP range support.

const IGV_SCRIPT_URL = "https://cdn.jsdelivr.net/npm/igv@2/dist/igv.min.js";
const IGV_ALIGNMENT_VISIBILITY_WINDOW = 5000000;
let _igvLoaded = null;       // Promise → resolves when window.igv exists
let _igvBrowser = null;       // current igv.js Browser instance
const _igvBams = [];          // [{label, path, sample_id, batch}]
let _igvSiblings = [];        // candidate add-ons from same batch
let _igvFolderBams = [];      // manual server-folder BAM candidates
let _igvFolderDirs = [];      // known server-side BAM folders
let _igvLocus = "";
let _igvSampleId = "";
let _igvVariant = null;
let _igvIsSry = false;

function _loadIgvScript() {
  if (_igvLoaded) return _igvLoaded;
  _igvLoaded = new Promise((res, rej) => {
    const s = document.createElement("script");
    s.src = IGV_SCRIPT_URL; s.async = true;
    s.onload = () => res(window.igv);
    s.onerror = () => rej(new Error("failed to load igv.js"));
    document.head.appendChild(s);
  });
  return _igvLoaded;
}

function _bamUrl(path) {
  return `${API_BASE}/igv/file?path=${encodeURIComponent(path)}`;
}

function _bamIndexUrl(bam) {
  const indexPath = bam?.index_path || (bam?.path ? `${bam.path}.bai` : "");
  return _bamUrl(indexPath);
}

// Pick a locus string ("chr1:12345-12545") with padding from a variant
// payload. SNV/Indel gets ~100bp context; CNV/SV gets 20% flanks so the
// reviewer can compare coverage around both breakpoints.
function _igvLocusFor(v) {
  if (!v) return "";
  const chrom = _normalizeChrom(v.CHROM || "");
  if (!chrom) return "";
  const start = Number(v.POS);
  const isCnvSv = v.END != null && v.REF == null && v.ALT == null;
  const end = isCnvSv ? Number(v.END) : start + Math.max(1, (v.REF || "A").length);
  if (!Number.isFinite(start)) return chrom;
  if (!Number.isFinite(end)) return chrom;
  let pad = isCnvSv ? Math.max(1, Math.ceil((end - start) * 0.2)) : 100;
  // Keep <=5 Mb CNV/SV events inside the alignment track's default
  // load window. Use as much flanking context as fits without forcing
  // the reviewer to zoom in before BAM coverage appears.
  if (isCnvSv && end - start <= IGV_ALIGNMENT_VISIBILITY_WINDOW) {
    pad = Math.min(pad, Math.max(0, Math.floor((IGV_ALIGNMENT_VISIBILITY_WINDOW - (end - start)) / 2)));
  }
  return `${chrom}:${Math.max(1, start - pad)}-${end + pad}`;
}

function _igvRoiFor(v) {
  if (!v || v.END == null || v.REF != null || v.ALT != null) return undefined;
  const chr = _normalizeChrom(v.CHROM || "");
  const start = Number(v.POS);
  const end = Number(v.END);
  if (!chr || !Number.isFinite(start) || !Number.isFinite(end)) return undefined;
  const source = (v.source || "SV").toUpperCase();
  const svType = v.sv_type || "";
  return [{
    name: `${source} ${svType} ${chr}:${start}-${end}`.trim(),
    color: "rgba(220, 38, 38, 0.18)",
    features: [{ chr, start: Math.max(0, start - 1), end }],
  }];
}

function _igvVariantTitleFor(v) {
  if (!v) return { label: "", coordinate: "" };
  const build = state.data?.genome_build || "hg38";
  if ((v.source || "").toLowerCase() === "mito" || String(v.CHROM || "").replace(/^chr/i, "").toUpperCase() === "M") {
    const hgvs = v.HGVS_M || `m.${v.POS || "?"}${v.REF || ""}>${v.ALT || ""}`;
    const id = v.id || `${v.CHROM}-${v.POS}-${v.REF}-${v.ALT}`;
    return {
      label: hgvs,
      coordinate: `[${build}] ${id}`,
    };
  }
  if (v.igv_label) {
    return {
      label: v.igv_label,
      coordinate: `[${build}] ${_normalizeChrom(v.CHROM || "?")}:${v.POS || "?"}-${v.END || "?"}`,
    };
  }
  if (v.REF != null && v.ALT != null) {
    const id = v.id || `${v.CHROM}-${v.POS}-${v.REF}-${v.ALT}`;
    return {
      label: displaySnvHgvs(v, id),
      coordinate: `[${build}] ${id}`,
    };
  }
  const source = (v.source || "SV").toUpperCase();
  const svType = v.sv_type || "";
  const annotSvId = v.AnnotSV_ID || v.id || "";
  return {
    label: [source, svType, annotSvId].filter(Boolean).join(" "),
    coordinate: `[${build}] ${_normalizeChrom(v.CHROM || "?")}:${v.POS || "?"}-${v.END || "?"}`,
  };
}

function displaySnvHgvs(v, fallback = "") {
  const raw = String(v?.HGVS || fallback || "").trim();
  const transcript = String(v?.transcript || "").trim();
  const parts = raw.split(":");
  const txIndex = parts.findIndex(p => /^ENST\d+(?:\.\d+)?$/i.test(p));
  const enst = txIndex >= 0 ? parts[txIndex] : transcript;
  if (!/^ENST/i.test(enst || "")) return raw || fallback;
  const enstBase = enst.split(".")[0].toUpperCase();
  const mane = Array.isArray(v?.MANE_ALL) ? v.MANE_ALL : [];
  const hit = mane.find(row => {
    const type = String(row?.transcript_type || "").toUpperCase();
    const rowEnst = String(row?.enst || "").split(".")[0].toUpperCase();
    const nm = String(row?.transcript || "");
    return type === "MANE_SELECT" && rowEnst === enstBase && /^NM_/i.test(nm);
  });
  const nm = String(hit?.transcript || "");
  if (!nm) return raw || fallback;
  if (txIndex >= 0) {
    parts[txIndex] = nm;
    return parts.map((part) => {
      if (/^ENSP\d+(?:\.\d+)?$/i.test(part)) return "";
      return part;
    }).filter(Boolean).join(":");
  }
  if (!raw) return nm;
  return raw.replace(enst, nm);
}

function _maneDisplayTranscript(row) {
  const nm = String(row?.transcript || "").trim();
  const enst = String(row?.enst || "").trim();
  if (/^NM_/i.test(nm)) return nm;
  if (/^NM_/i.test(enst)) return enst;
  return nm || enst;
}

function _maneDisplayHgvs(row, key) {
  const value = String(row?.[key] || "").trim();
  const nm = String(row?.refseq_transcript || row?.transcript || "").trim();
  const enst = String(row?.ensembl_transcript || row?.enst || "").trim();
  if (/^NM_/i.test(nm) && enst && value.startsWith(`${enst}:`)) {
    return `${nm}:${value.split(":", 2)[1] || ""}`;
  }
  return value;
}

function _selectedTranscriptOption(v, id) {
  const options = Array.isArray(v?.transcript_options) ? v.transcript_options : [];
  if (!options.length) return null;
  const selected = String(getEdit(id, "selected_transcript_key") || "").trim();
  if (selected) {
    const hit = options.find(opt => String(opt?.key || "") === selected);
    if (hit) return hit;
  }
  const def = String(v?.default_transcript_key || v?.selected_transcript_key || "").trim();
  if (def) {
    const hit = options.find(opt => String(opt?.key || "") === def);
    if (hit) return hit;
  }
  return options[0] || null;
}

function _variantWithSelectedTranscript(v, id) {
  const opt = _selectedTranscriptOption(v, id);
  if (!opt) return v;
  return {
    ...v,
    gene_symbol: opt.gene_symbol || v.gene_symbol,
    transcript: opt.transcript || v.transcript,
    ensembl_transcript: opt.ensembl_transcript || v.ensembl_transcript,
    refseq_transcript: opt.refseq_transcript || v.refseq_transcript,
    refseq_protein: opt.refseq_protein || v.refseq_protein,
    mane_status: opt.mane_status || v.mane_status,
    transcript_type: opt.transcript_type || v.transcript_type,
    HGVS_C: opt.HGVS_C || v.HGVS_C,
    HGVS_P: opt.HGVS_P || v.HGVS_P,
    HGVS: opt.HGVS || v.HGVS,
    Consequence: opt.Consequence || v.Consequence,
    impact: opt.impact || v.impact,
    exon: opt.exon || v.exon,
    intron: opt.intron || v.intron,
    hgnc_id: opt.hgnc_id || v.hgnc_id,
    selected_transcript_key: opt.key || v.selected_transcript_key,
  };
}

function _transcriptOptionLabel(opt) {
  const refseq = String(opt?.refseq_transcript || "").trim();
  const enst = String(opt?.ensembl_transcript || opt?.transcript || "").trim();
  const tx = refseq && enst && refseq !== enst ? `${refseq} / ${enst}` : (refseq || enst);
  const bits = [
    opt?.transcript_type || "",
    opt?.gene_symbol || "",
    tx,
    _maneDisplayHgvs(opt || {}, "HGVS_C") || opt?.HGVS_C || "",
    _maneDisplayHgvs(opt || {}, "HGVS_P") || opt?.HGVS_P || "",
    opt?.Consequence || "",
  ].filter(Boolean);
  return bits.join(" · ");
}

function renderTranscriptPicker(v, id) {
  const options = Array.isArray(v?.transcript_options) ? v.transcript_options : [];
  if (options.length <= 1) return "";
  const selected = _selectedTranscriptOption(v, id)?.key || "";
  const rows = options.map(opt => `
    <option value="${escapeAttr(opt.key || "")}" ${String(opt.key || "") === selected ? "selected" : ""}>
      ${escapeHtml(_transcriptOptionLabel(opt))}
    </option>`).join("");
  return `<span class="transcript-picker" title="選擇這張卡片與 DOCX 報告使用的 transcript">
    <span class="transcript-arrow" aria-hidden="true">▾</span>
    <select class="transcript-select" data-id="${escapeAttr(id)}">${rows}</select>
  </span>`;
}

function _sryIgvRegion() {
  // SRY gene coordinates from the reference assemblies used by the UI.
  // Keep this as a lightweight pseudo-variant so it follows the same
  // modal, BAM lookup, sibling-track and local-reference path as cards.
  return state.data?.genome_build === "hg19"
    ? { source: "gene", sv_type: "SRY", igv_label: "SRY", igv_mode: "sry", CHROM: "chrY", POS: 2654896, END: 2655723 }
    : { source: "gene", sv_type: "SRY", igv_label: "SRY", igv_mode: "sry", CHROM: "chrY", POS: 2786855, END: 2787682 };
}

async function openIgvModal(variant) {
  const modal = document.getElementById("igv-modal");
  if (!modal) return;
  const sid = state.data?.sample_id || state.currentLIS || "";
  _igvSampleId = sid;
  _igvVariant = variant;
  _igvIsSry = variant?.igv_mode === "sry";
  _igvLocus = _igvLocusFor(variant);
  const variantTitle = _igvVariantTitleFor(variant);
  document.getElementById("igv-title").textContent = `IGV — ${sid || "?"}`;
  document.getElementById("igv-locus").textContent = _igvLocus ? `(${_igvLocus})` : "";
  document.getElementById("igv-variant-title").textContent = variantTitle.label;
  document.getElementById("igv-variant-coordinate").textContent =
    variantTitle.coordinate ? `(${variantTitle.coordinate})` : "";
  document.getElementById("igv-bam-hint").textContent = "";
  document.getElementById("igv-load-status").textContent = "";
  document.getElementById("igv-host").innerHTML = "";
  document.getElementById("igv-bam-folder-panel").hidden = true;
  document.getElementById("igv-bam-folder-dir-select").innerHTML = '<option value="">— 選擇 BAM 資料夾 —</option>';
  document.getElementById("igv-bam-folder-select").innerHTML = '<option value="">— 選擇 primary BAM —</option>';
  const loadBtn = document.getElementById("igv-load-btn");
  loadBtn.disabled = true;
  loadBtn.textContent = "載入 IGV";
  _igvBrowser = null;
  _igvBams.length = 0;
  _igvSiblings = [];
  _igvFolderBams = [];
  _igvFolderDirs = [];
  modal.classList.remove("hidden");

  // Look up the primary BAM + sibling list for this sample.
  let bamIndex;
  try {
    bamIndex = await apiFetch(`/igv/bams?sample_id=${encodeURIComponent(sid)}`);
  } catch (e) {
    document.getElementById("igv-bam-hint").textContent = "BAM 查詢失敗：" + (e.message || e);
    _renderIgvBamList();
    return;
  }
  if (bamIndex?.primary) _igvBams.push(bamIndex.primary);
  if (_igvIsSry) {
    _igvSiblings = [...(bamIndex?.siblings || [])];
  } else {
    for (const sib of (bamIndex?.siblings || [])) _igvBams.push(sib);
  }
  const batch = bamIndex?.primary?.batch;
  if (batch) {
    try {
      const more = await apiFetch(`/igv/batch-samples?batch=${encodeURIComponent(batch)}`);
      const resolvedSid = bamIndex?.resolved_sample_id || sid;
      const candidates = [..._igvSiblings, ...(more?.samples || [])]
        .filter(s => s.sample_id !== resolvedSid);
      _igvSiblings = Array.from(new Map(candidates.map(s => [s.path, s])).values());
    } catch {
      if (!_igvIsSry) _igvSiblings = [];
    }
  }
  _renderIgvBamList();

  if (!_igvBams.length) {
    document.getElementById("igv-bam-hint").textContent =
      "找不到對應的 BAM（搜尋路徑：" + (bamIndex?.roots || []).join(", ") + "）";
    return;
  }
  loadBtn.disabled = false;
  document.getElementById("igv-load-status").textContent = "確認 BAM 列表後按「載入 IGV」開始載入";
}

function closeIgvModal() {
  document.getElementById("igv-modal")?.classList.add("hidden");
  if (_igvBrowser && window.igv) {
    try { window.igv.removeBrowser(_igvBrowser); } catch {}
  }
  _igvBrowser = null;
  _igvVariant = null;
  _igvIsSry = false;
  document.getElementById("igv-host").innerHTML = "";
}

async function _initIgvBrowser() {
  const igv = await _loadIgvScript();
  const host = document.getElementById("igv-host");
  host.innerHTML = "";
  const build = state.data?.genome_build === "hg19" ? "hg19" : "hg38";
  const coverageAutoscaleGroup = _igvRoiFor(_igvVariant)
    ? "ngs-ui-cnv-sv-coverage"
    : undefined;
  const variantSource = String(_igvVariant?.source || "").toLowerCase();
  const alignmentTrackHeight = variantSource === "cnv" || variantSource === "sv" ? 50 : 300;
  const tracks = _igvBams.map(b => ({
    name: b.label,
    type: "alignment",
    format: "bam",
    displayMode: "SQUISHED",
    height: alignmentTrackHeight,
    visibilityWindow: IGV_ALIGNMENT_VISIBILITY_WINDOW,
    autoscaleGroup: coverageAutoscaleGroup,
    autoscale: _igvIsSry ? false : undefined,
    max: _igvIsSry ? 100 : undefined,
    url: _bamUrl(b.path),
    indexURL: _bamIndexUrl(b),
  }));
  // igv.js's default hg38 reference is hosted on AWS S3, which is
  // blocked on the hospital intranet. Ask the backend for a custom
  // genome config pointing at our proxied local fasta + fai. hg19 is
  // not in scope yet — falls through to the default string config.
  let genomeOpt = build;
  if (build === "hg38") {
    try {
      const g = await apiFetch(`/igv/genome?build=hg38`);
      if (g?.ok && g.config) genomeOpt = g.config;
      else throw new Error(g?.reason || "no local hg38 reference");
    } catch (e) {
      throw new Error("無法載入本機 hg38 reference：" + (e.message || e)
        + "（在伺服器設 NGS_UI_IGV_REF_DIR 指向 hg38.fa + hg38.fa.fai）");
    }
  }
  _igvBrowser = await igv.createBrowser(host, {
    reference: typeof genomeOpt === "string" ? undefined : genomeOpt,
    genome:    typeof genomeOpt === "string" ? genomeOpt : undefined,
    locus:     _igvLocus || undefined,
    roi:       _igvRoiFor(_igvVariant) || undefined,
    tracks,
  });
}

function _renderIgvBamList() {
  const ul = document.getElementById("igv-bam-list");
  ul.innerHTML = "";
  _igvBams.forEach((b, i) => {
    const li = document.createElement("li");
    li.innerHTML = `
      <span class="igv-bam-label">${escapeHtml(b.label)}</span>
      <span class="igv-bam-path">${escapeHtml(b.path)}</span>
      <button type="button" class="igv-bam-remove" data-idx="${i}" title="移除">×</button>
    `;
    ul.appendChild(li);
  });
  // Refresh the add-dropdown excluding already-loaded BAMs.
  const sel = document.getElementById("igv-bam-add-select");
  const loaded = new Set(_igvBams.map(b => b.path));
  const opts = ["<option value=\"\">— 加入同 batch 的 sample —</option>"];
  for (const s of _igvSiblings) {
    if (loaded.has(s.path)) continue;
    opts.push(`<option value="${escapeAttr(s.path)}">${escapeHtml(s.label)}</option>`);
  }
  sel.innerHTML = opts.join("");
}

async function _loadIgvBamFolders() {
  const select = document.getElementById("igv-bam-folder-dir-select");
  const hint = document.getElementById("igv-bam-hint");
  hint.textContent = "讀取 BAM 路徑中…";
  try {
    const payload = await apiFetch("/igv/bam-folders");
    _igvFolderDirs = payload?.folders || [];
    if (!_igvFolderDirs.length) {
      select.innerHTML = '<option value="">— 沒有可用 BAM 資料夾 —</option>';
      hint.textContent = "找不到可用 BAM 資料夾";
      return;
    }
    select.innerHTML = ['<option value="">— 選擇 BAM 資料夾 —</option>']
      .concat(_igvFolderDirs.map((f) =>
        `<option value="${escapeAttr(f.dir)}">${escapeHtml(f.label)} — ${escapeHtml(f.dir)}</option>`
      )).join("");
    hint.textContent = payload?.truncated
      ? `顯示前 ${_igvFolderDirs.length} 個 BAM 資料夾`
      : `找到 ${_igvFolderDirs.length} 個 BAM 資料夾`;
  } catch (e) {
    _igvFolderDirs = [];
    select.innerHTML = '<option value="">— 選擇 BAM 資料夾 —</option>';
    hint.textContent = "讀取 BAM 路徑失敗：" + (e.message || e);
  }
}

async function _loadIgvBamFolder() {
  const dirSelect = document.getElementById("igv-bam-folder-dir-select");
  const select = document.getElementById("igv-bam-folder-select");
  const hint = document.getElementById("igv-bam-hint");
  const dir = (dirSelect?.value || "").trim();
  if (!dir) {
    hint.textContent = "請先選擇 BAM 資料夾";
    return;
  }
  hint.textContent = "讀取 BAM 清單中…";
  try {
    const payload = await apiFetch(`/igv/bam-folder?dir=${encodeURIComponent(dir)}`);
    _igvFolderBams = payload?.samples || [];
    if (!_igvFolderBams.length) {
      select.innerHTML = '<option value="">— 此資料夾沒有可用 BAM —</option>';
      hint.textContent = "此資料夾沒有可用 BAM（會排除 .repeats.bam）";
      return;
    }
    const sid = String(_igvSampleId || "").replace(/-(dragen|nckuh|inhouse|WES|WGS)$/i, "");
    const preferred = _igvFolderBams.find(b => b.sample_id === sid) || _igvFolderBams[0];
    select.innerHTML = ['<option value="">— 選擇 primary BAM —</option>']
      .concat(_igvFolderBams.map((b) =>
        `<option value="${escapeAttr(b.path)}" ${b.path === preferred.path ? "selected" : ""}>${escapeHtml(b.label)} — ${escapeHtml(b.path)}</option>`
      )).join("");
    hint.textContent = `找到 ${_igvFolderBams.length} 個 BAM`;
  } catch (e) {
    _igvFolderBams = [];
    select.innerHTML = '<option value="">— 選擇 primary BAM —</option>';
    hint.textContent = "讀取資料夾失敗：" + (e.message || e);
  }
}

function _useIgvBamFolderSelection() {
  const select = document.getElementById("igv-bam-folder-select");
  const hint = document.getElementById("igv-bam-hint");
  const path = select?.value || "";
  const primary = _igvFolderBams.find(b => b.path === path);
  if (!primary) {
    hint.textContent = "請先選擇 primary BAM";
    return;
  }
  _igvBams.length = 0;
  _igvBams.push(primary);
  _igvSiblings = _igvFolderBams.filter(b => b.path !== primary.path);
  _renderIgvBamList();
  document.getElementById("igv-load-btn").disabled = false;
  document.getElementById("igv-load-status").textContent = "確認 BAM 列表後按「載入 IGV」開始載入";
  hint.textContent = "已改用其他路徑 BAM；同 batch 清單來自同一個資料夾";
  if (_igvBrowser) _initIgvBrowser().catch(() => {});
}

document.addEventListener("click", (ev) => {
  if (ev.target.closest("#igv-bam-list .igv-bam-remove")) {
    const idx = Number(ev.target.dataset.idx);
    if (!Number.isFinite(idx)) return;
    _igvBams.splice(idx, 1);
    _renderIgvBamList();
    // If IGV is already loaded, reflect the change live; otherwise the
    // edit only affects the next 載入 IGV click.
    if (_igvBrowser) _initIgvBrowser().catch(() => {});
    return;
  }
  if (ev.target.id === "igv-bam-add-btn") {
    const sel = document.getElementById("igv-bam-add-select");
    const path = sel.value;
    if (!path) return;
    const found = _igvSiblings.find(s => s.path === path);
    if (!found || _igvBams.some(b => b.path === path)) return;
    _igvBams.push(found);
    _renderIgvBamList();
    if (_igvBrowser) _initIgvBrowser().catch(() => {});
    return;
  }
  if (ev.target.id === "igv-bam-folder-toggle") {
    const panel = document.getElementById("igv-bam-folder-panel");
    panel.hidden = !panel.hidden;
    if (!panel.hidden) {
      if (!_igvFolderDirs.length) _loadIgvBamFolders();
      document.getElementById("igv-bam-folder-dir-select")?.focus();
    }
    return;
  }
  if (ev.target.id === "igv-bam-folder-refresh-btn") {
    _loadIgvBamFolders();
    return;
  }
  if (ev.target.id === "igv-bam-folder-use-btn") {
    _useIgvBamFolderSelection();
    return;
  }
  if (ev.target.id === "igv-load-btn") {
    const btn = ev.target;
    if (btn.disabled) return;
    btn.disabled = true;
    const status = document.getElementById("igv-load-status");
    status.textContent = "載入中…";
    _initIgvBrowser()
      .then(() => {
        status.textContent = "已載入 ✓";
        btn.textContent = "重新載入";
        btn.disabled = false;
      })
      .catch((e) => {
        document.getElementById("igv-host").textContent =
          "IGV 初始化失敗：" + (e.message || e);
        status.textContent = "";
        btn.disabled = false;
      });
    return;
  }
  if (ev.target.id === "btn-igv-sry") {
    openIgvModal(_sryIgvRegion());
    return;
  }
  // IGV launch buttons on variant / CNV-SV cards.
  const igvBtn = ev.target.closest(".btn-igv");
  if (igvBtn) {
    ev.preventDefault();
    const id = igvBtn.dataset.id;
    const variant = _findVariantById(id);
    openIgvModal(variant);
  }
});

document.addEventListener("keydown", (ev) => {
  if (ev.target.id === "igv-bam-folder-dir-select" && ev.key === "Enter") {
    _loadIgvBamFolder();
  }
});

document.addEventListener("change", (ev) => {
  if (ev.target.id === "igv-bam-folder-dir-select") {
    _loadIgvBamFolder();
  }
});

function _findVariantById(id) {
  // Keep IGV launch lookup aligned with the report renderers. CNV/SV
  // staged loading merges flat `cnv_variants` / `sv_variants` maps
  // into state.data; older nested cnv_sv maps are no longer used.
  return lookupAnyVariant(id).v;
}

function variantUrls(v) {
  const tag = `${v.CHROM}-${v.POS}-${v.REF}-${v.ALT}`;
  // Route Varsome / Franklin / GeneBe to the matching genome build so
  // hg19 samples don't land on an hg38 coordinate page. Build info
  // comes from the R webdata writer (state.data.genome_build); falls
  // back to hg38 for older samples that predate the field.
  const build = state.data?.genome_build === "hg19" ? "hg19" : "hg38";
  const omimId = String(v.OMIM_id || v.omim_ids || "")
    .split(/[;,|\s]+/)
    .map(x => x.trim())
    .find(Boolean);
  const gene = String(v.gene_symbol || "").trim();
  return {
    varsome:  `https://varsome.com/variant/${build}/${tag}`,
    franklin: `https://franklin.genoox.com/clinical-db/variant/snp/${tag}-${build}`,
    genebe:   `https://genebe.net/variant/${build}/${tag}`,
    omim:     v.OMIM_link
      || (omimId ? `https://www.omim.org/entry/${encodeURIComponent(omimId)}` : null)
      || (gene ? `https://www.omim.org/search?index=geneMap&search=${encodeURIComponent(gene)}` : null),
  };
}

// ---------- Render: sample header / phenotype ----------------------

function renderPloidySexStatus(reportedSex) {
  const ploidyCall = String(state.data.ploidy?.karyotype || "").trim().toUpperCase();
  const sex = String(reportedSex || "").trim().toUpperCase();
  const sexControl = document.getElementById("m-sex-control");
  const ploidyLabel = document.getElementById("m-ploidy-call");
  const matches = (
    (sex === "M" && ploidyCall === "XY") ||
    (sex === "F" && ploidyCall === "XX")
  );
  sexControl?.classList.toggle("ploidy-match", !!ploidyCall && matches);
  sexControl?.classList.toggle("ploidy-mismatch", !!ploidyCall && !matches);
  if (ploidyLabel) {
    ploidyLabel.textContent = ploidyCall ? `ploidy VCF: ${ploidyCall}` : "";
    ploidyLabel.hidden = !ploidyCall || matches;
  }
}

function renderSampleMeta() {
  const m = state.data.meta || {};
  document.getElementById("m-lis").textContent       = m.LIS_ID || "—";
  document.getElementById("m-name").textContent      = m.Name || "—";
  document.getElementById("m-mrn").textContent       = m.MRN || "—";
  // "日期" prefers LIS 簽收時間 from the patient list xlsx; falls back
  // to the pipeline-generated timestamp for samples uploaded before the
  // roster carried that column. Strip the time-of-day suffix — the
  // reviewer only cares about the date.
  const dateRaw = m.SignReceivedAt || state.data.generated_at || "";
  document.getElementById("m-generated").textContent = dateRaw ? dateRaw.slice(0, 10) : "—";

  // Copy buttons next to LIS_ID / Name / MRN. Hide when the value is
  // missing so the icon doesn't dangle next to an em-dash.
  const setCopy = (btnId, value) => {
    const btn = document.getElementById(btnId);
    if (!btn) return;
    if (value) {
      btn.dataset.copy = value;
      btn.innerHTML = COPY_ICON_SVG;
      btn.hidden = false;
    } else {
      delete btn.dataset.copy;
      btn.hidden = true;
    }
  };
  setCopy("m-lis-copy",  m.LIS_ID);
  setCopy("m-name-copy", m.Name);
  setCopy("m-mrn-copy",  m.MRN);

  // EMR link is hospital-internal; only build it when MRN is present.
  const emr = document.getElementById("m-emr-link");
  if (m.MRN) {
    emr.href = `http://hisweb.hosp.ncku/Emrquery/autologin.aspx?chartno=${encodeURIComponent(m.MRN)}`;
    emr.hidden = false;
  } else {
    emr.removeAttribute("href"); emr.hidden = true;
  }
  // 🔄 EMR sync button only appears when the server has a client_id
  // configured AND the sample carries an MRN to look up.
  const sync = document.getElementById("btn-emr-sync");
  if (sync) sync.hidden = !(state.emrEnabled && m.MRN);

  // Editable selects backed by sample_metadata.json
  document.getElementById("m-test").value  = m.Test || "";
  document.getElementById("m-build").value = state.data.genome_build || "";
  document.getElementById("m-sex").value   = m.Sex || "";
  renderPloidySexStatus(m.Sex || "");
  document.getElementById("m-sry-confirmed").checked = !!state.reports.sry_confirmed;
  // 科別 / 開單醫師 — <select> populated from /api/patient_list/options
  // (cached). Reviewer picks an existing value or "＋ 新增…" which
  // prompts for a free-text value; the new value is added in-place
  // and persisted via _saveSampleMeta.
  _populateRosterOptions(m.Department || "", m.Physician || "");

  const sel = document.getElementById("m-category");
  const opts = (state.options && state.options.category_options) || [];
  const current = state.reports.category ?? m.Category ?? "";
  const all = Array.from(new Set(["", ...opts, current].filter(x => x !== undefined && x !== null)));
  sel.innerHTML = all.map(o => {
    const label = o === "" ? "—" : o;
    return `<option value="${escapeAttr(o)}" ${o === current ? "selected" : ""}>${escapeHtml(label)}</option>`;
  }).join("");

  document.getElementById("sample-card").classList.remove("hidden");
  applyDiagnosticAnalysisVisibility();
  renderQcWarnings();
}

// QC blacklist banner. The pipeline emits qc_summary.json with at most
// a top-level `blacklist` array of {gene, level, reason}; we hide the
// whole card when it's empty.
function renderQcWarnings() {
  const qc = state.data.qc_summary || {};
  const items = Array.isArray(qc.blacklist) ? qc.blacklist : [];
  const card = document.getElementById("qc-card");
  const ul   = document.getElementById("qc-warnings");
  if (!items.length) { card.classList.add("hidden"); ul.innerHTML = ""; return; }
  ul.innerHTML = items.map(w => `
    <li class="qc-warning qc-warning-${escapeAttr(w.level || "")}">
      <span class="qc-gene">${escapeHtml(w.gene || "?")}</span>
      <span class="qc-level">${escapeHtml(w.level || "")}</span>
      <span class="qc-reason">${escapeHtml(w.reason || "")}</span>
    </li>`).join("");
  card.classList.remove("hidden");
}

// Sample-metadata edit: save on change for Test / Build (Category goes
// via the legacy reports flow). Debounced so rapid keypresses don't
// produce a flurry of writes.
let _metaSaveTimer = null;
function _saveSampleMeta(patch) {
  if (!state.currentLIS) return;
  const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
  const sid = row?.sample_id || state.currentLIS;
  const hint = document.getElementById("m-meta-hint");
  clearTimeout(_metaSaveTimer);
  hint.textContent = "儲存中…";
  _metaSaveTimer = setTimeout(async () => {
    try {
      await apiPut(`/samples/${encodeURIComponent(sid)}/metadata`, patch);
      hint.textContent = `已儲存 ${new Date().toLocaleTimeString()}`;
    } catch (e) {
      hint.textContent = "儲存失敗：" + e.message;
    }
  }, 300);
}

document.addEventListener("change", ev => {
  if (ev.target.id === "m-test") {
    const meta = state.data?.meta || {};
    const testType = normalizeSampleTestType(
      ev.target.value,
      meta.LIS_ID || meta.lis_id || state.currentLIS || "",
    );
    ev.target.value = testType;
    if (state.data?.meta) state.data.meta.Test = testType;
    const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
    if (row) row.Test = testType;
    state.diagnosticAnalysisVisible = testType !== "TITAN-WGS";
    applyDiagnosticAnalysisVisibility();
    _saveSampleMeta({ test_type: testType });
  }
  if (ev.target.id === "m-build") _saveSampleMeta({ genome_build: ev.target.value });
  if (ev.target.id === "m-sex") {
    if (state.data.meta) state.data.meta.Sex = ev.target.value;
    renderPloidySexStatus(ev.target.value);
    _saveSampleMeta({ sex: ev.target.value });
  }
  if (ev.target.id === "m-department") _onRosterSelectChange("department", "科別");
  if (ev.target.id === "m-physician")  _onRosterSelectChange("physician",  "開單醫師");
});

function _onRosterSelectChange(field, label) {
  const sel = document.getElementById(`m-${field}`);
  if (!sel) return;
  if (sel.value === "__new__") {
    const v = (prompt(`新增${label}：`, "") || "").trim();
    if (!v) {
      // Cancelled — revert to whatever was set when the sample loaded.
      const m = state.data?.meta || {};
      sel.value = field === "department" ? (m.Department || "") : (m.Physician || "");
      return;
    }
    // Add the new value to the cached options + the live <select> so
    // it's pickable next time without a roundtrip.
    const bucket = field === "department" ? "departments" : "physicians";
    _rosterOptions = _rosterOptions || { departments: [], physicians: [] };
    if (!_rosterOptions[bucket].includes(v)) {
      _rosterOptions[bucket] = [..._rosterOptions[bucket], v].sort();
    }
    const opt = document.createElement("option");
    opt.value = v; opt.text = v;
    sel.insertBefore(opt, sel.querySelector('option[value="__new__"]'));
    sel.value = v;
    _saveSampleMeta({ [field]: v });
  } else {
    _saveSampleMeta({ [field]: sel.value });
  }
}

// /api/patient_list/options returns {departments, physicians} unioned
// across the roster. Cache on first fetch — refresh on each upload
// (handled by setting _rosterOptions = null in the upload flow).
let _rosterOptions = null;
async function _populateRosterOptions(curDept, curPhys) {
  const deps = document.getElementById("m-department");
  const docs = document.getElementById("m-physician");
  if (!deps || !docs) return;
  if (_rosterOptions == null) {
    try { _rosterOptions = await apiFetch("/patient_list/options") || {}; }
    catch { _rosterOptions = {}; }
  }
  const fill = (sel, opts, current) => {
    // Always include the current value (even if it isn't in roster) so
    // a manually-entered value survives a re-render.
    const seen = new Set();
    const items = ["", ...(opts || [])];
    if (current && !items.includes(current)) items.push(current);
    sel.innerHTML = items.map(v => {
      if (seen.has(v)) return "";
      seen.add(v);
      const label = v === "" ? "—" : escapeHtml(v);
      const sel_ = v === current ? "selected" : "";
      return `<option value="${escapeAttr(v)}" ${sel_}>${label}</option>`;
    }).join("") + `<option value="__new__">＋ 新增…</option>`;
  };
  fill(deps, _rosterOptions.departments, curDept || "");
  fill(docs, _rosterOptions.physicians,  curPhys || "");
}

// Generic renderer for collapsible free-text cards. User-toggled state is
// remembered across re-renders via toggledBlocks.
function renderCollapsibleCard(cardId, headerId, bodyId, taId, value, defaultOpen = false) {
  const card   = document.getElementById(cardId);
  const header = document.getElementById(headerId);
  const body   = document.getElementById(bodyId);
  const ta     = document.getElementById(taId);

  ta.value = value || "";
  const open = toggledBlocks.has(cardId)
    ? card.dataset.wasOpen === "1"
    : defaultOpen;
  card.dataset.wasOpen = open ? "1" : "0";
  header.classList.toggle("open", open);
  body.classList.toggle("open", open);
  card.classList.remove("hidden");
  // After the body becomes display:block, run autoGrow so the
  // textarea matches the loaded content. Doing this synchronously
  // while the body is still display:none would yield scrollHeight=0.
  if (open && (taId === "clinical-text" || taId === "counseling-text")) {
    requestAnimationFrame(() => autoGrow(ta));
  }
}

function renderClinicalDescription() {
  renderCollapsibleCard("clinical-card", "clinical-header", "clinical-body",
                        "clinical-text", state.reports.clinical_description);
}

function renderGeneticCounseling() {
  // Counseling text lives in state.reports.genetic_counseling (also
  // mirrored in state.data.genetic_counseling on load). The header
  // shows the last EMR sync timestamp so reviewers know whether the
  // text is auto-pulled or hand-edited.
  const value = state.reports.genetic_counseling
              ?? state.data.genetic_counseling ?? "";
  renderCollapsibleCard("counseling-card", "counseling-header", "counseling-body",
                        "counseling-text", value);
  const syncedEl = document.getElementById("counseling-synced");
  const synced = state.data.emr_synced_at || "";
  if (syncedEl) {
    syncedEl.textContent = synced ? `EMR synced: ${synced.slice(0, 10)}` : "";
  }
}

function renderComment() {
  renderCollapsibleCard("comment-card", "comment-header", "comment-body",
                        "comment-text", state.reports.comment, true);
  renderTagPicker();
}

function renderDeadZoneCard() {
  const card = document.getElementById("dead-zone-card");
  const body = document.getElementById("dead-zone-body");
  const thrEl = document.getElementById("dead-zone-threshold");
  const headerToggle = document.getElementById("dead-zone-header-toggle");
  if (!card || !body) return;
  const dz = state.data?.dead_zone || {};
  const threshold = dz.threshold || "";
  const entries = sortDeadZoneEntries(Array.isArray(dz.entries) ? dz.entries : []);
  if (thrEl) {
    thrEl.textContent = threshold ? `clinical threshold: ${threshold}X` : "";
  }
  card.classList.remove("hidden");
  if (!entries.length) {
    body.innerHTML = `<div class="muted">目前 HPO / panel 基因沒有 cohort dead-zone 註記。</div>`;
    body.dataset.expanded = "0";
    body.dataset.deadZoneKey = "";
    if (headerToggle) {
      headerToggle.classList.add("hidden");
      headerToggle.textContent = "";
      headerToggle.onclick = null;
    }
    return;
  }
  const key = entries.map(e => `${e.gene || ""}:${e.exons_label || ""}:${e.cds_dead_pct || ""}`).join("|");
  if (body.dataset.deadZoneKey !== key) {
    body.dataset.deadZoneKey = key;
    body.dataset.expanded = "0";
  }
  const expanded = body.dataset.expanded === "1";
  const collapsedEntries = entries.filter(e => Number(e?.cds_dead_pct || 0) >= 50);
  const visibleEntries = expanded ? entries : collapsedEntries;
  const collapsedHiddenCount = Math.max(0, entries.length - collapsedEntries.length);
  const hasToggle = collapsedHiddenCount > 0;
  const toggleLabel = expanded ? "收合" : `展開全部（另 ${collapsedHiddenCount} 列）`;
  const toggleHtml = hasToggle ? `
    <button id="dead-zone-toggle" class="dead-zone-toggle" type="button" aria-expanded="${expanded ? "true" : "false"}">
      <span class="dead-zone-toggle-arrow">${expanded ? "▾" : "▸"}</span>
      ${toggleLabel}
    </button>` : "";
  if (headerToggle) {
    headerToggle.classList.toggle("hidden", !hasToggle);
    headerToggle.setAttribute("aria-expanded", expanded ? "true" : "false");
    headerToggle.innerHTML = hasToggle
      ? `<span class="dead-zone-toggle-arrow">${expanded ? "▾" : "▸"}</span>${toggleLabel}`
      : "";
    headerToggle.onclick = null;
  }
  body.innerHTML = `<ul class="dead-zone-list">${visibleEntries.map(e => {
    const gene = e.gene || "";
    const label = e.exons_label || (Array.isArray(e.exons) ? e.exons.join(", ") : "");
    const pct = Number(e.cds_dead_pct || 0);
    const pctClass = deadZonePctClass(pct);
    const pctLabel = Number.isFinite(pct) ? `${pct.toFixed(1).replace(/\\.0$/, "")}%` : "";
    const pheno = Number(e.pheno_score || 0);
    const phenoLabel = Number.isFinite(pheno) && pheno > 0
      ? `<span class="dead-zone-pheno">Pheno ${pheno.toFixed(1).replace(/\\.0$/, "")}</span>`
      : "";
    return `<li class="${pctClass}"><span class="dead-zone-gene">${escapeHtml(gene)}</span>
      <span class="dead-zone-exons">exon ${escapeHtml(label)}</span>
      ${pctLabel ? `<span class="dead-zone-cds">CDS ${escapeHtml(pctLabel)}</span>` : ""}
      ${phenoLabel}</li>`;
  }).join("")}</ul>${!expanded && !visibleEntries.length && entries.length
    ? `<div class="muted">沒有 CDS ≥50% 的 dead-zone；可展開全部查看其他 ${entries.length} 列。</div>`
    : ""}${toggleHtml}`;
  const toggleExpanded = () => {
    body.dataset.expanded = expanded ? "0" : "1";
    renderDeadZoneCard();
  };
  document.getElementById("dead-zone-toggle")?.addEventListener("click", toggleExpanded);
  if (headerToggle && hasToggle) headerToggle.onclick = toggleExpanded;
}

function deadZoneBucket(pct) {
  const n = Number(pct || 0);
  if (n >= 70) return 0;
  if (n >= 50) return 1;
  if (n >= 30) return 2;
  return 3;
}

function deadZonePctClass(pct) {
  const n = Number(pct || 0);
  if (n >= 70) return "dead-zone-pct-high";
  if (n >= 50) return "dead-zone-pct-mid";
  if (n >= 30) return "dead-zone-pct-low";
  return "dead-zone-pct-base";
}

function sortDeadZoneEntries(entries) {
  return (entries || []).slice().sort((a, b) => {
    const pctA = Number(a?.cds_dead_pct || 0);
    const pctB = Number(b?.cds_dead_pct || 0);
    const phenoA = Number(a?.pheno_score || 0);
    const phenoB = Number(b?.pheno_score || 0);
    return (deadZoneBucket(pctA) - deadZoneBucket(pctB))
      || (phenoB - phenoA)
      || (pctB - pctA)
      || String(a?.gene || "").localeCompare(String(b?.gene || ""));
  });
}

function _reportedOutOfDiseaseAssociatedSnvs() {
  const variants = state.data?.variants || {};
  const status = state.reports?.status || {};
  const rows = [];
  for (const [id, rawStatus] of Object.entries(status)) {
    const vals = _statusValues(rawStatus).filter(v => v === "1" || v === "2" || v === "C");
    if (!vals.length) continue;
    const v = variants[id];
    if (!v || v.disease_associated) continue;
    rows.push({ id, gene: v.gene_symbol || "?", status: vals.join("/") });
  }
  rows.sort((a, b) => a.gene.localeCompare(b.gene) || a.id.localeCompare(b.id));
  return rows;
}

function renderDiseaseAssociatedReportWarning() {
  const el = document.getElementById("disease-associated-warning");
  if (!el) return;
  const rows = _reportedOutOfDiseaseAssociatedSnvs();
  el.classList.toggle("hidden", rows.length === 0);
  if (!rows.length) {
    el.innerHTML = "";
    return;
  }
  const genes = Array.from(new Set(rows.map(r => r.gene))).join(", ");
  el.innerHTML = `<strong>提醒：</strong>已標記的 SNV/Indel 有基因不在 disease-associated gene list，DOCX 的 §五.4 基因清單不會列出這些基因：${escapeHtml(genes)}`;
}

// Known tag suggestions = tags pulled from the Tag column of every loaded
// NGS_list row, plus anything the user has added during this session.
const sessionTags = new Set();
function getAllKnownTags() {
  const set = new Set(sessionTags);
  if (Array.isArray(state.index)) {
    for (const r of state.index) {
      const raw = r.Tag ?? r.tag ?? "";
      if (!raw) continue;
      String(raw).split(/[,;]\s*/).forEach(t => {
        const v = t.trim();
        if (v) set.add(v);
      });
    }
  }
  return Array.from(set).sort((a, b) => a.localeCompare(b));
}

function renderTagPicker() {
  const wrap = document.getElementById("tag-picker");
  if (!wrap) return;
  const tags = state.reports.tags || [];
  // Track selected ones in sessionTags so they remain auto-completable
  // when the user moves between samples without saving in between.
  tags.forEach(t => sessionTags.add(t));
  const dlOpts = getAllKnownTags()
    .filter(t => !tags.includes(t))
    .map(t => `<option value="${escapeAttr(t)}"></option>`)
    .join("");
  wrap.innerHTML = `
    <div class="tag-label">Tag</div>
    <div class="tag-chips">
      ${tags.map(t => `
        <span class="tag-chip">${escapeHtml(t)}<button class="tag-remove" data-tag="${escapeAttr(t)}" type="button" title="移除">×</button></span>
      `).join("")}
      <input class="tag-input" list="tag-options-dl" placeholder="新增…" autocomplete="off" />
    </div>
    <datalist id="tag-options-dl">${dlOpts}</datalist>
  `;
}

function addTag(value) {
  const v = String(value || "").trim();
  if (!v) return;
  if (!Array.isArray(state.reports.tags)) state.reports.tags = [];
  if (state.reports.tags.includes(v)) return;
  state.reports.tags.push(v);
  sessionTags.add(v);
  state.dirty = true;
  renderTagPicker();
  updateSaveHint();
}

function removeTag(value) {
  if (!Array.isArray(state.reports.tags)) return;
  state.reports.tags = state.reports.tags.filter(x => x !== value);
  state.dirty = true;
  renderTagPicker();
  updateSaveHint();
}

// ---------- Phenotype + panel editor (Phase A/B) ------------------

// Editable working copy. We don't mutate state.data.patient_phenotype
// directly so a "reload sample" cleanly resets the form to whatever
// the server currently has on disk.
const phenoEdit = {
  hpo:    [],   // [{phenotype, label, weight}]
  panels: [],   // [panel_name]
};

function renderPhenotype() {
  // Seed working copy from sample payload. Panels persist as {name, weight}
  // dicts; legacy server payloads where it was a flat list of strings are
  // upgraded to weight=1.
  phenoEdit.hpo = (state.data.patient_phenotype || []).map(r => ({
    phenotype: r.phenotype || "",
    label:     r.label || "",
    weight:    Number.isFinite(Number(r.weight)) ? Number(r.weight) : 1,
  }));
  const rawPanels = Array.isArray(state.data.selected_panels)
    ? state.data.selected_panels
    : [];
  phenoEdit.panels = rawPanels.map(p => typeof p === "string"
    ? { name: p, weight: 1 }
    : { name: p.name, weight: Number(p.weight) || 1 });

  renderHpoChips();
  renderPanelChips();
  document.getElementById("phenotype-stats").textContent = "";
  document.getElementById("phenotype-hint").textContent  = "";
  document.getElementById("phenotype-top10").classList.add("hidden");
  document.getElementById("phenotype-card").classList.remove("hidden");
  _initPhenoPanelTabs();
  renderFixedPanelHosts();  // async, but the chip area updates on its own.
  // If there's a queued/running job for this sample, pick up polling.
  _resumeJobPollingIfAny();
}

async function _resumeJobPollingIfAny() {
  if (!state.currentLIS) return;
  clearInterval(_jobPollTimer);
  _setJobStatus("");
  if (!Array.isArray(phenoEdit.hpo) || !phenoEdit.hpo.length) return;
  const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
  const sid = row?.sample_id || state.currentLIS;
  try {
    const jobs = await apiFetch(`/samples/${encodeURIComponent(sid)}/jobs`) || [];
    const live = jobs.find(j => j.status === "queued" || j.status === "running");
    if (live) {
      _activeJobId = live.job_id;
      const tool = _stepTool(live.step);
      _setJobStatus(tool ? `${live.status} ${tool}` : live.status, true);
      _startJobPolling(sid, live.job_id);
    } else if (jobs.length) {
      _setJobStatus(jobs[0].status);
    }
  } catch (e) { /* ignore */ }
}

function renderHpoChips() {
  const ul = document.getElementById("phenotype-list");
  ul.innerHTML = "";
  phenoEdit.hpo.forEach((row, idx) => {
    const li = document.createElement("li");
    li.className = "chip chip-hpo";
    li.innerHTML = `
      <span class="chip-label">${escapeHtml(row.label || row.phenotype)}</span>
      <span class="chip-id">${escapeHtml(row.phenotype)}</span>
      <select class="chip-weight" data-idx="${idx}" title="Weight">
        ${[1,2,3,4,5].map(n => `<option value="${n}" ${n===row.weight?"selected":""}>w=${n}</option>`).join("")}
      </select>
      <button class="chip-remove" data-idx="${idx}" type="button" title="移除">×</button>`;
    ul.appendChild(li);
  });
}

function renderPanelChips() {
  const ul = document.getElementById("panel-chips");
  ul.innerHTML = "";
  phenoEdit.panels.forEach((row, idx) => {
    const li = document.createElement("li");
    li.className = "chip chip-panel";
    // Fixed-panel keys are wide; show the pretty "WES-I · 皮膚科 · EB"
    // form when we know about it. Free-text panels render as-is.
    const label = _fixedPanelMeta.has(row.name)
      ? _fixedPanelDisplayName(row.name) : row.name;
    li.innerHTML = `
      <span class="chip-label" title="${escapeAttr(row.name)}">${escapeHtml(label)}</span>
      <select class="chip-weight" data-panel-idx="${idx}" title="Weight">
        ${[1,2,3,4,5].map(n => `<option value="${n}" ${n===row.weight?"selected":""}>w=${n}</option>`).join("")}
      </select>
      <button class="chip-remove" data-panel-idx="${idx}" type="button" title="移除">×</button>`;
    ul.appendChild(li);
  });
  // Fixed-panel chip-checkbox state mirrors phenoEdit.panels — keep in sync.
  syncFixedPanelChipState();
}

// Cached panel list from /api/panels — fetched once per session.
let _panelOptions = null;
async function loadPanelOptions() {
  if (_panelOptions) return _panelOptions;
  _panelOptions = await apiFetch("/panels") || [];
  return _panelOptions;
}

// ── Fixed gene-panel tabs (WES-I / WES-II / WGS / Other) ─────────
// Index file written by scripts/import_fixed_panels.py. Each fixed
// panel's `key` matches a file in GENE_PANELS_DIR, so phenotype_scorer
// consumes them via the same `phenoEdit.panels[].name` path as
// reviewer-typed panels — toggling a chip just adds/removes the key
// from phenoEdit.panels with weight=1.
let _fixedPanelIndex = null;        // { series: [{key, label, groups: [...]}] }
const _fixedPanelKeys = new Set();
const _fixedPanelMeta = new Map();  // key → {series, category, name}

async function loadFixedPanelIndex() {
  if (_fixedPanelIndex) return _fixedPanelIndex;
  try {
    const resp = await fetch(`${API_BASE}/phenotype-tool/fixed-panels`,
                             { credentials: "same-origin" });
    _fixedPanelIndex = resp.ok ? await resp.json() : { series: [] };
  } catch { _fixedPanelIndex = { series: [] }; }
  _fixedPanelKeys.clear();
  _fixedPanelMeta.clear();
  for (const s of (_fixedPanelIndex.series || [])) {
    for (const g of (s.groups || [])) {
      for (const p of (g.panels || [])) {
        _fixedPanelKeys.add(p.key);
        _fixedPanelMeta.set(p.key, { series: s.key, category: g.category, name: p.name });
      }
    }
  }
  return _fixedPanelIndex;
}

function _initPhenoPanelTabs() {
  document.querySelectorAll("#phenotype-card .panel-tab").forEach((btn) => {
    if (btn.dataset.wired === "1") return;
    btn.dataset.wired = "1";
    btn.addEventListener("click", () => {
      const target = btn.dataset.phenoTab;
      const isOpen = btn.classList.contains("is-active");
      // Toggle: click the open tab to collapse, click any tab to switch.
      document.querySelectorAll("#phenotype-card .panel-tab").forEach((b) =>
        b.classList.toggle("is-active", !isOpen && b === btn));
      document.querySelectorAll("#phenotype-card .panel-tab-body").forEach((body) =>
        body.classList.toggle("is-active", !isOpen && body.dataset.phenoTabBody === target));
    });
  });
}

function _fixedPanelHostHtml(series) {
  if (!series || !(series.groups || []).length) {
    return '<div class="muted">尚未匯入此系列的 panel</div>';
  }
  return series.groups.map((g) => `
    <div class="fp-group">
      <div class="fp-group-title">${escapeHtml(g.category)}</div>
      <div class="fp-chips">
        ${(g.panels || []).map((p) => `
          <label class="fp-chip" data-key="${escapeAttr(p.key)}" title="${escapeAttr(p.key)}">
            <input type="checkbox" class="fp-chip-cb" value="${escapeAttr(p.key)}">
            <span class="fp-chip-label">${escapeHtml(p.name)}</span>
            <span class="fp-chip-count">(${p.gene_count || 0})</span>
          </label>
        `).join("")}
      </div>
    </div>
  `).join("");
}

async function renderFixedPanelHosts() {
  await loadFixedPanelIndex();
  const seriesByKey = {};
  for (const s of (_fixedPanelIndex.series || [])) seriesByKey[s.key] = s;
  document.querySelectorAll("#phenotype-card .fixed-panel-host").forEach((host) => {
    host.innerHTML = _fixedPanelHostHtml(seriesByKey[host.dataset.series]);
  });
  document.querySelectorAll("#phenotype-card .fp-chip-cb").forEach((cb) => {
    cb.addEventListener("change", () => {
      const key = cb.value;
      if (cb.checked) addPanel(key);
      else {
        const idx = phenoEdit.panels.findIndex(p => p.name === key);
        if (idx >= 0) removePanel(idx);
      }
      // renderPanelChips already called inside add/remove — re-sync chip state.
      syncFixedPanelChipState();
    });
  });
  syncFixedPanelChipState();
  // renderPhenotype() may have drawn selected chips before the async
  // fixed-panel index arrived. Refresh them now so reviewer-facing
  // labels use the corrected display name instead of the stored key.
  renderPanelChips();
}

function syncFixedPanelChipState() {
  const picked = new Set(phenoEdit.panels.map(p => p.name));
  document.querySelectorAll("#phenotype-card .fp-chip-cb").forEach((cb) => {
    const on = picked.has(cb.value);
    cb.checked = on;
    cb.closest(".fp-chip").classList.toggle("is-selected", on);
  });
}

function _initNewCasePanelTabs() {
  document.querySelectorAll("#new-case-form .panel-tab").forEach((btn) => {
    if (btn.dataset.wired === "1") return;
    btn.dataset.wired = "1";
    btn.addEventListener("click", () => {
      const target = btn.dataset.ncPhenoTab;
      const isOpen = btn.classList.contains("is-active");
      document.querySelectorAll("#new-case-form .panel-tab").forEach((b) =>
        b.classList.toggle("is-active", !isOpen && b === btn));
      document.querySelectorAll("#new-case-form .panel-tab-body").forEach((body) =>
        body.classList.toggle("is-active", !isOpen && body.dataset.ncPhenoTabBody === target));
    });
  });
}

function _resetNewCasePanelTabs() {
  document.querySelectorAll("#new-case-form .panel-tab").forEach((btn) =>
    btn.classList.toggle("is-active", btn.dataset.ncPhenoTab === "other"));
  document.querySelectorAll("#new-case-form .panel-tab-body").forEach((body) =>
    body.classList.toggle("is-active", body.dataset.ncPhenoTabBody === "other"));
}

async function renderNewCaseFixedPanelHosts() {
  await loadFixedPanelIndex();
  const seriesByKey = {};
  for (const s of (_fixedPanelIndex.series || [])) seriesByKey[s.key] = s;
  document.querySelectorAll("#new-case-form .fixed-panel-host").forEach((host) => {
    host.innerHTML = _fixedPanelHostHtml(seriesByKey[host.dataset.series]);
  });
  document.querySelectorAll("#new-case-form .fp-chip-cb").forEach((cb) => {
    cb.addEventListener("change", () => {
      const key = cb.value;
      const idx = newCaseEdit.panels.findIndex(p => p.name === key);
      if (cb.checked && idx < 0) newCaseEdit.panels.push({ name: key, weight: 1 });
      if (!cb.checked && idx >= 0) newCaseEdit.panels.splice(idx, 1);
      _markNewCaseEdited();
      renderNewCasePhenoEditor();
    });
  });
  syncNewCaseFixedPanelChipState();
  // Existing phenotype files can seed selected panels before the
  // index request completes. Re-render once metadata is available.
  renderNewCasePhenoEditor();
}

function syncNewCaseFixedPanelChipState() {
  const picked = new Set((newCaseEdit.panels || []).map(p => p.name));
  document.querySelectorAll("#new-case-form .fp-chip-cb").forEach((cb) => {
    const on = picked.has(cb.value);
    cb.checked = on;
    cb.closest(".fp-chip").classList.toggle("is-selected", on);
  });
}

function _fixedPanelDisplayName(key) {
  const m = _fixedPanelMeta.get(key);
  return m ? `${m.series} · ${m.category} · ${m.name}` : key;
}

// HPO + panel typeaheads use delegated listeners on `document` so they
// keep working even if the phenotype-card is re-rendered. Per-element
// addEventListener was racing the legacy global handlers; delegated
// dispatch sidesteps it entirely.
let _hpoSearchAbort = null;
let _hpoSearchTimer = null;
const _comboActive = new WeakMap();

function _hpoOpen()  { document.getElementById("hpo-search-dropdown")?.classList.remove("hidden"); }
function _hpoClose() {
  const drop = document.getElementById("hpo-search-dropdown");
  drop?.classList.add("hidden");
  if (drop) _comboClearActive(drop);
}
function _panelOpen()  { document.getElementById("panel-search-dropdown")?.classList.remove("hidden"); }
function _panelClose() {
  const drop = document.getElementById("panel-search-dropdown");
  drop?.classList.add("hidden");
  if (drop) _comboClearActive(drop);
}

function _comboOptions(dropdown) {
  return Array.from(dropdown?.querySelectorAll(".combobox-option") || []);
}

function _comboSetActive(dropdown, idx) {
  if (!dropdown) return;
  const opts = _comboOptions(dropdown);
  opts.forEach((opt, i) => opt.classList.toggle("active", i === idx));
  if (idx >= 0 && opts[idx]) {
    _comboActive.set(dropdown, idx);
    opts[idx].scrollIntoView({ block: "nearest" });
  } else {
    _comboActive.delete(dropdown);
  }
}

function _comboClearActive(dropdown) {
  _comboSetActive(dropdown, -1);
}

function _comboMove(dropdown, delta) {
  if (!dropdown) return;
  const opts = _comboOptions(dropdown);
  if (!opts.length) return;
  const cur = _comboActive.has(dropdown) ? _comboActive.get(dropdown) : -1;
  const next = Math.max(0, Math.min(opts.length - 1, cur + delta));
  _comboSetActive(dropdown, next);
}

function _pickHpoOption(opt) {
  if (!opt) return false;
  addHpo(opt.dataset.id, opt.dataset.name);
  const inp = document.getElementById("hpo-search");
  if (inp) inp.value = "";
  _hpoClose();
  return true;
}

function _pickPanelOption(opt) {
  if (!opt) return false;
  addPanel(opt.dataset.name);
  const inp = document.getElementById("panel-search");
  if (inp) inp.value = "";
  _panelClose();
  return true;
}

function _pickNewCaseHpoOption(opt) {
  if (!opt?.dataset.ncHpoPick) return false;
  const r = JSON.parse(opt.dataset.ncHpoPick);
  const id = r.hpo_id || r.phenotype || "";
  if (!id) return false;
  if (!newCaseEdit.hpo.some(h => h.phenotype === id)) {
    newCaseEdit.hpo.push({phenotype: id, label: r.name || id, weight: 1});
    _markNewCaseEdited();
  }
  document.getElementById("new-case-hpo-search").value = "";
  document.getElementById("new-case-hpo-search-dropdown").classList.add("hidden");
  _comboClearActive(document.getElementById("new-case-hpo-search-dropdown"));
  renderNewCasePhenoEditor();
  return true;
}

function _pickNewCasePanelOption(opt) {
  if (!opt?.dataset.ncPanelPick) return false;
  const r = JSON.parse(opt.dataset.ncPanelPick);
  const name = r.name || "";
  if (!name) return false;
  if (!newCaseEdit.panels.some(p => p.name === name)) {
    newCaseEdit.panels.push({name, weight: 1});
    _markNewCaseEdited();
  }
  document.getElementById("new-case-panel-search").value = "";
  document.getElementById("new-case-panel-search-dropdown").classList.add("hidden");
  _comboClearActive(document.getElementById("new-case-panel-search-dropdown"));
  renderNewCasePhenoEditor();
  return true;
}

function _comboPickActive(dropdown, fallbackFirst = true) {
  const opts = _comboOptions(dropdown);
  if (!opts.length) return false;
  const idx = _comboActive.has(dropdown) ? _comboActive.get(dropdown) : (fallbackFirst ? 0 : -1);
  const opt = opts[idx];
  if (!opt) return false;
  if (dropdown.id === "hpo-search-dropdown") return _pickHpoOption(opt);
  if (dropdown.id === "panel-search-dropdown") return _pickPanelOption(opt);
  if (dropdown.id === "new-case-lis-id-dropdown") return _pickNewCaseLisOption(opt);
  if (dropdown.id === "new-case-hpo-search-dropdown") return _pickNewCaseHpoOption(opt);
  if (dropdown.id === "new-case-panel-search-dropdown") return _pickNewCasePanelOption(opt);
  return false;
}

function _handleComboKeydown(ev, dropdown) {
  if (!dropdown || dropdown.classList.contains("hidden")) return false;
  if (ev.key === "ArrowDown") {
    ev.preventDefault();
    ev.stopPropagation();
    _comboMove(dropdown, 1);
    return true;
  }
  if (ev.key === "ArrowUp") {
    ev.preventDefault();
    ev.stopPropagation();
    _comboMove(dropdown, -1);
    return true;
  }
  if (ev.key === "Enter") {
    ev.preventDefault();
    ev.stopPropagation();
    _comboPickActive(dropdown, true);
    return true;
  }
  if (ev.key === "Escape") {
    ev.preventDefault();
    ev.stopPropagation();
    dropdown.classList.add("hidden");
    _comboClearActive(dropdown);
    return true;
  }
  return false;
}

function setupHpoSearchInput() {
  // No-op now; kept so boot() doesn't throw if the call site is still here.
  // Real work lives in the document-level handler at the bottom of this
  // file (search for "Phase B delegated typeahead").
}

function setupPanelSearchInput() {
  // No-op; see setupHpoSearchInput().
}

async function _runHpoSearch(q) {
  const dropdown = document.getElementById("hpo-search-dropdown");
  if (!dropdown) return;
  if (_hpoSearchAbort) _hpoSearchAbort.abort();
  _hpoSearchAbort = new AbortController();
  try {
    const url  = `${API_BASE}/hpo/search?q=${encodeURIComponent(q)}&limit=20`;
    const resp = await fetch(url, { signal: _hpoSearchAbort.signal });
    if (!resp.ok) { _hpoClose(); return; }
    const list = await resp.json();
    if (!Array.isArray(list) || !list.length) {
      dropdown.innerHTML = '<li class="muted" style="padding:6px 10px">（無結果）</li>';
    } else {
      dropdown.innerHTML = list.map(t => {
        // "12 genes" for an annotated HPO; nothing when 0 (the dim grey
        // ones are HPOs without any gene annotation, e.g. parent terms).
        const gc = Number.isFinite(Number(t.gene_count)) && t.gene_count > 0
          ? `${t.gene_count} genes`
          : "";
        return `
        <li class="combobox-option" data-id="${escapeAttr(t.hpo_id)}" data-name="${escapeAttr(t.name)}">
          <span class="opt-lis">${escapeHtml(t.hpo_id)}</span>
          <span class="opt-name">${escapeHtml(t.name)}</span>
          <span class="opt-mrn">${escapeHtml(gc)}</span>
        </li>`;
      }).join("");
    }
    _comboClearActive(dropdown);
    _hpoOpen();
  } catch (e) {
    if (e.name !== "AbortError") console.error("HPO search failed", e);
  }
}

async function _runPanelSearch(q) {
  const dropdown = document.getElementById("panel-search-dropdown");
  if (!dropdown) return;
  const opts = await loadPanelOptions();
  await loadFixedPanelIndex();   // populates _fixedPanelKeys
  const ql = (q || "").trim().toLowerCase();
  const picked = new Set(phenoEdit.panels.map(p => p.name));
  const matches = opts
    .filter(p => !picked.has(p.name) && !_fixedPanelKeys.has(p.name)
                 && (!ql || p.name.toLowerCase().includes(ql)))
    .slice(0, 30);
  if (!matches.length) {
    dropdown.innerHTML = '<li class="muted" style="padding:6px 10px">（無結果）</li>';
  } else {
    dropdown.innerHTML = matches.map(p => `
      <li class="combobox-option" data-name="${escapeAttr(p.name)}">
        <span class="opt-lis">${escapeHtml(p.name)}</span>
        <span class="opt-name">${p.gene_count} genes</span>
        <span class="opt-mrn"></span>
      </li>`).join("");
  }
  _comboClearActive(dropdown);
  _panelOpen();
}

// Phase B delegated typeahead: catch input/focus/blur on the two
// search boxes and clicks on their dropdown options at the document
// level so we don't depend on element-specific listeners.
document.addEventListener("input", ev => {
  const t = ev.target;
  if (t.id === "hpo-search") {
    clearTimeout(_hpoSearchTimer);
    const q = t.value.trim();
    if (!q) { _hpoClose(); return; }
    _hpoSearchTimer = setTimeout(() => _runHpoSearch(q), 200);
  } else if (t.id === "panel-search") {
    _runPanelSearch(t.value);
  }
});

document.addEventListener("focusin", ev => {
  if (ev.target.id === "hpo-search") {
    if (ev.target.value.trim()) _runHpoSearch(ev.target.value.trim());
  } else if (ev.target.id === "panel-search") {
    _runPanelSearch(ev.target.value);
  }
});

document.addEventListener("focusout", ev => {
  if (ev.target.id === "hpo-search")     setTimeout(_hpoClose,   150);
  if (ev.target.id === "panel-search")   setTimeout(_panelClose, 150);
});

document.addEventListener("keydown", ev => {
  if (ev.target.id === "hpo-search") {
    _handleComboKeydown(ev, document.getElementById("hpo-search-dropdown"));
  } else if (ev.target.id === "panel-search") {
    _handleComboKeydown(ev, document.getElementById("panel-search-dropdown"));
  } else if (ev.target.id === "new-case-lis-id-search") {
    _handleComboKeydown(ev, document.getElementById("new-case-lis-id-dropdown"));
  } else if (ev.target.id === "new-case-hpo-search") {
    _handleComboKeydown(ev, document.getElementById("new-case-hpo-search-dropdown"));
  } else if (ev.target.id === "new-case-panel-search") {
    _handleComboKeydown(ev, document.getElementById("new-case-panel-search-dropdown"));
  }
});

// Mousedown so the option fires before the input's blur kills the dropdown.
document.addEventListener("mousedown", ev => {
  const opt = ev.target.closest(".combobox-option");
  if (!opt) return;
  const dropdown = opt.parentElement;
  if (dropdown?.id === "hpo-search-dropdown") {
    ev.preventDefault();
    _pickHpoOption(opt);
  } else if (dropdown?.id === "panel-search-dropdown") {
    ev.preventDefault();
    _pickPanelOption(opt);
  }
});

function addHpo(id, label) {
  if (!id) return;
  if (phenoEdit.hpo.some(r => r.phenotype === id)) return;
  // Default weight = 1; user adjusts via the chip's own select.
  phenoEdit.hpo.push({ phenotype: id, label: label || id, weight: 1 });
  renderHpoChips();
}

function removeHpo(idx) {
  phenoEdit.hpo.splice(idx, 1);
  renderHpoChips();
}

function setHpoWeight(idx, weight) {
  if (phenoEdit.hpo[idx]) {
    phenoEdit.hpo[idx].weight = Number(weight) || 1;
  }
}

function addPanel(name) {
  if (!name || phenoEdit.panels.some(p => p.name === name)) return;
  phenoEdit.panels.push({ name, weight: 1 });
  renderPanelChips();
}

function removePanel(idx) {
  phenoEdit.panels.splice(idx, 1);
  renderPanelChips();
}

function setPanelWeight(idx, weight) {
  if (phenoEdit.panels[idx]) {
    phenoEdit.panels[idx].weight = Number(weight) || 1;
  }
}

async function apiPost(path, body) {
  const resp = await fetch(`${API_BASE}${path}`, {
    method: "POST",
    credentials: "same-origin",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (resp.status === 401) { showLoginModal(); throw new Error("not authenticated"); }
  if (!resp.ok) throw new Error(`${resp.status} ${resp.statusText} on ${path}`);
  return await resp.json();
}

// ---------- Auth flow ----------------------------------------------

function showLoginModal(msg) {
  const m = document.getElementById("login-modal");
  if (!m) return;
  const err = document.getElementById("login-error");
  if (msg) { err.textContent = msg; err.classList.remove("hidden"); }
  else     { err.classList.add("hidden"); }
  m.classList.remove("hidden");
  document.getElementById("login-username").focus();
}
function hideLoginModal() {
  document.getElementById("login-modal")?.classList.add("hidden");
}

// ---------- Gene search modal -------------------------------------
//
// The SNV / Indel and CNV / SV card headers each carry a gene search
// box. Typing one or more symbols + Enter opens a modal listing every
// matching variant as cards — the same renderers (renderVariantCard /
// renderCnvSvCard) the tier tables use, so the cards are fully
// interactive (status dropdown, disease list, comment, …). The
// modal's own input lets the reviewer pivot to another gene without
// closing it.

function _parseGeneSearch(raw) {
  return Array.from(new Set(
    String(raw || "").toUpperCase().split(/[,，、或\s]+/)
      .map(g => g.trim()).filter(Boolean)
  ));
}

function _numericValue(raw) {
  const s = String(raw ?? "").trim();
  if (!s || s === "." || s.toUpperCase() === "NA" || s.toUpperCase() === "N/A") return null;
  const n = Number(s);
  return Number.isFinite(n) ? n : null;
}

function _passesGnomadAfFilter(v) {
  const af = _numericValue(v?.AF);
  return af == null || af < 0.01;
}

function _isReferenceZygosity(v) {
  return ["ref", "hom_ref", "hom-ref", "0/0", "0|0"].includes(
    String(v?.zygosity || "").trim().toLowerCase()
  );
}

async function _geneSearchSnv(genesUpper, { filterGnomad = true } = {}) {
  const sid = state.data?.sample_id || state.currentLIS;
  if (sid && genesUpper.length) {
    const params = new URLSearchParams({ genes: genesUpper.join(",") });
    const payload = await apiFetch(`/samples/${encodeURIComponent(sid)}/snv-search?${params}`);
    if (payload?.variants) {
      state.snvSearchVariants = state.snvSearchVariants || {};
      Object.assign(state.snvSearchVariants, payload.variants);
    }
  }
  const genes = new Set(genesUpper);
  const all = {
    ...(state.data?.variants || {}),
    ...(state.snvSearchVariants || {}),
  };
  const matches = Object.entries(all)
    .filter(([, v]) => genes.has((v.gene_symbol || "").toUpperCase()))
    .filter(([, v]) => !filterGnomad || _passesGnomadAfFilter(v));
  matches.sort((a, b) => {
    const sa = Number(a[1].total_score), sb = Number(b[1].total_score);
    return (Number.isFinite(sb) ? sb : -Infinity) - (Number.isFinite(sa) ? sa : -Infinity);
  });
  return matches;
}

function _geneSearchCnvSv(genesUpper) {
  const genes = new Set(genesUpper);
  const all = [
    ...Object.entries(state.data?.cnv_variants || {}),
    ...Object.entries(state.data?.sv_variants  || {}),
  ];
  const matches = all.filter(([, v]) =>
    (v.gene_list || []).some(g => genes.has(String(g).toUpperCase()))
  );
  matches.sort((a, b) => {
    const ra = Number(a[1].cnv_sv_sort_score), rb = Number(b[1].cnv_sv_sort_score);
    return (Number.isFinite(rb) ? rb : -Infinity) - (Number.isFinite(ra) ? ra : -Infinity);
  });
  return matches;
}

let _geneSearchToken = 0;
async function renderGeneSearchResults(kind, rawGenes) {
  const token = ++_geneSearchToken;
  const titleEl = document.getElementById("gene-search-title");
  const host    = document.getElementById("gene-search-results");
  if (!titleEl || !host) return;
  host.innerHTML = "";
  const genesUpper = _parseGeneSearch(rawGenes);
  const geneLabel = genesUpper.join("、");
  if (!genesUpper.length) {
    titleEl.textContent = "基因變異搜尋";
    host.innerHTML = `<div class="muted" style="padding:12px">輸入基因名稱以搜尋。</div>`;
    return;
  }
  const filterGnomad = document.getElementById("gene-search-filter-gnomad-af")?.checked ?? true;
  if (kind === "all") {
    host.innerHTML = `<div class="muted" style="padding:12px">搜尋中…</div>`;
    let snvMatches;
    try {
      snvMatches = await _geneSearchSnv(genesUpper, { filterGnomad });
    } catch (e) {
      if (token === _geneSearchToken) {
        host.innerHTML = `<div class="muted" style="padding:12px">搜尋失敗：${escapeHtml(e.message || String(e))}</div>`;
      }
      return;
    }
    if (token !== _geneSearchToken) return;
    host.innerHTML = "";
    const cnvMatches = _geneSearchCnvSv(genesUpper);
    titleEl.textContent = `${geneLabel} 的所有變異（SNV/Indel: ${snvMatches.length}，CNV/SV: ${cnvMatches.length}）`;
    if (!snvMatches.length && !cnvMatches.length) {
      host.innerHTML = `<div class="muted" style="padding:12px">找不到 ${escapeHtml(geneLabel)} 的變異。</div>`;
      return;
    }
    if (snvMatches.length) {
      const h = document.createElement("h3");
      h.className = "gene-search-section";
      h.textContent = `SNV / Indel（${snvMatches.length}）`;
      host.appendChild(h);
      snvMatches.forEach(([id, v], i) =>
        host.appendChild(renderVariantCard(v, id, "candidate", { index: i + 1, diseaseCheckbox: true })));
    }
    if (cnvMatches.length) {
      const h = document.createElement("h3");
      h.className = "gene-search-section";
      h.textContent = `CNV / SV（${cnvMatches.length}）`;
      host.appendChild(h);
      cnvMatches.forEach(([id, v], i) =>
        host.appendChild(renderCnvSvCard(v, id, { index: i + 1 })));
    }
    return;
  }
  const label = kind === "snv" ? "SNV/Indel" : "CNV/SV";
  if (kind === "snv") {
    host.innerHTML = `<div class="muted" style="padding:12px">搜尋中…</div>`;
    let matches;
    try {
      matches = await _geneSearchSnv(genesUpper, { filterGnomad });
    } catch (e) {
      if (token === _geneSearchToken) {
        host.innerHTML = `<div class="muted" style="padding:12px">搜尋失敗：${escapeHtml(e.message || String(e))}</div>`;
      }
      return;
    }
    if (token !== _geneSearchToken) return;
    host.innerHTML = "";
    titleEl.textContent = `${geneLabel} 的 ${label} 變異（${matches.length}）`;
    if (!matches.length) {
      host.innerHTML = `<div class="muted" style="padding:12px">找不到 ${escapeHtml(geneLabel)} 的變異。</div>`;
      return;
    }
    matches.forEach(([id, v], i) => {
      host.appendChild(renderVariantCard(v, id, "candidate", { index: i + 1, diseaseCheckbox: true }));
    });
  } else {
    const matches = _geneSearchCnvSv(genesUpper);
    titleEl.textContent = `${geneLabel} 的 ${label} 變異（${matches.length}）`;
    if (!matches.length) {
      host.innerHTML = `<div class="muted" style="padding:12px">找不到涵蓋 ${escapeHtml(geneLabel)} 的 CNV/SV。</div>`;
      return;
    }
    matches.forEach(([id, v], i) => {
      host.appendChild(renderCnvSvCard(v, id, { index: i + 1 }));
    });
  }
}

// Combined SNV/Indel + CNV/SV search for one gene — the "搜尋同基因"
// button on each variant card. Used to spot compound-het / mixed-mode
// hits while reviewing an AR candidate.
function openSameGeneModal(gene) {
  openGeneSearchModal("all", gene);
}

function openGeneSearchModal(kind, gene) {
  const inp = document.getElementById("gene-search-modal-input");
  document.getElementById("gene-search-filter-row")
    ?.classList.toggle("hidden", kind === "cnv-sv");
  document.getElementById("gene-search-filter-gnomad-label")
    ?.classList.remove("hidden");
  if (inp) {
    inp.style.display = "";
    inp.value = gene || "";
    inp.dataset.kind = kind;
  }
  renderGeneSearchResults(kind, gene || "");
  showModal("gene-search-modal");
  inp?.focus();
}

// LIRICAL / Exomiser top-20 list, sharing the gene-search modal
// shell. Variants with a per-variant rank 1–20 for the chosen tool,
// sorted by rank ascending; the card's #N marker shows the rank.
function openToolRankModal(tool) {
  const titleEl = document.getElementById("gene-search-title");
  const host    = document.getElementById("gene-search-results");
  const inp      = document.getElementById("gene-search-modal-input");
  document.getElementById("gene-search-filter-row")?.classList.remove("hidden");
  document.getElementById("gene-search-filter-gnomad-label")?.classList.add("hidden");
  if (inp) inp.style.display = "none";   // re-search input is meaningless in rank mode
  const rankKey  = tool === "lirical" ? "rank_lirical_variant" : "rank_exomiser_variant";
  const toolName = tool === "lirical" ? "LIRICAL" : "Exomiser";
  const matches = Object.entries(state.data?.variants || {})
    .map(([id, v]) => [id, v, Number(v[rankKey])])
    .filter(([, , r]) => Number.isFinite(r) && r >= 1 && r <= 20)
    .sort((a, b) => a[2] - b[2]);
  if (titleEl) titleEl.textContent = `${toolName} rank 1–20（${matches.length}）`;
  if (host) {
    host.innerHTML = "";
    if (!matches.length) {
      host.innerHTML = `<div class="muted" style="padding:12px">沒有 ${toolName} rank 1–20 的變異（可能還沒跑分析）。</div>`;
    } else {
      matches.forEach(([id, v, rank]) => {
        host.appendChild(renderVariantCard(v, id, "candidate", { index: rank, diseaseCheckbox: true }));
      });
    }
  }
  showModal("gene-search-modal");
}

function setupGeneSearch() {
  // Delegated: per-card "搜尋同基因" buttons (SNV + CNV/SV cards live
  // both in the main view and in the modal, plus they're re-rendered
  // a lot, so capture-phase delegation is simpler than per-render
  // binding).
  document.addEventListener("click", ev => {
    const btn = ev.target.closest(".same-gene-btn");
    if (!btn) return;
    const g = btn.getAttribute("data-gene");
    if (g) openSameGeneModal(g);
  });
  document.querySelectorAll(".gene-search-input").forEach(inp => {
    inp.addEventListener("keydown", ev => {
      if (ev.key !== "Enter") return;
      ev.preventDefault();
      const g = inp.value.trim();
      if (!g) return;
      openGeneSearchModal(inp.dataset.kind || "snv", g);
    });
  });
  document.querySelectorAll(".ac-tool-btn").forEach(btn => {
    btn.addEventListener("click", () => openToolRankModal(btn.dataset.tool || "exomiser"));
  });
  document.getElementById("gene-search-modal-input")?.addEventListener("keydown", ev => {
    if (ev.key !== "Enter") return;
    ev.preventDefault();
    const inp = ev.currentTarget;
    renderGeneSearchResults(inp.dataset.kind || "snv", inp.value);
  });
  document.getElementById("gene-search-filter-gnomad-af")?.addEventListener("change", () => {
    const inp = document.getElementById("gene-search-modal-input");
    if (inp && inp.style.display !== "none") {
      renderGeneSearchResults(inp.dataset.kind || "snv", inp.value);
    }
  });
}

// 上傳個案清單: opens a modal listing past uploads + a button that
// picks an xlsx → POST /api/patient_list → re-renders history. The
// roster it builds is what the 載入新個案 modal reads to auto-fill
// MRN / 姓名 / Test type.
function _fmtUploadTime(iso) {
  if (!iso) return "—";
  try { return new Date(iso).toLocaleString(); } catch { return iso; }
}

async function _renderPatientListHistory() {
  const host = document.getElementById("patient-list-history");
  if (!host) return;
  host.innerHTML = `<div class="muted">載入中…</div>`;
  try {
    const rows = await apiFetch("/patient_list/uploads") || [];
    if (!rows.length) {
      host.innerHTML = `<div class="muted" style="padding:10px">（尚無上傳記錄）</div>`;
      return;
    }
    const head = `
      <tr>
        <th>時間</th>
        <th>檔名</th>
        <th class="num">解析</th>
        <th class="num">新增</th>
        <th class="num">更新</th>
        <th class="num">roster</th>
        <th>封存檔</th>
      </tr>`;
    const body = rows.map(r => `
      <tr>
        <td>${escapeHtml(_fmtUploadTime(r.uploaded_at))}</td>
        <td class="fn" title="${escapeAttr(r.original_filename || "")}">${escapeHtml(r.original_filename || "—")}</td>
        <td class="num">${r.parsed ?? "—"}</td>
        <td class="num">${r.added ?? "—"}</td>
        <td class="num">${r.updated ?? "—"}</td>
        <td class="num">${r.total_after ?? "—"}</td>
        <td class="fn muted" title="${escapeAttr(r.archive_name || "")}">${escapeHtml(r.archive_name || "—")}</td>
      </tr>`).join("");
    host.innerHTML = `<table>${head}${body}</table>`;
  } catch (e) {
    host.innerHTML = `<div class="muted" style="padding:10px">載入失敗：${escapeHtml(String(e))}</div>`;
  }
}

function setupPatientListUpload() {
  const btn  = document.getElementById("btn-upload-list");
  const file = document.getElementById("upload-list-file");
  const pick = document.getElementById("patient-list-pick-btn");
  const status = document.getElementById("patient-list-status");
  if (!btn || !file) return;
  btn.addEventListener("click", async () => {
    showModal("patient-list-modal");
    if (status) status.textContent = "";
    await _renderPatientListHistory();
  });
  pick?.addEventListener("click", () => { file.value = ""; file.click(); });
  file.addEventListener("change", async () => {
    const f = file.files && file.files[0];
    if (!f) return;
    if (pick) { pick.disabled = true; pick.textContent = "上傳中…"; }
    if (status) status.textContent = "";
    try {
      const fd = new FormData();
      fd.append("file", f, f.name);
      const resp = await fetch(`${API_BASE}/patient_list`, {
        method: "POST",
        credentials: "same-origin",
        body: fd,
      });
      if (resp.status === 401) { showLoginModal(); throw new Error("尚未登入"); }
      const body = await resp.json().catch(() => ({}));
      if (!resp.ok) throw new Error(body.detail || `${resp.status} ${resp.statusText}`);
      if (status) status.textContent = `✓ ${f.name}：解析 ${body.parsed} · 新增 ${body.added} · 更新 ${body.updated} · roster ${body.total}`;
      _unregisteredById = {};
      _rosterOptions = null;          // refresh 科別/開單醫師 datalists
      await _renderPatientListHistory();
    } catch (e) {
      if (status) status.textContent = "上傳失敗：" + (e.message || e);
    } finally {
      if (pick) { pick.disabled = false; pick.textContent = "選擇 xlsx 上傳"; }
    }
  });
}

// 個案清單: lists registered NGS-UI samples and lets reviewers remove
// registration/report state while keeping local tertiary outputs available.
let _caseListRows = [];
const _caseListTestFilters = new Set(SAMPLE_TEST_TYPES);

async function loadCaseListRows() {
  const rows = await apiFetch("/samples/case-summary") || [];
  return rows.map(row => ({
    ...row,
    test_type: normalizeSampleTestType(
      row.test_type || "",
      row.lis_id || row.sample_id || "",
    ),
  }));
}

function _caseListVisibleRows() {
  const query = (document.getElementById("case-list-search")?.value || "").trim().toLowerCase();
  return _caseListRows.filter(row => {
    const testType = (row.test_type || "").toUpperCase();
    if (!_caseListTestFilters.has(testType)) return false;
    if (!query) return true;
    return Object.values(row).some(value => String(value ?? "").toLowerCase().includes(query));
  });
}

async function _renderCaseList({ refresh = true } = {}) {
  const host = document.getElementById("case-list-table");
  const status = document.getElementById("case-list-status");
  if (!host) return;
  if (refresh) host.innerHTML = `<div class="muted" style="padding:10px">載入中…</div>`;
  try {
    if (refresh) _caseListRows = await loadCaseListRows();
    const rows = _caseListVisibleRows();
    if (!_caseListRows.length) {
      host.innerHTML = `<div class="muted" style="padding:10px">（尚無已載入個案）</div>`;
      return;
    }
    if (!rows.length) {
      host.innerHTML = `<div class="muted" style="padding:10px">（沒有符合篩選條件的個案）</div>`;
      return;
    }
    const head = `
      <tr>
        <th>LIS_ID</th>
        <th>姓名</th>
        <th>病歷號</th>
        <th>Test type</th>
        <th>Phenotype</th>
        <th>Causative variant</th>
        <th>Disease</th>
        <th>Other variant</th>
        <th>Comment</th>
        <th>簽收時間</th>
        <th>載入時間</th>
        <th></th>
      </tr>`;
    const body = rows.map(r => `
      <tr>
        <td class="case-list-id">${escapeHtml(r.lis_id || r.sample_id || "")}</td>
        <td class="case-list-id">${escapeHtml(r.name || "—")}</td>
        <td class="case-list-id">${escapeHtml(r.mrn || "—")}</td>
        <td class="case-list-test">${escapeHtml(r.test_type || "—")}</td>
        <td class="case-list-long case-list-phenotype">${escapeHtml(r.phenotype_summary || "—")}</td>
        <td class="case-list-long case-list-variant">${escapeHtml(r.causative_variants || "—")}</td>
        <td class="case-list-long case-list-disease">${escapeHtml(r.diseases || "—")}</td>
        <td class="case-list-long case-list-variant">${escapeHtml(r.other_variants || "—")}</td>
        <td class="case-list-long">${escapeHtml(r.comment || "—")}</td>
        <td class="case-list-date">${escapeHtml(_fmtUploadTime(r.sign_received_at))}</td>
        <td class="case-list-date">${escapeHtml(_fmtUploadTime(r.created_at))}</td>
        <td><button type="button" class="btn btn-danger case-list-delete"
          data-sample-id="${escapeAttr(r.lis_id || r.sample_id || "")}">刪除</button></td>
      </tr>`).join("");
    host.innerHTML = `<table>${head}${body}</table>`;
    if (status) status.textContent = "";
  } catch (e) {
    host.innerHTML = `<div class="muted" style="padding:10px">載入失敗：${escapeHtml(String(e))}</div>`;
  }
}

function setupCaseList() {
  const btn = document.getElementById("btn-case-list");
  const host = document.getElementById("case-list-table");
  const status = document.getElementById("case-list-status");
  const search = document.getElementById("case-list-search");
  let pendingSid = "";
  let pendingDeleteButton = null;
  if (!btn || !host) return;
  btn.addEventListener("click", async () => {
    showModal("case-list-modal");
    if (status) status.textContent = "";
    if (!await flushPendingSave()) {
      if (status) status.textContent = `儲存失敗：${_saveError || "尚未完成儲存"}`;
      return;
    }
    await _renderCaseList();
  });
  search?.addEventListener("input", () => _renderCaseList({ refresh: false }));
  document.querySelectorAll(".case-list-test-filters input").forEach(filter => {
    filter.addEventListener("change", () => {
      if (filter.checked) _caseListTestFilters.add(filter.value);
      else _caseListTestFilters.delete(filter.value);
      _renderCaseList({ refresh: false });
    });
  });
  document.querySelectorAll(".case-list-test-filters .sample-test-only").forEach(button => {
    button.addEventListener("click", () => {
      const selected = button.dataset.testType || "";
      _caseListTestFilters.clear();
      if (selected) _caseListTestFilters.add(selected);
      document.querySelectorAll(".case-list-test-filters input").forEach(filter => {
        filter.checked = filter.value === selected;
      });
      _renderCaseList({ refresh: false });
    });
  });
  host.addEventListener("click", async ev => {
    const del = ev.target.closest?.(".case-list-delete");
    if (!del) return;
    const sid = del.dataset.sampleId || "";
    if (!sid) return;
    pendingSid = sid;
    pendingDeleteButton = del;
    document.getElementById("case-delete-sid").textContent = sid;
    document.getElementById("case-delete-ui-path").textContent = `NGS_UI/tertiary_output/${sid}/sample_metadata.json`;
    showModal("case-delete-modal");
  });

  const cancelDelete = () => {
    pendingSid = "";
    pendingDeleteButton = null;
    hideModal("case-delete-modal");
  };

  async function deletePendingCase() {
    const sid = pendingSid;
    const del = pendingDeleteButton;
    if (!sid) return;
    hideModal("case-delete-modal");
    pendingSid = "";
    pendingDeleteButton = null;
    if (del) del.disabled = true;
    if (status) status.textContent = `刪除 ${sid} 中…`;
    try {
      const resp = await fetch(
        `${API_BASE}/samples/${encodeURIComponent(sid)}?delete_pipeline_output=false`,
        { method: "DELETE", credentials: "same-origin" },
      );
      if (resp.status === 401) { showLoginModal(); throw new Error("尚未登入"); }
      const body = await resp.json().catch(() => ({}));
      if (!resp.ok) throw new Error(body.detail || `${resp.status} ${resp.statusText}`);
      if (state.currentLIS === sid) {
        window.location.reload();
        return;
      }
      await loadIndex();
      _unregisteredCache.loadedAt = 0;
      _unregisteredCache.list = null;
      _caseListRows = await loadCaseListRows();
      await _renderCaseList({ refresh: false });
      if (status) {
        status.textContent = `已刪除 ${sid} 的載入狀態，三級分析輸出檔案已保留，可重新載入。`;
      }
    } catch (e) {
      if (status) status.textContent = `刪除失敗：${e.message || e}`;
      if (del) del.disabled = false;
    }
  }

  document.getElementById("case-delete-cancel")?.addEventListener("click", cancelDelete);
  document.getElementById("case-delete-confirm")?.addEventListener("click", () => {
    deletePendingCase();
  });
}

function setLoggedInUser(username) {
  const span = document.getElementById("topbar-user");
  const btn  = document.getElementById("btn-logout");   // doubles as the 登入 button when signed out
  if (span) { span.textContent = username; span.hidden = !username; }
  if (btn) {
    btn.hidden = false;                                  // always visible — toggles label/action
    btn.textContent = username ? "登出" : "登入";
    btn.dataset.loggedIn = username ? "1" : "0";
  }
  const up = document.getElementById("btn-upload-list");
  if (up) up.hidden = !username;
  const cases = document.getElementById("btn-case-list");
  if (cases) cases.hidden = !username;
  const secondary = document.getElementById("btn-secondary-launch");
  if (secondary) secondary.hidden = !username;
  const dr = document.getElementById("btn-dragen-launch");
  if (dr) dr.hidden = !username;
  const pipelineList = document.getElementById("btn-pipeline-list");
  if (pipelineList) pipelineList.hidden = !username;
  // #btn-phenotype-tool is intentionally always visible — the HPO/panel
  // tool needs no login (it runs on the intranet), so the link stays
  // reachable even before sign-in.
}

async function handleLogin(ev) {
  ev?.preventDefault();
  const u = document.getElementById("login-username").value.trim();
  const p = document.getElementById("login-password").value;
  try {
    const me = await fetch(`${API_BASE}/auth/login`, {
      method: "POST",
      credentials: "same-origin",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ username: u, password: p }),
    }).then(async r => {
      if (r.status === 401) throw new Error("帳號或密碼錯誤");
      if (!r.ok) throw new Error(`登入失敗 (${r.status})`);
      return r.json();
    });
    document.getElementById("login-password").value = "";
    hideLoginModal();
    setLoggedInUser(me.username);
    await bootAfterAuth();
  } catch (e) {
    showLoginModal(e.message);
  }
}

// Download diagnostic DOCX. First asks the reviewer how 「本次檢測
// 基因包括」(§五.4) should be rendered — per-panel/HPO sections, or
// a single deduped flat list. The backend saves an archived copy under
// NGS_UI/report/ and also streams it back as a download.
async function exportDiagnosticDocx() {
  if (!state.currentLIS) return;
  const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
  const sid = row?.sample_id || state.currentLIS;

  const mode = await _pickGeneListMode();
  if (!mode) return;   // cancelled

  try {
    if (!await flushPendingSave()) throw new Error(_saveError || "尚未完成儲存");
    const url = `${API_BASE}/samples/${encodeURIComponent(sid)}/report.docx?gene_list_mode=${mode}`;
    const resp = await fetch(url, { credentials: "same-origin" });
    if (resp.status === 401) { showLoginModal(); return; }
    if (!resp.ok) throw new Error(`匯出失敗 (${resp.status})`);
    const blob = await resp.blob();
    downloadBlob(blob, `${sid}_diagnosis.docx`);
  } catch (e) {
    alert("匯出失敗：" + e.message);
  }
}

function _pickHealthReportSections() {
  const options = [
    { key: "acmg_sf", title: "ACMG 疾病風險基因", checked: true },
    { key: "stroke", title: "中風相關基因", checked: false },
    { key: "carrier", title: "帶因者篩查", checked: false },
    { key: "pgx", title: "藥物基因體學", checked: true },
  ];
  return new Promise((resolve) => {
    const wrap = document.createElement("div");
    wrap.className = "modal";
    wrap.innerHTML = `
      <div class="modal-card health-report-card">
        <h2>匯出健檢報告</h2>
        <p class="health-report-hint">選擇這次要放入報告的項目。</p>
        <div class="health-report-options">
          ${options.map(opt => `
            <label class="health-report-option">
              <input type="checkbox" value="${escapeAttr(opt.key)}"${opt.checked ? " checked" : ""} />
              <span>${escapeHtml(opt.title)}</span>
            </label>
          `).join("")}
        </div>
        <div class="modal-actions">
          <button type="button" class="btn btn-ghost" data-act="cancel">取消</button>
          <button type="button" class="btn btn-primary" data-act="ok">匯出</button>
        </div>
      </div>`;
    document.body.appendChild(wrap);
    const finish = (val) => { document.body.removeChild(wrap); resolve(val); };
    wrap.addEventListener("click", (ev) => {
      if (ev.target === wrap) finish(null);
      const act = ev.target?.dataset?.act;
      if (act === "cancel") finish(null);
      if (act === "ok") {
        const selected = Array.from(wrap.querySelectorAll('input[type="checkbox"]:checked'))
          .map(input => input.value);
        if (!selected.length) {
          alert("請至少選擇一個報告項目。");
          return;
        }
        finish(selected);
      }
    });
  });
}

async function exportHealthDocx() {
  if (!state.currentLIS) return;
  const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
  const sid = row?.sample_id || state.currentLIS;
  const sections = await _pickHealthReportSections();
  if (!sections) return;

  const btns = document.querySelectorAll(".btn-export-health");
  const setBusy = b => btns.forEach(x => { x.disabled = b; });
  const hint = msg => {
    document.querySelectorAll(".js-save-hint").forEach(el => { el.textContent = msg; });
  };
  setBusy(true);
  hint("產生健檢報告…");
  try {
    if (!await flushPendingSave()) throw new Error(_saveError || "尚未完成儲存");
    const qs = encodeURIComponent(sections.join(","));
    const url = `${API_BASE}/samples/${encodeURIComponent(sid)}/health-report.docx?sections=${qs}`;
    const resp = await fetch(url, { credentials: "same-origin" });
    if (resp.status === 401) { showLoginModal(); return; }
    if (!resp.ok) throw new Error(`匯出失敗 (${resp.status})`);
    const blob = await resp.blob();
    downloadBlob(blob, `${sid}_health.docx`);
    hint("健檢報告已下載");
  } catch (e) {
    hint("");
    alert("匯出健檢報告失敗：" + e.message);
  } finally {
    setBusy(false);
  }
}

function _reportPrintBlock(sectionId, { omitWhenEmpty = false } = {}) {
  const source = document.getElementById(sectionId);
  if (!source) return "";
  const clone = source.cloneNode(true);
  const cards = clone.querySelectorAll(".variant-card");
  if (omitWhenEmpty && !cards.length) return "";

  // The PDF is a reviewer-facing card summary, not an editable copy of
  // the UI. Drop controls, comments, external links, and collapsed
  // drill-down content before serialising the print window.
  clone.querySelectorAll([
    ".status-radio", ".ext-links", ".btn-copy", ".btn-more",
    ".more-extras", ".comment-row", ".cnv-sv-comment", ".transcript-picker",
    ".btn-add-manual", ".btn-remove-manual", ".same-gene-btn",
    ".disease-detail", ".disease-collapse",
    ".disease-needs-description-star", ".disease-source-badges",
    ".cnv-sv-reasoning", ".cnv-sv-gene-overflow",
    ".block-header .count",
  ].join(",")).forEach(el => el.remove());
  clone.querySelectorAll(".cnv-sv-section-title").forEach(title => {
    if (["已知致病區域重疊", "已知良性區域重疊"].includes(title.textContent.trim())) {
      title.closest(".cnv-sv-section")?.remove();
    }
  });
  clone.querySelectorAll(".disease-summary-text").forEach(summary => {
    const text = summary.textContent || "";
    const inheritance = text.match(
      /\((?:AD|AR|XLD|XLR|XL|YL|MT|MI|DR|DD|SMU|MU|ISOL)(?:\s*[/,;]\s*(?:AD|AR|XLD|XLR|XL|YL|MT|MI|DR|DD|SMU|MU|ISOL))*\)/i,
    );
    if (inheritance) {
      summary.textContent = text.slice(0, inheritance.index + inheritance[0].length);
    }
  });
  clone.querySelectorAll(".manual-row").forEach(row => {
    if (row.querySelector(".manual-comment")) row.remove();
  });
  clone.querySelectorAll("button, a").forEach(el => el.remove());
  clone.querySelectorAll("input, select, textarea").forEach(control => {
    if (control.matches('[type="checkbox"], [type="radio"]')) {
      control.remove();
      return;
    }
    const text = control.tagName === "SELECT"
      ? control.options[control.selectedIndex]?.textContent || ""
      : control.value || "";
    const span = document.createElement("span");
    span.className = "print-field";
    span.textContent = text;
    control.replaceWith(span);
  });
  clone.querySelectorAll("details").forEach(details => {
    details.removeAttribute("open");
    details.querySelectorAll(":scope > :not(summary)").forEach(el => el.remove());
  });
  return clone.outerHTML;
}

function _printReportTimestamp(now = new Date()) {
  const pad = n => String(n).padStart(2, "0");
  return `${now.getFullYear()}/${now.getMonth() + 1}/${now.getDate()} ${pad(now.getHours())}:${pad(now.getMinutes())}`;
}

function _reportGeneListPrintBlock(mode, data = {}) {
  const sections = Array.isArray(data.grouped) ? data.grouped : [];
  const merged = Array.isArray(data.merged) ? data.merged : [];
  if (!sections.length && !merged.length) {
    return `<section class="print-gene-list"><h2>本次檢測基因包括</h2><p>（未設定 HPO / panel — 無檢測基因清單）</p></section>`;
  }
  if (mode === "merged") {
    return `<section class="print-gene-list"><h2>本次檢測基因包括</h2><p>${escapeHtml(merged.join(", "))}</p></section>`;
  }
  const html = sections.map(section => `
    <div class="print-gene-section">
      <h3>${escapeHtml(section.name || "")}</h3>
      <p>${escapeHtml((section.genes || []).length ? section.genes.join(", ") : "（無對應基因）")}</p>
    </div>`).join("");
  return `<section class="print-gene-list"><h2>本次檢測基因包括</h2>${html}</section>`;
}

async function printReportCards() {
  if (!state.currentLIS || !state.data) return;
  const mode = await _pickGeneListMode({ title: "輸出 PDF", okLabel: "輸出 PDF" });
  if (!mode) return;
  const sampleRow = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
  const sampleId = sampleRow?.sample_id || state.currentLIS;
  let reportGeneList = {};
  try {
    reportGeneList = await apiFetch(`/samples/${encodeURIComponent(sampleId)}/report-gene-list`) || {};
  } catch (e) {
    alert("讀取基因清單失敗：" + e.message);
    return;
  }
  const popup = window.open("", "_blank");
  if (!popup) {
    alert("無法開啟列印視窗，請允許本站開啟彈出式視窗後再試一次。");
    return;
  }

  // Refresh the report sections so the print copy includes the newest
  // status choices even when auto-save has not fired yet.
  renderReportSections();
  const meta = state.data.meta || {};
  const sid = meta.LIS_ID || state.currentLIS;
  const printHeader = `${sid}_report_${_printReportTimestamp()}`;
  const printHeaderCss = JSON.stringify(printHeader);
  const stylesheet = document.querySelector('link[rel="stylesheet"]')?.href || "";
  const sections = [
    _reportPrintBlock("sec-causative"),
    _reportPrintBlock("sec-other"),
    _reportPrintBlock("sec-candidate", { omitWhenEmpty: true }),
  ].join("");
  const geneList = _reportGeneListPrintBlock(mode, reportGeneList);

  popup.document.open();
  popup.document.write(`<!doctype html>
<html lang="zh-Hant">
<head>
  <meta charset="utf-8" />
  <title>${escapeHtml(sid)}_report_${todayYmd()}</title>
  ${stylesheet ? `<link rel="stylesheet" href="${escapeAttr(stylesheet)}" />` : ""}
  <style>
    body { margin: 0; background: #fff; color: #24292f; font-family: Arial, "Noto Sans TC", sans-serif; }
    .print-page { max-width: 1180px; margin: 0 auto; padding: 16px 14px; }
    .print-toolbar { display: flex; justify-content: flex-end; margin-bottom: 12px; }
    .print-title { margin: 0 0 14px; font-size: 24px; }
    .report-block { margin-bottom: 12px; }
    .block-header, .block-body { display: block !important; }
    .block-header { pointer-events: none; }
    .block-header .arrow { display: none; }
    .variant-card { break-inside: avoid; page-break-inside: avoid; }
    .print-gene-list { break-inside: avoid; page-break-inside: avoid; margin-top: 18px; }
    .print-gene-list h2 { font-size: 16px; margin: 0 0 8px; }
    .print-gene-section { margin: 0 0 10px; }
    .print-gene-section h3 { font-size: 13px; margin: 0 0 3px; }
    .print-gene-list p { margin: 0; font-size: 11px; line-height: 1.45; overflow-wrap: anywhere; }
    .print-field { display: inline-block; white-space: pre-wrap; overflow-wrap: anywhere; }
    .disease-row summary { line-height: 1.15; padding-top: 1px; padding-bottom: 1px; }
    details > summary { list-style: none; }
    details > summary::before { display: none !important; }
    @media print {
      @page {
        size: A4;
        margin: 12mm 5mm 10mm;
        @top-right { content: ${printHeaderCss}; color: #57606a; font-size: 9px; }
        @bottom-right { content: counter(page); color: #57606a; font-size: 9px; }
      }
      .print-page { max-width: none; padding: 3mm 0 0; }
      .print-toolbar { display: none; }
      .variant-card { box-shadow: none; }
    }
  </style>
</head>
<body>
  <main class="print-page">
    <div class="print-toolbar"><button type="button" onclick="window.print()">列印 / 儲存 PDF</button></div>
    <h1 class="print-title">${escapeHtml(sid)} Report</h1>
    ${sections}
    ${geneList}
  </main>
  <script>window.addEventListener("load", () => setTimeout(() => window.print(), 250));<\/script>
</body>
</html>`);
  popup.document.close();
}

// Small inline modal asking how the 檢測基因清單 should appear in §五.4.
// Returns "grouped" | "merged" | null (cancelled).
function _pickGeneListMode(opts = {}) {
  return new Promise((resolve) => {
    // One-shot modal — built on the fly, removed after pick.
    const wrap = document.createElement("div");
    wrap.className = "modal";
    wrap.innerHTML = `
      <div class="modal-card" style="max-width:520px">
        <h2>${escapeHtml(opts.title || "匯出診斷報告")}</h2>
        <p style="margin:8px 0 12px;font-size:13px;color:#555">
          「本次檢測基因包括」要怎麼呈現？
        </p>
        <label class="gene-list-mode-opt" style="display:block;margin:8px 0;padding:8px;border:1px solid var(--border);border-radius:6px;cursor:pointer">
          <input type="radio" name="gene-list-mode" value="grouped" checked />
          <strong>按 panel / HPO 分組</strong>
          <div style="font-size:12px;color:#777;margin-left:22px">每一個 HPO term 或 panel 各自一段，標題後接該組基因清單</div>
        </label>
        <label class="gene-list-mode-opt" style="display:block;margin:8px 0;padding:8px;border:1px solid var(--border);border-radius:6px;cursor:pointer">
          <input type="radio" name="gene-list-mode" value="merged" />
          <strong>全部合併（去重）</strong>
          <div style="font-size:12px;color:#777;margin-left:22px">所有 HPO + panel 的基因取聯集，去重後列成一行</div>
        </label>
        <div class="dragen-actions" style="margin-top:14px;display:flex;gap:8px;justify-content:flex-end">
          <button type="button" class="btn btn-ghost" data-act="cancel">取消</button>
          <button type="button" class="btn btn-primary" data-act="ok">${escapeHtml(opts.okLabel || "匯出")}</button>
        </div>
      </div>`;
    document.body.appendChild(wrap);
    const finish = (val) => { document.body.removeChild(wrap); resolve(val); };
    wrap.addEventListener("click", (ev) => {
      if (ev.target === wrap) finish(null);
      const act = ev.target?.dataset?.act;
      if (act === "cancel") finish(null);
      if (act === "ok") {
        const sel = wrap.querySelector('input[name="gene-list-mode"]:checked');
        finish(sel?.value || "grouped");
      }
    });
  });
}

async function handleLogout() {
  await apiPost("/auth/logout", {}).catch(() => {});
  setLoggedInUser("");
  // Easiest reset: full reload returns to a clean state with the modal up.
  location.reload();
}

// ---- Phase C: Exomiser/LIRICAL rerun --------------------------------

let _jobPollTimer = null;
let _activeJobId  = null;

// Map a worker step like "exomiser:run" / "lirical:render" to the
// short tool name shown in the status pill. Anything else (parse,
// queued, done, …) becomes empty so the pill just shows the status.
function _stepTool(step) {
  const s = String(step || "");
  if (s.startsWith("exomiser")) return "exomiser";
  if (s.startsWith("lirical"))  return "lirical";
  return "";
}

function _setJobStatus(text, busy = false) {
  const el = document.getElementById("job-status");
  if (el) el.textContent = text || "";
  const btn = document.getElementById("btn-rerun-tools");
  if (btn) btn.disabled = !!busy;
}

// 「開始分析」: instant in-house pheno_score; Exomiser/LIRICAL is queued
// only when the analysis contains at least one HPO term.
// The phenotype POST returns immediately so cards can apply the fresh
// pheno_score.tsv-derived in-panel state; Exomiser/LIRICAL runs in the
// background and the polling loop refreshes the sample once it lands.
//
// `opts.version` selects which analysis version to write into; `opts.mode`
// is "overwrite" (clear sidecars first) or "new" (create fresh version).
// Defaults: overwrite the currently-active version, no clear (legacy).
async function startAnalysis(opts = {}) {
  if (!state.currentLIS) return;
  const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
  const sid = row?.sample_id || state.currentLIS;
  const hint = document.getElementById("phenotype-hint");
  const version = opts.version || state.data?.active_analysis || "default";
  const mode    = opts.mode    || "overwrite";
  const hasHpo  = Array.isArray(phenoEdit.hpo) && phenoEdit.hpo.length > 0;

  _setJobStatus("送出工作中…", true);
  hint.textContent = "計算中…";

  // Make sure the chosen version exists and is the active one before we
  // write to it. Overwriting an existing version wipes its sidecars so
  // a partial re-run can't blend old + new outputs.
  try {
    await apiPost(`/samples/${encodeURIComponent(sid)}/analyses`, {
      name:           version,
      hpo:            phenoEdit.hpo,
      panels:         phenoEdit.panels,
      set_active:     true,
      clear_sidecars: mode === "overwrite",
    });
  } catch (e) {
    hint.textContent = "失敗：" + e.message;
    _setJobStatus("失敗", false);
    return;
  }

  try {
    const result = await apiPost(`/samples/${encodeURIComponent(sid)}/phenotype`, {
      hpo:    phenoEdit.hpo,
      panels: phenoEdit.panels,
      version,
    });
    document.getElementById("phenotype-stats").textContent =
      `${result.n_hpo} HPO + ${result.n_panels} panels → ${result.n_in_panel_genes} genes in panel · top ${(result.top_score ?? 0).toFixed(0)}`;
    const top10El = document.getElementById("phenotype-top10");
    const top10Ul = document.getElementById("phenotype-top10-list");
    top10Ul.innerHTML = (result.top10 || []).map(x =>
      `<li><span class="mane-tx">${escapeHtml(x.gene)}</span> &nbsp; ${x.score.toFixed(2)}</li>`
    ).join("");
    top10El.classList.toggle("hidden", !(result.top10 && result.top10.length));
    hint.textContent = `已重算 (${new Date().toLocaleTimeString()})`;
    // Refresh so cards see the freshly written pheno_score.tsv right
    // away, before the slower Exomiser job finishes.
    await loadSample(state.currentLIS);
    renderAll();
  } catch (e) {
    hint.textContent = "失敗：" + e.message;
    _setJobStatus("失敗", false);
    return;
  }
  if (!hasHpo) {
    _setJobStatus("", false);
    return;
  }
  try {
    const job = await apiPost(`/samples/${encodeURIComponent(sid)}/jobs/exomiser_lirical`, {
      version,
    });
    _activeJobId = job.job_id;
    _setJobStatus(`已排入：${job.job_id}`, true);
    _startJobPolling(sid, job.job_id);
  } catch (e) {
    _setJobStatus("Exomiser 排入失敗：" + e.message, false);
  }
}

function _startJobPolling(sid, jobId) {
  clearInterval(_jobPollTimer);
  _jobPollTimer = setInterval(async () => {
    try {
      const j = await apiFetch(`/jobs/${encodeURIComponent(jobId)}`);
      if (!j) return;
      const status = j.status || j.rq_status || "?";
      const tool   = _stepTool(j.step);
      _setJobStatus(tool ? `${status} ${tool}` : status,
                    status === "queued" || status === "running");
      if (status === "succeeded" || status === "failed") {
        clearInterval(_jobPollTimer);
        if (status === "succeeded") {
          _setJobStatus(`完成 · Exomiser ${j.n_exomiser_variants ?? 0} / LIRICAL ${j.n_lirical_variants ?? 0} variants`, false);
          // Reload sample so cards pick up the new score columns,
          // without blocking the reviewer behind the global loading modal.
          await loadSample(state.currentLIS, { showLoading: false });
          renderAll();
        } else {
          _setJobStatus(`失敗 (${j.step || ""}) — 看 analysis_files/rerun.log`, false);
        }
      }
    } catch (e) {
      // Network blip — keep polling.
    }
  }, 5000);
}

// ---------- Render: variant card -----------------------------------

function statusOptions(kind) {
  // kind: "candidate" → 1 / 2 / C / 0
  //                     1 = Causative · 2 = Other · C = Candidate · 0 = reviewed
  //       "panel"     → secondary finding panels → selected ✓ only
  // X was dropped — reviewers asked to mark-only, not hide.
  if (kind === "panel") return ["✓"];
  return ["1", "2", "C", "0"];
}

// Render the four status options as radio chips (clickable circles)
// instead of a dropdown — saves a click per status change and the
// available choices are visible at a glance. `panelAttr` carries the
// optional data-panel="..." string for panel-scoped statuses.
let _statusRadioSeq = 0;
function _renderStatusRadio(id, curStatus, opts, panelAttr = "") {
  // A variant can be rendered in both the analysis area and report
  // sections. Candidate widgets use checkboxes because C and 0 may
  // coexist; secondary-finding widgets are also checkboxes because
  // reviewers must be able to unselect a default ClinVar P/LP pick.
  const groupName = `status-${id}-${++_statusRadioSeq}${panelAttr ? "-" + panelAttr.replace(/[^A-Za-z0-9]/g, "") : ""}`;
  const inputType = "checkbox";
  return `<span class="status-radio" data-id="${escapeAttr(id)}" ${panelAttr}>` +
    opts.map(o => {
      const checked = _statusHas(curStatus, o) ? " checked" : "";
      const clsKey = o === "✓" ? "check" : o.toLowerCase();
      const cls = `status-radio-chip status-radio-${clsKey}`;
      return `<label class="${cls}"><input type="${inputType}" name="${escapeAttr(groupName)}" value="${escapeAttr(o)}"${checked} /><span>${escapeHtml(o)}</span></label>`;
    }).join("") +
  `</span>`;
}

function _statusValues(raw) {
  if (Array.isArray(raw)) return raw.map(String);
  const text = String(raw || "").trim();
  if (!text) return [];
  return text.split(",").map(s => s.trim()).filter(Boolean);
}

function _statusHas(raw, option) {
  return _statusValues(raw).includes(option);
}

function getStatus(id) {
  return (state.reports.status && state.reports.status[id]) || "";
}

function _syncStatusRadios(id, panel, val) {
  document.querySelectorAll(".status-radio").forEach(wrap => {
    if (wrap.dataset.id !== id) return;
    if (panel ? !wrap.dataset.panel : !!wrap.dataset.panel) return;
    wrap.querySelectorAll('input[type="radio"], input[type="checkbox"]').forEach(input => {
      input.checked = _statusHas(val, input.value);
    });
  });
}

function setStatus(id, val) {
  state.reports.status = state.reports.status || {};
  if (val) state.reports.status[id] = val;
  else     delete state.reports.status[id];
  state.dirty = true;
  // Scoped re-render: status only affects which Causative / Other /
  // Candidate section the variant lands in. Skip the full renderAll
  // (which re-builds sample meta, phenotype, every tier panel, etc.)
  // so flipping a status pill feels instant instead of taking ~1 s.
  renderReportSections();
  renderDiseaseAssociatedReportWarning();
  _syncStatusRadios(id, "", val);
  updateSaveHint();
}

function toggleStatus(id, option, checked) {
  let next = "";
  if (option === "1" || option === "2") {
    next = checked ? option : "";
  } else {
    const values = new Set(_statusValues(getStatus(id)).filter(v => v === "C" || v === "0"));
    if (checked) values.add(option);
    else values.delete(option);
    next = ["C", "0"].filter(v => values.has(v)).join(",");
  }
  setStatus(id, next);
}

function _isPlpClass(text) {
  return ["sig-p", "sig-lp"].includes(classifySignificance(text));
}

function _isClinvarPlp(v) {
  return _isPlpClass(v?.CLNSIG || "");
}

function _isSecondaryEligible(id) {
  const v = state.data?.variants?.[id];
  // Secondary analysis mirrors the main SNV retrieval buckets: retain all
  // existing ClinVar P/LP calls plus 1A/1B/1C (LOFTEE HC, ACMG points >=4,
  // P-KNN LLR >=1, and the other predictor triggers). isSecondarySelected()
  // below still defaults to ClinVar P/LP only, so the broader candidates are
  // not silently added to the report.
  return _isClinvarPlp(v) || ["1A", "1B", "1C"].includes(String(v?.tier || "").toUpperCase());
}

function _secondarySection(panel) {
  state.reports.secondary_findings = state.reports.secondary_findings || {};
  const section = state.reports.secondary_findings[panel] || {};
  const selected = Array.isArray(section.selected) ? section.selected.map(String) : [];
  const dismissed = Array.isArray(section.dismissed) ? section.dismissed.map(String) : [];
  return { selected, dismissed };
}

function _legacyPanelSelected(id, panel) {
  const m = state.reports.panels && state.reports.panels[id];
  return (m && m[panel]) === "V";
}

function _legacyPanelDismissed(id, panel) {
  const m = state.reports.panels && state.reports.panels[id];
  return (m && m[panel]) === "0";
}

function _secondaryPanelsForVariant(id) {
  const categories = state.data?.categories || {};
  return SECONDARY_PANEL_DEFS
    .map(def => def.key)
    .filter(panel => (categories[panel] || []).map(String).includes(String(id)));
}

function isSecondarySelected(id, panel) {
  const panels = _secondaryPanelsForVariant(id);
  if (!panels.length && panel) panels.push(panel);

  // A secondary finding has one global review state even when its gene is
  // present in multiple panels. Explicit dismissal wins so previously saved
  // ACMG/stroke disagreements resolve to the reviewer's opt-out.
  if (panels.some(key => {
    const section = _secondarySection(key);
    return section.dismissed.includes(String(id)) || _legacyPanelDismissed(id, key);
  })) return false;
  if (panels.some(key => {
    const section = _secondarySection(key);
    return section.selected.includes(String(id)) || _legacyPanelSelected(id, key);
  })) return true;
  const v = state.data?.variants?.[id];
  return _isClinvarPlp(v);
}

// Secondary status is variant-specific across all retained panels. New state
// lives in reports.secondary_findings; reports.panels is read only for legacy
// V/0 values.
function getPanelStatus(id, panel) {
  return isSecondarySelected(id, panel) ? "✓" : "";
}

function setPanelStatus(id, panel, val) {
  state.reports.secondary_findings = state.reports.secondary_findings || {};
  const panels = _secondaryPanelsForVariant(id);
  if (!panels.length && panel) panels.push(panel);
  panels.forEach(key => {
    const current = _secondarySection(key);
    const selected = new Set(current.selected);
    const dismissed = new Set(current.dismissed);
    if (val) {
      selected.add(String(id));
      dismissed.delete(String(id));
    } else {
      selected.delete(String(id));
      dismissed.add(String(id));
    }
    state.reports.secondary_findings[key] = {
      selected: Array.from(selected).sort(),
      dismissed: Array.from(dismissed).sort(),
    };
  });
  state.dirty = true;
  renderReportSections();
  renderCandidateSections();
  _syncStatusRadios(id, panel, val ? "✓" : "");
  updateSaveHint();
}

function getEdit(id, field) {
  const e = (state.reports.edits && state.reports.edits[id]) || {};
  return e[field];
}

function setEdit(id, field, val) {
  state.reports.edits = state.reports.edits || {};
  state.reports.edits[id] = state.reports.edits[id] || {};
  state.reports.edits[id][field] = val;
  state.dirty = true;
}

function _shortAcmgClass(value) {
  const cls = classifySignificance(value);
  return ({ "sig-p": "P", "sig-lp": "LP", "sig-vus": "VUS", "sig-lb": "LB", "sig-b": "B" })[cls]
      || String(value || "—");
}

function _acmgSourceHint(id, v) {
  const manual = getEdit(id, "ACMG_classification");
  const manualScore = getEdit(id, "ACMG_score");
  const manualCriteria = getEdit(id, "ACMG_criteria");
  const geneBe = v.genebe_acmg_class;
  const inHouse = v.ACMG_classification;
  let source = "in-house";
  const lines = [];
  const hasManualEdit = [manual, manualScore, manualCriteria]
    .some(value => value !== null && value !== undefined && value !== "");
  if (hasManualEdit) source = "manual";
  else if (geneBe !== null && geneBe !== undefined && geneBe !== "") source = "GeneBe";
  lines.push(`source: ${source}`);
  if (source === "manual") {
    lines.push(`in-house: ${_shortAcmgClass(inHouse)}`);
    lines.push(`GeneBe: ${_shortAcmgClass(geneBe)}`);
  } else if (source === "GeneBe") {
    lines.push(`in-house: ${_shortAcmgClass(inHouse)}`);
  } else {
    lines.push(`GeneBe: ${_shortAcmgClass(geneBe)}`);
  }
  return _annotationHint("ACMG source", lines, null, {
    className: "acmg-source-hint", dataId: id,
  });
}

function _refreshAcmgSourceHints(id) {
  const v = state.data?.variants?.[id];
  if (!v) return;
  const replacement = _acmgSourceHint(id, v);
  document.querySelectorAll(`.acmg-source-hint[data-id="${CSS.escape(id)}"]`).forEach(el => {
    el.outerHTML = replacement;
  });
}

function _syncEditControls(id, field, val, source = null) {
  const selectorByField = {
    comment: ".variant-comment",
    ACMG_score: ".acmg-score",
    ACMG_criteria: ".acmg-crit",
  };
  const selector = selectorByField[field];
  if (!selector) return;
  document.querySelectorAll(`${selector}[data-id="${CSS.escape(id)}"]`).forEach(el => {
    if (el === source) return;
    el.value = val ?? "";
  });
}

function _syncVariantCheckboxes(selector, id, idx, checked, source = null) {
  document.querySelectorAll(
    `${selector}[data-id="${CSS.escape(id)}"][data-idx="${CSS.escape(String(idx))}"]`
  ).forEach(input => {
    if (input !== source) input.checked = checked;
  });
}

function renderVariantCard(v, id, dropdownKind, opts = {}) {
  const isPanel    = dropdownKind === "panel";
  const panelKey   = isPanel ? (opts.category || "") : "";
  const panelAttr  = isPanel ? `data-panel="${escapeAttr(panelKey)}"` : "";
  const curStatus  = isPanel ? getPanelStatus(id, panelKey) : getStatus(id);
  const options    = statusOptions(dropdownKind);
  const idxTxt     = opts.index ? `#${opts.index}` : "";

  if (!v) {
    // Marked variant no longer present in current variants payload
    const card = document.createElement("div");
    card.className = "variant-card missing";
    card.innerHTML = `
      <div class="variant-head">
        ${idxTxt ? `<span class="card-idx">${idxTxt}</span>` : ""}
        <span class="muted">⚠️ 此 variant 在最新分析結果中不存在</span>
        <span class="hgvs">${escapeHtml(id)}</span>
        ${_renderStatusRadio(id, curStatus, options, panelAttr)}
      </div>`;
    return card;
  }

  const urls = variantUrls(v);
  const card = document.createElement("div");
  card.className = "variant-card";
  card.dataset.inPanel = v.in_panel ? "true" : "false";
  v = _variantWithSelectedTranscript(v, id);

  const links = [
    `<a href="#" class="btn-igv" data-id="${escapeAttr(id)}" title="在 IGV 內檢視">IGV</a>`,
    `<a href="${urls.varsome}"  target="_blank" rel="noopener">Varsome</a>`,
    `<a href="${urls.franklin}" target="_blank" rel="noopener">Franklin</a>`,
    `<a href="${urls.genebe}"   target="_blank" rel="noopener">GeneBe</a>`,
    urls.omim ? `<a href="${urls.omim}"     target="_blank" rel="noopener">OMIM</a>` : "",
  ].join("");

  // ACMG priority: reviewer override > GeneBe (second-opinion post-processing)
  // > pipeline. GeneBe ACMG_CLASS / SCORE / CRITERIA tend to be tighter
  // calibrated than the pipeline's per-rule classifier, so use them
  // when present and fall back to pipeline otherwise.
  const firstNonBlank = (...values) => values.find(x => x !== null && x !== undefined && x !== "") ?? "";
  const editAcmgClass = firstNonBlank(getEdit(id, "ACMG_classification"),
                                      v.genebe_acmg_class, v.ACMG_classification);
  const editAcmgCrit  = firstNonBlank(getEdit(id, "ACMG_criteria"),
                                      v.genebe_acmg_criteria, v.ACMG_criteria);
  const editAcmgScore = firstNonBlank(getEdit(id, "ACMG_score"),
                                      v.genebe_acmg_score, v.ACMG_score);
  const editComment   = getEdit(id, "comment")             ?? "";

  const clinvarDate = formatClinvarDate(state.data?.clinvar_date);
  const clinvarHint = _annotationHint("ClinVar", clinvarDate
    ? [`version date: ${clinvarDate}`]
    : ["version date: 三級輸出未提供"]);
  const scoreHint = _annotationHint("Score", [
    "Total score = Variant score + Phenotype score",
    "Variant score：ACMG score 轉換成 0–100",
    "Phenotype score：100 ×（此 gene 命中的 HPO/panel 權重總和）÷（全部有效輸入的 HPO/panel 權重總和）",
  ]);
  const hgvsLabel = displaySnvHgvs(v, id);

  // Extras shown only when the user clicks the "More" button. Each row is
  // pushed only when the underlying field is present in the webdata, so a
  // sample with none of these fields has the More button hidden too.
  const extras = [];
  if (v.in_silico_prediction != null && v.in_silico_prediction !== "") {
    extras.push({ key: "In silico prediction",
                  html: fmtInSilico(v.in_silico_prediction),
                  hint: _annotationHint("In silico prediction", [
                    "格式：pathogenic tool count - VUS tool count - benign tool count",
                    "此為工具計數摘要，不等同 PP3/BP4 evidence strength；請以各工具註解的 calibrated threshold 為準",
                  ]) });
  }
  if (v.LoGoFunc != null && v.LoGoFunc !== "" && v.LoGoFunc !== "NA") {
    extras.push({ key: "LoGoFunc", text: String(v.LoGoFunc),
                  cls: classifyLoGoFunc(v.LoGoFunc),
                  hint: _annotationHint("LoGoFunc", [
                    "GOF/LOF: model-predicted functional direction",
                    "*：probability 超過作者定義的 gene-level significance cutoff",
                    "尚無通用 PP3/BP4 strength calibration",
                  ], IN_SILICO_REFERENCES.logofunc) });
  }
  if (v.MaxEntScan_diff != null && v.MaxEntScan_diff !== "") {
    extras.push({ key: "MaxEntScan", text: fmtNum(v.MaxEntScan_diff),
                  cls: classifyMaxEntScan(v.MaxEntScan_diff),
                  hint: _annotationHint("MaxEntScan difference", [
                    "顯示 alternate 與 reference splice-site score 的差值",
                    "尚無可靠的通用 PP3/BP4 strength calibration，因此以黃色標示為 contextual evidence",
                  ], IN_SILICO_REFERENCES.maxentscan) });
  }
  if (v.PDIVAS_score != null && v.PDIVAS_score !== "") {
    extras.push({ key: "PDIVAS", text: fmtNum(v.PDIVAS_score),
                  cls: classifyPDIVAS(v.PDIVAS_score),
                  hint: _annotationHint("PDIVAS", [
                    "80%-sensitivity benchmark: ≥ 0.501",
                    "95%-sensitivity benchmark: ≥ 0.082",
                    "適用 deep-intronic splice-altering variants；不是 PP3/BP4 calibration",
                  ], IN_SILICO_REFERENCES.pdivas) });
  }
  // The in-silico predictor row follows the reviewer-defined order in
  // IN_SILICO_TOOLS; its first three populated values are shown directly and
  // any remaining predictors are rendered under More.
  if (v.loftee_hc || v.loftee_filter || v.loftee_flags) {
    const parts = [v.loftee_hc, v.loftee_filter, v.loftee_flags]
      .filter(Boolean).join(" / ");
    extras.push({ key: "LOFTEE", text: parts,
                  hint: _annotationHint("LOFTEE", [
                    "HC 表示 high-confidence predicted loss-of-function",
                    "此欄是 LoF transcript QC/filter，不是 PP3/BP4 score",
                  ]) });
  }
  const extrasHtml = extras.map(x => {
    const valHtml = x.html != null ? x.html : escapeHtml(x.text);
    return `<span class="k">${escapeHtml(x.key)}${x.hint || ""}</span>`
         + `<span class="v ${x.cls || ''}">${valHtml}</span>`;
  }).join("");

  // Score line: total (variant + pheno). Variant score is ACMG_POINTS
  // clamped to [-10, 10] then mapped to 0–100 (in backend); pheno
  // score is the in-house gene-level score (0–100). Total = sum, may
  // exceed 100 by design — the parens make it clear it's a composition.
  const _hasNum = x => x !== null && x !== undefined && x !== "" && Number.isFinite(Number(x));
  const _i = x => _hasNum(x) ? fmtInt(x) : "—";
  const scoreLine = (() => {
    const t = v.total_score, g = v.geno_score, p = v.pheno_score;
    if (![t, g, p].some(_hasNum)) return "—";
    return `${_i(t)} (${_i(g)} + ${_i(p)})`;
  })();
  const fmtScoreRank = (score, rank) => {
    if (!_hasNum(score) && !_hasNum(rank)) return "—";
    const s = _hasNum(score) ? fmtInt(score) : "—";
    return _hasNum(rank) ? `${s} (rank ${Number(rank)})` : s;
  };

  card.innerHTML = `
    <div class="variant-head">
      ${idxTxt ? `<span class="card-idx">${idxTxt}</span>` : ""}
      ${_renderStatusRadio(id, curStatus, options, panelAttr)}
      <span class="hgvs">${v.clinvar_upgrade ? `<span class="clinvar-upgrade-arrow" title="ClinVar 升級">${escapeHtml(v.clinvar_upgrade)}</span> ` : ""}${escapeHtml(hgvsLabel)}<button class="btn-copy" data-copy="${escapeAttr(hgvsLabel)}" title="複製 HGVS">${COPY_ICON_SVG}</button>${renderTranscriptPicker(v, id)} <span class="variant-tag">([${escapeHtml(state.data?.genome_build || "hg38")}] ${escapeHtml(id)}<button class="btn-copy" data-copy="${escapeAttr(id)}" title="複製 chr-pos-ref-alt">${COPY_ICON_SVG}</button>)</span></span>
      <span class="ext-links">${links}</span>
    </div>
    ${renderVariantBadges(v)}
    <div class="comment-row">
      <label>Comment:
        <input class="variant-comment" data-id="${escapeAttr(id)}" type="text" value="${escapeAttr(editComment)}" />
      </label>
    </div>
    <div class="info-grid">
      <div>
        <span class="k">Score${scoreHint}</span><span class="v">${escapeHtml(scoreLine)}</span>
        <span class="k">Exomiser</span><span class="v">${escapeHtml(fmtScoreRank(v.total_score_exomiser_variant, v.rank_exomiser_variant))}</span>
        <span class="k">LIRICAL</span><span class="v">${escapeHtml(fmtScoreRank(v.lirical_variant_score, v.rank_lirical_variant))}</span>
      </div>
      <div>
        <span class="k">Zygosity</span><span class="v">${fmtTxt(v.zygosity)}</span>
        <span class="k">Read depth (VAF)</span><span class="v ${v.low_depth ? "sig-lp" : ""}" title="${v.low_depth ? `Low DP (DP ${v.depth || "?"}) — 建議 IGV / Sanger 確認` : ""}">${escapeHtml(fmtAdVaf(v.AD, v.alt_af))}</span>
        <span class="k">Consequence</span><span class="v">${_renderConsequenceCell(v.Consequence)}</span>
        <div class="more-extras hidden">
          <span class="k">Exon / Intron</span><span class="v">${fmtExonIntron(v)}</span>
          <span class="k">Phase</span><span class="v">${fmtPhase(v)}</span>
        </div>
      </div>
      <div>
        <span class="k">ClinVar${clinvarHint}</span><span class="v ${classifySignificance(v.CLNSIG) || ""}">${escapeHtml(formatClinvar(v.CLNSIG, v.CLNSIGCONF, v.clinvar_stars))}${v.clinvar_upgrade && v.CLNSIG_old ? ` <span class="clinvar-old" title="原 ClinVar 分類">(was: ${escapeHtml(formatClinvar(v.CLNSIG_old, v.CLNSIGCONF_old, v.clinvar_stars_old))})</span>` : ""}</span>
        <span class="k">ACMG${_acmgSourceHint(id, v)}</span>
        <span class="acmg-class-row">
          <select class="acmg-class ${classifySignificance(editAcmgClass) || ""}" data-id="${escapeAttr(id)}">
            <option value=""                       ${editAcmgClass === ""                      ? "selected" : ""}>—</option>
            <option value="Pathogenic"             ${editAcmgClass === "Pathogenic"            ? "selected" : ""}>Pathogenic</option>
            <option value="Likely pathogenic"      ${editAcmgClass === "Likely pathogenic"     ? "selected" : ""}>Likely pathogenic</option>
            <option value="Uncertain significance" ${editAcmgClass === "Uncertain significance"? "selected" : ""}>VUS</option>
            <option value="Likely benign"          ${editAcmgClass === "Likely benign"         ? "selected" : ""}>Likely benign</option>
            <option value="Benign"                 ${editAcmgClass === "Benign"                ? "selected" : ""}>Benign</option>
          </select>
          <span class="acmg-paren">(</span>
          <input class="acmg-score" data-id="${escapeAttr(id)}" type="text" value="${escapeAttr(editAcmgScore)}" />
          <span class="acmg-paren">)</span>
        </span>
        <textarea class="acmg-crit" data-id="${escapeAttr(id)}" rows="2">${escapeHtml(editAcmgCrit)}</textarea>
      </div>
      <div>
        ${(() => {
          // Only render tools that actually have a numeric value;
          // first N → primary, rest → More.
          const populated = IN_SILICO_TOOLS.filter(t => _hasNum(v[t.scoreField]));
          const primary = populated.slice(0, IN_SILICO_PRIMARY_COUNT)
            .map(t => _renderInSilicoCell(v, t)).join("");
          const secondary = populated.slice(IN_SILICO_PRIMARY_COUNT)
            .map(t => _renderInSilicoCell(v, t)).join("");
          const more = secondary + extrasHtml;
          return primary + (more ? `<div class="more-extras hidden">${more}</div>` : "");
        })()}
      </div>
      <div>
        <span class="k">AF</span><span class="v">${fmtNum(v.AF, 5)}</span>
        <span class="k">AF_eas</span><span class="v">${fmtNum(v.AF_eas, 5)}</span>
        <span class="k">AF_nckuh</span><span class="v">${fmtNum(v.inhouse_af, 5)}${
          (v.inhouse_af != null && v.inhouse_ac != null && v.inhouse_an != null)
            ? ` (${v.inhouse_ac}/${v.inhouse_an})` : ""}</span>
      </div>
    </div>
    <button class="btn-more" type="button">▾ More</button>
    <div class="more-extras hidden">${Number.isFinite(Number(v.TG_eas_af))
      ? `<div><span class="k">1000G EAS</span><span class="v">${fmtNum(v.TG_eas_af, 5)}</span></div>`
      : ""}${renderManeAll(v)}</div>
    ${renderDiseaseList(v, id, !!opts.diseaseCheckbox)}
  `;

  return card;
}

// ---- helpers used by renderVariantCard (Phase 4) -----------------

// GIAB stratification label → badge text + tooltip. Labels come from the
// TSV GIAB_STRATA column (written by scripts/annotate_giab_strata.py, which
// maps each BED to a label via strata_manifest.json). Unknown labels fall
// back to the raw label so adding a new stratum needs no frontend change.
const GIAB_STRATA_DISPLAY = {
  homopolymer:     { display: "Homopolymer",    tip: "GIAB: homopolymer region — low-confidence indel calls" },
  tandem_repeat:   { display: "Tandem repeat",  tip: "GIAB: tandem / simple repeat region" },
  segdup:          { display: "Segdup",         tip: "GIAB: segmental duplication — possible mismapping / false positives" },
  low_mappability: { display: "Low mappability",tip: "GIAB: low-mappability region — reads hard to place uniquely" },
  gc_extreme:      { display: "GC extreme",     tip: "GIAB: extreme GC content — coverage / calling bias" },
  other_difficult: { display: "Other difficult",tip: "GIAB: other difficult region (MHC/KIR/VDJ, false dup, gaps, ...)" },
};

// Top-of-card chip row: TRANSCRIPT_TYPE / CALLERS / panel/ROH/blacklist
// hits / LOFTEE HC / GIAB strata. Empty-string entries are filtered out so
// the row hides itself when nothing is worth showing.
function renderVariantBadges(v) {
  const chips = [];
  if (v.transcript_type) {
    const cls = "badge-tx badge-" + v.transcript_type.toLowerCase().replace(/_/g, "-");
    chips.push(`<span class="badge ${cls}" title="Transcript type">${escapeHtml(v.transcript_type)}</span>`);
  }
  if (v.callers) {
    const cls = v.callers === "DV+HC" ? "badge-callers-both"
              : v.callers === "DV"    ? "badge-callers-dv"
              : v.callers === "HC"    ? "badge-callers-hc"
              :                          "badge-callers";
    chips.push(`<span class="badge ${cls}" title="Variant callers">${escapeHtml(v.callers)}</span>`);
  }
  if (v.in_panel)     chips.push(`<span class="badge badge-panel"     title="Gene is in the requested panel">In panel</span>`);
  if (v.in_roh)       chips.push(`<span class="badge badge-roh"       title="Variant falls inside an ROH region">In ROH</span>`);
  if (v.in_blacklist) chips.push(`<span class="badge badge-blacklist" title="Variant or gene flagged on the QC blacklist">⚠ Blacklist</span>`);
  if (v.loftee_hc === "HC") {
    chips.push(`<span class="badge badge-loftee-hc" title="LOFTEE high-confidence LoF">LOFTEE HC</span>`);
  }
  // GIAB genome-stratification flags — difficult regions where short-read
  // calls are less reliable (homopolymers, repeats, segdups, ...). One
  // amber badge per label so reviewers treat the call with care.
  for (const label of (v.giab_strata || [])) {
    const meta = GIAB_STRATA_DISPLAY[label] || { display: label, tip: "GIAB difficult region" };
    chips.push(`<span class="badge badge-giab" title="${escapeAttr(meta.tip)}">${escapeHtml(meta.display)}</span>`);
  }
  // Right-aligned "搜尋同基因" — lists every SNV/Indel + CNV/SV that
  // touches this gene. Mainly for spotting compound-het / mixed-mode
  // hits when the AR diagnosis is on the table.
  const sameGeneBtn = v.gene_symbol
    ? `<button class="same-gene-btn" data-gene="${escapeAttr(v.gene_symbol)}" type="button" title="列出此基因的所有 SNV/Indel + CNV/SV 變異">搜尋同基因</button>`
    : "";
  if (!chips.length && !sameGeneBtn) return "";
  return `<div class="variant-badges">
    <span class="variant-badges-chips">${chips.join("")}</span>
    ${sameGeneBtn}
  </div>`;
}

// "trans / cis / unphased" — show phase group too when present so the
// user can see which co-segregating variants share the same haplotype.
// VEP EXON/INTRON come through as "current/total" (e.g. "6/10"); show
// whichever the variant falls in, "—" when both are blank.

// Consequence cell — VEP joins multiple SO terms with "&"
// (e.g. "missense_variant&splice_region_variant"). Show only the
// first term by default + a tiny ▾ button to expand the rest;
// reviewers don't usually need to see every joined consequence and
// the column is narrow.
function _renderConsequenceCell(raw) {
  const s = (raw || "").toString().trim();
  if (!s) return "—";
  const parts = s.split("&").map(p => p.trim()).filter(Boolean);
  if (parts.length <= 1) return escapeHtml(s);
  const first = parts[0];
  const rest  = parts.slice(1).join(" & ");
  return `<span class="consequence-multi">
    <span class="consequence-first">${escapeHtml(first)}</span>
    <button class="consequence-toggle" type="button" title="展開其餘 ${parts.length - 1} 個 consequence">▾</button>
    <span class="consequence-rest hidden"> &amp; ${escapeHtml(rest)}</span>
  </span>`;
}

// "舊格式，請重跑新版 pipeline" banner — shown inside the SNV card
// whenever sample_loader rejected the TSV (old-format detection
// raised OldFormatError). Reviewer state (status / edits / phenotype)
// is preserved on the backend; just the variant list is empty.
function _renderSnvTsvErrorBanner() {
  const card = document.getElementById("card-snv");
  if (!card) return;
  const existing = card.querySelector(".snv-tsv-error");
  const msg = (state.data?.snv_tsv_error || "").trim();
  if (!msg) {
    if (existing) existing.remove();
    return;
  }
  if (existing) { existing.textContent = msg; return; }
  const banner = document.createElement("div");
  banner.className = "snv-tsv-error";
  banner.innerHTML = `<strong>⚠ ${escapeHtml(msg)}</strong>
    <div class="muted" style="margin-top:4px">reviewer 編輯狀態已保留；重跑後直接 reload 即可。</div>`;
  card.insertBefore(banner, card.firstChild);
}

function fmtExonIntron(v) {
  const cleanRank = value => {
    const text = String(value || "").trim();
    return ["", ".", "-", "NA", "N/A"].includes(text.toUpperCase()) ? "" : text;
  };
  const e = cleanRank(v.exon);
  if (e) return `exon ${e}`;
  const i = cleanRank(v.intron);
  if (i) return `intron ${i}`;
  return "—";
}

function fmtPhase(v) {
  const result = (v.phase_result || "").trim();
  const group  = (v.phase_group  || "").trim();
  if (!result && !group) return "—";
  if (!result) return group;
  if (!group)  return result;
  return `${result} <span class="muted">(PG=${escapeHtml(group)})</span>`;
}

// Inline <details> block listing transcript annotations carried by the
// grouped TSV rows (new pipeline) or legacy MANE_ALL payload.
function renderManeAll(v) {
  const rows = Array.isArray(v.transcript_options) && v.transcript_options.length
    ? v.transcript_options
    : (Array.isArray(v.MANE_ALL) ? v.MANE_ALL : []);
  if (!rows.length) return "";
  const cells = rows.map(r => `
    <tr>
      <td><span class="badge badge-${(r.transcript_type || "").toLowerCase().replace(/_/g, "-")}">${escapeHtml(r.transcript_type || "")}</span></td>
      <td class="mane-tx">${escapeHtml(_maneDisplayTranscript(r))}</td>
      <td>${escapeHtml(_maneDisplayHgvs(r, "HGVS_C") || _maneDisplayHgvs(r, "hgvs_c"))}</td>
      <td>${escapeHtml(_maneDisplayHgvs(r, "HGVS_P") || _maneDisplayHgvs(r, "hgvs_p"))}</td>
      <td>${escapeHtml(r.Consequence || r.consequence || "")}</td>
    </tr>`).join("");
  return `
    <details class="mane-all">
      <summary>Transcripts (${rows.length})</summary>
      <table class="mane-table">
        <thead><tr><th>Type</th><th>Transcript</th><th>HGVS.c</th><th>HGVS.p</th><th>Consequence</th></tr></thead>
        <tbody>${cells}</tbody>
      </table>
    </details>`;
}

// Manual variant cards live in state.reports.manual_variants instead of
// state.data.variants — they have no upstream call, just three free-text
// fields the user fills in (typically for CNVs that the SNV pipeline
// doesn't touch). Sync targets:
//   position → Causative_variant or Other_variant cell in xlsx
//   disease  → appended into the Disease cell
//   comment  → kept on the report JSON only (helps the user, not exported)
function renderManualVariantCard(m) {
  const card = document.createElement("div");
  card.className = "variant-card variant-card-manual";
  card.dataset.mid = m.id;
  const pos = m.position || "";
  card.innerHTML = `
    <div class="manual-row manual-row-pos">
      <input class="manual-position" data-mid="${escapeAttr(m.id)}"
             placeholder="點位（如 chr2:123456-654321 del）"
             value="${escapeAttr(pos)}" />
      <button class="btn-copy" data-copy="${escapeAttr(pos)}" title="複製點位">${COPY_ICON_SVG}</button>
      <a class="btn-link" href="https://www.deciphergenomics.org/" target="_blank" rel="noopener">Decipher</a>
      <button class="btn-remove-manual" data-mid="${escapeAttr(m.id)}" title="刪除這個 variant" type="button">×</button>
    </div>
    <div class="manual-row">
      <label>Comment:
        <input class="manual-comment" data-mid="${escapeAttr(m.id)}"
               placeholder="備註"
               value="${escapeAttr(m.comment || "")}" />
      </label>
    </div>
    <div class="manual-row">
      <label>Disease:
        <input class="manual-disease" data-mid="${escapeAttr(m.id)}"
               placeholder="疾病名稱（可包含遺傳模式 e.g. (AD)）"
               value="${escapeAttr(m.disease || "")}" />
      </label>
    </div>
  `;
  return card;
}

function addManualVariant(status) {
  if (!Array.isArray(state.reports.manual_variants)) state.reports.manual_variants = [];
  const id = "m_" + Date.now() + "_" + Math.floor(Math.random() * 1e6);
  state.reports.manual_variants.push({
    id, status, position: "", comment: "", disease: "",
  });
  state.dirty = true;
  renderAll();
}

function removeManualVariant(mid) {
  if (!Array.isArray(state.reports.manual_variants)) return;
  state.reports.manual_variants = state.reports.manual_variants.filter(m => m.id !== mid);
  state.dirty = true;
  renderAll();
}

function updateManualVariant(mid, field, value) {
  const m = (state.reports.manual_variants || []).find(x => x.id === mid);
  if (!m) return;
  m[field] = value;
  state.dirty = true;
}

function diseaseAssociationSummary(a) {
  if (!a) return "";
  const name = a.display_name || "";
  const mim = a.phenotype_mim ? ` (${a.phenotype_mim})` : "";
  const inh = a.inheritance ? `(${a.inheritance})` : "";
  return `${name}${mim}${inh}`.trim();
}

function diseaseAssociationDetail(a) {
  if (!a) return "";
  if (a.detail) return String(a.detail);
  const lines = [diseaseAssociationSummary(a)].filter(Boolean);
  const evidence = Array.isArray(a.evidence) ? a.evidence.filter(Boolean) : [];
  const sources = Array.isArray(a.sources) ? a.sources.filter(Boolean) : [];
  if (evidence.length) lines.push(`來源：${evidence.join("；")}`);
  else if (sources.length) lines.push(`來源：${sources.join("；")}`);
  if (a.mondo_id) lines.push(`MONDO：${a.mondo_id}`);
  return lines.join("\n");
}

function diseaseSourceBadges(a) {
  const labels = Array.isArray(a?.evidence) && a.evidence.length
    ? a.evidence
    : (Array.isArray(a?.sources) ? a.sources : []);
  return labels.filter(Boolean).slice(0, 3).map(label =>
    `<span class="disease-source-badge">${escapeHtml(label)}</span>`
  ).join("");
}

function hasOmimDescriptionText(a, d) {
  if (a && Object.prototype.hasOwnProperty.call(a, "needs_description")) {
    return !a.needs_description;
  }
  return String(d || "").split("\n").some((line, idx) => idx > 0 && line.trim());
}

function renderDiseaseList(v, id, withCheckbox) {
  const rows = [];
  const picked = (getEdit(id, "report_diseases") || {});
  const associations = Array.isArray(v.disease_associations) ? v.disease_associations : null;

  if (associations && associations.length) {
    const orderedAssociations = associations
      .map((association, originalIndex) => ({ association, originalIndex }))
      .sort((left, right) => {
        const leftSupplemental = left.association?.source_kind === "omim" ? 0 : 1;
        const rightSupplemental = right.association?.source_kind === "omim" ? 0 : 1;
        return leftSupplemental - rightSupplemental || left.originalIndex - right.originalIndex;
      })
      .map(({ association }) => association);
    let supplementalStarted = false;
    for (const a of orderedAssociations) {
      const d = diseaseAssociationDetail(a);
      if (!d || d === "NA") continue;
      const summary = (diseaseAssociationSummary(a) || String(d).split("\n")[0] || "").slice(0, 120);
      const idx = Number(a.omim_slot || 0);
      const canPick = Boolean(idx && a.source_kind === "omim");
      const checked = canPick && picked[idx] ? "checked" : "";
      const checkbox = withCheckbox && canPick
        ? `<input type="checkbox" class="disease-pick" data-id="${escapeAttr(id)}" data-idx="${idx}" ${checked} title="報告要發這個疾病" />`
        : "";
      const needsDescription = a.source_kind === "omim" && !hasOmimDescriptionText(a, d);
      const star = needsDescription ? `<span class="disease-needs-description-star" aria-label="OMIM description missing">*</span>` : "";
      const badges = a.source_kind === "omim" ? diseaseSourceBadges(a).replace('<span class="disease-source-badge">OMIM</span>', "") : diseaseSourceBadges(a);
      const isOmim = a.source_kind === "omim";
      const isFirstSupplemental = !isOmim && !supplementalStarted;
      if (!isOmim) supplementalStarted = true;
      const extraClass = isOmim
        ? " disease-row-omim"
        : ` disease-row-supplemental${isFirstSupplemental ? " disease-row-supplemental-first" : ""}`;
      rows.push(`
        <details class="disease-row${extraClass}">
          <summary>${checkbox}<span class="disease-summary-text">${escapeHtml(summary)}${star}</span><span class="disease-source-badges">${badges}</span></summary>
          <div class="disease-detail">${escapeHtml(String(d))}<button type="button" class="disease-collapse">▴ 收合</button></div>
        </details>`);
    }
  } else {
    for (let i = 1; i <= 5; i++) {
      const d = v[`Disease${i}`];
      if (!d || d === "NA") continue;
      const summary = (String(d).split("\n")[0] || "").slice(0, 120);
      const checked = picked[i] ? "checked" : "";
      const checkbox = withCheckbox
        ? `<input type="checkbox" class="disease-pick" data-id="${escapeAttr(id)}" data-idx="${i}" ${checked} title="報告要發這個疾病" />`
        : "";
      const star = hasOmimDescriptionText(null, d) ? "" : `<span class="disease-needs-description-star" aria-label="OMIM description missing">*</span>`;
      rows.push(`
        <details class="disease-row disease-row-omim">
          <summary>${checkbox}<span class="disease-summary-text">${escapeHtml(summary)}${star}</span></summary>
          <div class="disease-detail">${escapeHtml(String(d))}<button type="button" class="disease-collapse">▴ 收合</button></div>
        </details>`);
    }
  }
  if (!rows.length) return "";
  return `<div class="disease-list">${rows.join("")}</div>`;
}

// ---------- Render: sections ---------------------------------------

const SECONDARY_PANEL_DEFS = [
  { key: "acmg_sf",            title: "ACMG SF" },
  { key: "stroke",             title: "中風相關基因" },
  { key: "carrier",            title: "Carrier screening" },
];

const REPORT_SECTION_DEFS = [
  { el: "sec-causative", title: "Causative variants", match: id => getStatus(id) === "1", dropdown: "candidate", defaultOpen: true, diseaseCheckbox: true, manualStatus: "1" },
  { el: "sec-other",     title: "Other variants",     match: id => getStatus(id) === "2", dropdown: "candidate", defaultOpen: true, diseaseCheckbox: true, manualStatus: "2" },
  { el: "sec-candidate", title: "Candidate variants", match: id => _statusHas(getStatus(id), "C"), dropdown: "candidate", defaultOpen: true, diseaseCheckbox: true, manualStatus: "C" },
  // ACMG SF / health-screening panels / PharmCat all live inside the
  // Secondary findings collapsible group in the HTML; they render
  // the same way as before, just nested in a different container.
  ...SECONDARY_PANEL_DEFS.map(def => ({
    el: `sec-${def.key.replace(/_/g, "-")}`,
    title: def.title,
    category: def.key,
    dropdown: "panel",
    defaultOpen: def.key === "acmg_sf",
    diseaseCheckbox: true,
  })),
];

// Tier sections per 三級輸出計畫.md §2.3. Backend categorises each variant
// into 1A / 1B / 2 / 3 / 4 / 5 based on ClinVar / LOFTEE / ACMG points.
const CANDIDATE_SECTION_DEFS = [
  { el: "cat-tier-1a", title: "1A — ClinVar P/LP ≥ 1★",        category: "1A", dropdown: "candidate", tier: "1A", defaultOpen: true },
  { el: "cat-tier-1b", title: "1B — Frameshift / Nonsense (LOFTEE HC)", category: "1B", dropdown: "candidate", tier: "1B", defaultOpen: true },
  { el: "cat-tier-1c", title: "1C — Predicted suspect",        category: "1C", dropdown: "candidate", tier: "1C", defaultOpen: true },
  { el: "cat-tier-2",  title: "2 — Other",                     category: "2",  dropdown: "candidate", tier: "2"  },
  ...SECONDARY_PANEL_DEFS.map(def => ({
    el: `cat-${def.key.replace(/_/g, "-")}-c`,
    title: def.title,
    category: def.key,
    dropdown: "panel",
    defaultOpen: def.key === "acmg_sf",
  })),
];

// Look up a variant id across all four maps (SNV / Mito / CNV / SV)
// so the Causative / Other report sections can render variants no
// matter which card the reviewer flipped the status on. Returns
// {v, kind} or {v: null, kind: null}.
function lookupAnyVariant(id) {
  const d = state.data || {};
  // Mito variants can share the same chrM-pos-ref-alt id with raw SNV
  // rows. Report sections must keep the mito-specific card and m.HGVS
  // display instead of accidentally resolving to the SNV transcript row.
  let v = (d.mito_variants || {})[id];
  if (v) return { v, kind: "mito" };
  v = (d.variants || {})[id];
  if (!v) v = (state.snvSearchVariants || {})[id];
  if (v) return { v, kind: "snv" };
  v = _cnvSvVariantById(id);
  if (v) return { v, kind: v.source === "sv" ? "sv" : "cnv" };
  return { v: null, kind: null };
}

function _annotSvSortScore(v) {
  const n = Number(v?.cnv_sv_sort_score);
  return Number.isFinite(n) ? n : -Infinity;
}

function _reportVariantSortScore(id) {
  const { v, kind } = lookupAnyVariant(id);
  if (kind === "cnv" || kind === "sv") return _annotSvSortScore(v);
  const n = Number(v?.total_score);
  return Number.isFinite(n) ? n : -Infinity;
}

function _reportVariantKindRank(id) {
  const kind = lookupAnyVariant(id).kind || "";
  if (kind === "snv") return 0;
  if (kind === "cnv" || kind === "sv") return 1;
  if (kind === "mito") return 2;
  return 3;
}

function idsForReportSection(def) {
  const d = state.data || {};
  const known = [
    ...Object.keys(d.variants      || {}),
    ...Object.keys(d.mito_variants || {}),
    ...Object.keys(d.cnv_variants  || {}),
    ...Object.keys(d.sv_variants   || {}),
  ];
  const reported      = Object.keys(state.reports.status || {});
  const panelReported = Object.keys(state.reports.panels || {});
  const secondaryReported = Object.values(state.reports.secondary_findings || {})
    .flatMap(section => [
      ...(Array.isArray(section?.selected) ? section.selected : []),
      ...(Array.isArray(section?.dismissed) ? section.dismissed : []),
    ]);
  const all = Array.from(new Set([...known, ...reported, ...panelReported, ...secondaryReported]));

  if (def.match) {
    // Causative / Other / Candidate report sections: keep the visual order
    // stable by variant type (SNV/Indel → CNV/SV → Mito), then use each
    // adapter's score inside that type. Same-gene variants cluster within
    // their type.
    // The gene with the highest-scored variant leads; its lower-
    // scored siblings get pulled up directly behind it instead of
    // scattering down the list. Manual entries (no gene_symbol)
    // stay put as singleton clusters.
    const sorted = all.filter(def.match).sort((a, b) => {
      const rankDiff = _reportVariantKindRank(a) - _reportVariantKindRank(b);
      if (rankDiff) return rankDiff;
      return _reportVariantSortScore(b) - _reportVariantSortScore(a);
    });
    const groups = new Map();
    for (const id of sorted) {
      const { v, kind } = lookupAnyVariant(id);
      const key = `${kind || "unknown"}:${v?.gene_symbol || `__${id}`}`;
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(id);
    }
    return Array.from(groups.values()).flat();
  }
  if (def.category) {
    const inCat = new Set(state.data.categories?.[def.category] || []);
    return all.filter(id =>
      inCat.has(id)
      && _isSecondaryEligible(id)
      && isSecondarySelected(id, def.category)
    );
  }
  return [];
}

function idsForCandidateSection(def, { ignoreInPanelOnly = false } = {}) {
  const ids = state.data.categories?.[def.category] || [];
  if (def.dropdown === "panel" && def.category) {
    return ids.filter(id => _isSecondaryEligible(id));
  }
  return ids.filter(id => _passesMainSnvDisplayFilters(
    state.data.variants?.[id],
    { ignoreInPanelOnly },
  ));
}

function candidateIdsForSection(def) {
  const countIds = idsForCandidateSection(def, { ignoreInPanelOnly: true });
  if (def.dropdown === "panel" && def.category) {
    return { displayIds: countIds, countIds };
  }
  const displayIds = countIds.filter(id =>
    _passesMainSnvDisplayFilters(state.data.variants?.[id])
  );
  return { displayIds, countIds };
}

function _passesMainSnvDisplayFilters(v, { ignoreInPanelOnly = false, ignoreDiseaseAssociated = false } = {}) {
  if (!v) return false;
  if (!ignoreDiseaseAssociated
      && document.getElementById("filter-disease-associated")?.checked
      && !v.disease_associated) return false;
  if (!ignoreInPanelOnly
      && document.getElementById("filter-in-panel-only")?.checked
      && !v.in_panel) return false;
  if (!document.getElementById("filter-vaf")?.checked) {
    const vaf = _numericValue(v.alt_af);
    if ((vaf != null && vaf < 0.2) || _isReferenceZygosity(v)) return false;
  }
  if (!document.getElementById("filter-impact-modifier")?.checked
      && String(v.impact || "").toUpperCase() === "MODIFIER"
      && !_isClinvarPlp(v)) return false;
  return true;
}

function renderBlock(def, ids, openKey, countIds = ids) {
  const host = document.getElementById(def.el);
  host.innerHTML = "";
  host.dataset.openKey = openKey;

  const isPanel = def.dropdown === "panel" && !!def.category;
  // Skip X-marked variants — panel sections use panel-specific X, others use global
  const visibleIds = ids.filter(id => isPanel
    ? getPanelStatus(id, def.category) !== "X"
    : getStatus(id) !== "X");
  const countVisibleIds = countIds.filter(id => isPanel
    ? getPanelStatus(id, def.category) !== "X"
    : getStatus(id) !== "X");
  const wasOpen = toggledBlocks.has(def.el)
    ? host.dataset.wasOpen === "1"
    : (!!def.defaultOpen || (isPanel && currentSampleTestType() === "TITAN-WGS"));
  host.dataset.wasOpen = wasOpen ? "1" : "0";

  const header = document.createElement("div");
  header.className = "block-header" + (wasOpen ? " open" : "");
  // Counts: "In panel X / Total Y" so reviewers can see at a glance how
  // much of the section overlaps the requested panel.
  const inPanelCount = countVisibleIds.filter(
    id => state.data.variants?.[id]?.in_panel
  ).length;
  const countLabel = isPanel
    ? `Total ${countVisibleIds.length}`
    : `In panel ${inPanelCount} / Total ${countVisibleIds.length}`;
  header.innerHTML = `
    <span><span class="arrow"></span><span class="title">${escapeHtml(def.title)}</span></span>
    <span class="count">${escapeHtml(countLabel)}</span>`;
  host.appendChild(header);

  const body = document.createElement("div");
  body.className = "block-body" + (wasOpen ? " open" : "");
  const manuals = def.manualStatus
    ? (state.reports.manual_variants || []).filter(m => m.status === def.manualStatus)
    : [];
  if (isPanel && state.data?.secondary_pending) {
    body.innerHTML = `<div class="muted">載入中…</div>`;
  } else if (!visibleIds.length && !manuals.length && !def.manualStatus) {
    body.innerHTML = `<div class="muted">（無符合點位）</div>`;
  } else {
    visibleIds.forEach((id, i) => {
      // Dispatch by variant kind so mito and CNV/SV variants tagged
      // status=1/2 render with their proper cards (heteroplasmy /
      // copy-number / cytoband, etc.) instead of the SNV-shaped
      // "missing variant" placeholder. Falls back to the SNV card's
      // missing-variant warning when the id is reported but not in
      // any current payload.
      const { v, kind } = lookupAnyVariant(id);
      const opts = {
        category: def.category,
        index: i + 1,
        diseaseCheckbox: !!def.diseaseCheckbox,
      };
      let card;
      if (v && kind === "mito") {
        card = renderMitoCard(v, id, opts);
      } else if (v && (kind === "cnv" || kind === "sv")) {
        card = renderCnvSvCard(v, id, opts);
      } else {
        card = renderVariantCard(v, id, def.dropdown, opts);
      }
      body.appendChild(card);
    });
    manuals.forEach(m => body.appendChild(renderManualVariantCard(m)));
    if (def.manualStatus) {
      const addBtn = document.createElement("button");
      addBtn.type = "button";
      addBtn.className = "btn-add-manual";
      addBtn.dataset.status = def.manualStatus;
      addBtn.textContent = "＋ 新增 variant";
      body.appendChild(addBtn);
    }
    if (!visibleIds.length && !manuals.length) {
      const empty = document.createElement("div");
      empty.className = "muted";
      empty.textContent = "（無符合點位）";
      body.insertBefore(empty, body.firstChild);
    }
  }
  host.appendChild(body);

  header.addEventListener("click", () => {
    const open = body.classList.toggle("open");
    header.classList.toggle("open", open);
    host.dataset.wasOpen = open ? "1" : "0";
    toggledBlocks.add(def.el);
  });
}

function renderReportSections() {
  for (const def of REPORT_SECTION_DEFS) {
    renderBlock(def, idsForReportSection(def), def.el);
  }
  renderPharmcatBlock("sec-pharmcat");
  document.getElementById("report-sections").classList.remove("hidden");
  document.getElementById("save-row-mid")?.classList.remove("hidden");
}

function renderCandidateSections() {
  // Select the active tab before building cards. Hidden SNV tiers stay
  // empty until clicked, avoiding thousands of unnecessary DOM nodes.
  renderTierTabBar();
  for (const def of CANDIDATE_SECTION_DEFS) {
    if (def.tier && def.tier !== activeTierTab) {
      const host = document.getElementById(def.el);
      if (host) host.innerHTML = "";
      continue;
    }
    const { displayIds, countIds } = candidateIdsForSection(def);
    renderBlock(
      def,
      displayIds,
      def.el,
      countIds,
    );
  }
  renderPharmcatBlock("cat-pharmcat-c");
  document.getElementById("category-sections").classList.remove("hidden");
  updateInPanelCount();
  renderCnvSvTabBar();
  renderMitoTabBar();
  renderStrTabBar();
}

function renderActiveSnvTier() {
  const def = CANDIDATE_SECTION_DEFS.find(d => d.tier === activeTierTab);
  if (!def) {
    applyTierTabActive();
    return;
  }
  const { displayIds, countIds } = candidateIdsForSection(def);
  renderBlock(def, displayIds, def.el, countIds);
  applyTierTabActive();
}

// Build the SNV/Indel tier tab bar from the same defs / counts as the
// panels themselves. Each tab carries the tier title + 'In panel X /
// Total Y' so the collapsed view still surfaces those numbers. The
// active panel is whatever was active before, falling back to the
// first tier with any visible variants, falling back to 1A.
const TIER_ORDER = ["1A", "1B", "1C", "2"];
let activeTierTab = null;

function renderTierTabBar() {
  const bar = document.getElementById("tier-tab-bar");
  if (!bar) return;
  // Old-format TSV banner — sample_loader couldn't parse the SNV TSV
  // (lacks ACMG_CRITERIA). Render a single message in the SNV card so
  // the reviewer knows what to do; other side-channels (mito, CNV, …)
  // still render normally.
  _renderSnvTsvErrorBanner();

  const counts = {};
  for (const tier of TIER_ORDER) {
    const def = CANDIDATE_SECTION_DEFS.find(d => d.tier === tier);
    const { displayIds, countIds } = def
      ? candidateIdsForSection(def)
      : { displayIds: [], countIds: [] };
    const visible = countIds.filter(id => getStatus(id) !== "X");
    const displayTotal = displayIds.filter(id => getStatus(id) !== "X").length;
    const inPanel = visible.filter(
      id => state.data.variants?.[id]?.in_panel
    ).length;
    counts[tier] = { total: visible.length, inPanel, displayTotal };
  }

  // Pick the active tier: keep what was active if it still exists,
  // else first tier with variants, else 1A so the bar is never blank.
  if (!TIER_ORDER.includes(activeTierTab)) activeTierTab = null;
  if (!activeTierTab) {
    activeTierTab = TIER_ORDER.find(t => counts[t].displayTotal > 0) || "1A";
  }

  const titles = {
    "1A": "1A — ClinVar P/LP ≥ 1★",
    "1B": "1B — Frameshift / Nonsense",
    "1C": "1C — Predicted suspect",
    "2":  "2 — Other",
  };
  bar.innerHTML = TIER_ORDER.map(t => {
    const c = counts[t];
    const cls = "tier-" + t.toLowerCase();
    const active = t === activeTierTab ? " active" : "";
    return `<button type="button" class="tier-tab ${cls}${active}" data-tier="${t}">
              <span class="tier-tab-title">${escapeHtml(titles[t])}</span>
              <span class="tier-tab-count">In panel ${c.inPanel} / Total ${c.total}</span>
            </button>`;
  }).join("");

  applyTierTabActive();
}

function applyTierTabActive() {
  const panels = document.getElementById("tier-tab-panels");
  if (!panels) return;
  panels.querySelectorAll(".tier-panel").forEach(p => {
    p.classList.toggle("active", p.dataset.tier === activeTierTab);
  });
}

// CNV/SV tab bar: same tab UX as SNV, but the backend doesn't produce
// these tiers yet so every panel renders an empty placeholder. The
// structure stays so the next pipeline pass can drop variants in
// without touching the UI.
const CNV_SV_TIER_ORDER = ["CNV-1A", "CNV-1B", "SV-2A", "SV-2B"];
const CNV_SV_TITLES = {
  "CNV-1A": "1A CNV Clinical",
  "CNV-1B": "1B CNV Pathogenic",
  "SV-2A":  "2A SV Clinical",
  "SV-2B":  "2B SV Pathogenic",
};
const CNV_SV_TIER_CLASS = {
  "CNV-1A": "tier-cnv",
  "CNV-1B": "tier-cnv",
  "SV-2A":  "tier-sv",
  "SV-2B":  "tier-sv",
};
let activeCnvSvTab = null;

// Reads cnv_variants/sv_variants/cnv_categories/sv_categories from
// state.data and dispatches each variant id to the right tier panel
// renderer. Tier counts on the tab bar reflect the actual list size.
function _cnvSvVariantById(id) {
  return _cnvSvBaseVariantById(id)
      || _cnvSvVirtualParents()[id]
      || null;
}

function _cnvSvBaseVariantById(id) {
  return (state.data?.cnv_variants?.[id])
      || (state.data?.sv_variants?.[id])
      || null;
}

function _confirmedCnvSvMerges() {
  return Array.isArray(state.reports?.cnv_sv_merges) ? state.reports.cnv_sv_merges : [];
}

function _cnvSvMergeId(source, chrom, start, end, svType) {
  return `MERGED-${String(source || "cnv").toUpperCase()}-${chrom}-${start}-${end}-${String(svType || "").toUpperCase()}`;
}

function _cnvSvBuildParent(merge) {
  const segments = (merge.member_ids || []).map(_cnvSvBaseVariantById).filter(Boolean);
  if (segments.length < 2) return null;
  const chrom = segments[0].CHROM || "";
  const svType = String(segments[0].sv_type || "").toUpperCase();
  const source = String(merge.source || segments[0].source || "cnv").toLowerCase();
  if (!chrom || !["DEL", "DUP"].includes(svType)) return null;
  if (segments.some(v => v.CHROM !== chrom || String(v.sv_type || "").toUpperCase() !== svType)) return null;
  const start = Math.min(...segments.map(v => Number(v.POS)));
  const end = Math.max(...segments.map(v => Number(v.END)));
  if (!Number.isFinite(start) || !Number.isFinite(end)) return null;
  const rep = segments.slice().sort((a, b) => Number(b.cnv_sv_sort_score || -999) - Number(a.cnv_sv_sort_score || -999))[0];
  const bestRank = Math.max(...segments.map(v => Number(v.ranking_score)).filter(Number.isFinite));
  const bestScaledRank = Math.max(0, ...segments.map(v => Number(v.ranking_score_scaled)).filter(Number.isFinite));
  const bestPheno = Math.max(0, ...segments.map(v => Number(v.max_pheno_score)).filter(Number.isFinite));
  const bestCombined = Math.max(...segments.map(v => Number(v.cnv_sv_sort_score)).filter(Number.isFinite));
  const genes = [];
  const seenGenes = new Set();
  segments.forEach(seg => [...(seg.genes || []), ...(seg.genes_overflow || [])].forEach(g => {
    if (!g?.gene || seenGenes.has(g.gene)) return;
    seenGenes.add(g.gene);
    genes.push(g);
  }));
  return {
    ...rep,
    id: merge.id || _cnvSvMergeId(source, chrom, start, end, svType),
    source, CHROM: chrom, POS: start, END: end, length: end - start,
    gene_count: genes.length, genes, genes_overflow: [], genes_compact: [],
    ranking_score: Number.isFinite(bestRank) ? bestRank : rep.ranking_score,
    ranking_score_scaled: bestScaledRank,
    max_pheno_score: bestPheno || rep.max_pheno_score,
    cnv_sv_sort_score: Number.isFinite(bestCombined) ? bestCombined : rep.cnv_sv_sort_score,
    genes_total: genes.length, merged_segment_ids: segments.map(v => v.id),
    is_merged_parent: true,
  };
}

let _cnvSvMergeViewCache = null;

function _cnvSvMergeView() {
  const cnvVariants = state.data?.cnv_variants || null;
  const svVariants = state.data?.sv_variants || null;
  const savedMerges = state.reports?.cnv_sv_merges || null;
  if (_cnvSvMergeViewCache
      && _cnvSvMergeViewCache.cnvVariants === cnvVariants
      && _cnvSvMergeViewCache.svVariants === svVariants
      && _cnvSvMergeViewCache.savedMerges === savedMerges) {
    return _cnvSvMergeViewCache;
  }

  const merges = [];
  const consumed = new Set();
  [..._automaticCnvSvMerges(), ..._confirmedCnvSvMerges()].forEach(merge => {
    if ((merge.member_ids || []).some(id => consumed.has(id))) return;
    merges.push(merge);
    (merge.member_ids || []).forEach(id => consumed.add(id));
  });
  const parents = {};
  merges.forEach(merge => {
    const parent = _cnvSvBuildParent(merge);
    if (parent) parents[parent.id] = parent;
  });
  _cnvSvMergeViewCache = { cnvVariants, svVariants, savedMerges, merges, parents };
  return _cnvSvMergeViewCache;
}

function _cnvSvVirtualParents() {
  return _cnvSvMergeView().parents;
}

function _cnvSvIdsForTier(tier) {
  const cats = tier.startsWith("CNV-")
    ? state.data?.cnv_categories
    : state.data?.sv_categories;
  const ids = [...((cats && cats[tier]) || [])];
  const replacements = new Map();
  const suppressed = new Set();
  _effectiveCnvSvMerges().forEach(merge => {
    const tierMembers = (merge.member_ids || []).filter(id => ids.includes(id));
    if (!tierMembers.length) return;
    const parent = _cnvSvBuildParent(merge);
    if (!parent) return;
    tierMembers.forEach(id => suppressed.add(id));
    const anchor = tierMembers.reduce((best, id) => ids.indexOf(id) < ids.indexOf(best) ? id : best);
    replacements.set(anchor, parent.id);
  });
  const out = [];
  ids.forEach(id => {
    if (replacements.has(id)) out.push(replacements.get(id));
    if (!suppressed.has(id)) out.push(id);
  });
  return out.sort((a, b) => {
    const diff = _annotSvSortScore(_cnvSvVariantById(b)) - _annotSvSortScore(_cnvSvVariantById(a));
    return diff || String(a).localeCompare(String(b));
  });
}

// ---------- CNV/SV near-duplicate clustering ----------------------
//
// SV callers (Manta / DELLY / GRIDSS / cnvkit) routinely report the
// same biological event as several nearly-identical SVs whose
// breakpoints differ by a few hundred bp. This is a pure-visual
// display collapse: pick a representative per cluster (the first
// id, which carries the highest ranking_score because tiers are
// pre-sorted) and tuck the others inside an expandable "顯示同位點
// N 個近似 SV" detail block on the rep card. No edits / no data are
// lost — the alternatives are still full, interactive cards.
//
// Long-term fix is upstream: pipeline should run Truvari / SURVIVOR
// before AnnotSV so AnnotSV sees one SV per real event. When that
// ships this clustering naturally becomes a no-op.

const CNV_SV_CLUSTER_OVERLAP_THRESHOLD = 0.8;
const CNV_SV_MERGE_GAP_THRESHOLD = 250000;

function _cnvSvSpan(v) {
  const s = Number(v?.POS);
  const e = Number(v?.END);
  if (!Number.isFinite(s) || !Number.isFinite(e) || e <= s) return null;
  return [s, e];
}

function _cnvSvReciprocalOverlap(a, b) {
  const sa = _cnvSvSpan(a);
  const sb = _cnvSvSpan(b);
  if (!sa || !sb) return 0;
  const ov = Math.max(0, Math.min(sa[1], sb[1]) - Math.max(sa[0], sb[0]));
  return ov / Math.max(sa[1] - sa[0], sb[1] - sb[0]);
}

function _cnvSvClusterIds(ids) {
  // {reps: [repId, ...], members: {repId: [otherIds]}} preserving input
  // order so the first id of each cluster (= highest ranking_score
  // because tiers are pre-sorted) becomes the representative.
  const reps = [];
  const members = {};
  for (const id of ids) {
    const v = _cnvSvVariantById(id);
    if (!v) continue;
    let foundRep = null;
    for (const repId of reps) {
      const rv = _cnvSvVariantById(repId);
      if (!rv) continue;
      if ((rv.CHROM || "") !== (v.CHROM || "")) continue;
      if ((rv.sv_type || "") !== (v.sv_type || "")) continue;
      if (_cnvSvReciprocalOverlap(rv, v) >= CNV_SV_CLUSTER_OVERLAP_THRESHOLD) {
        foundRep = repId;
        break;
      }
    }
    if (foundRep) {
      members[foundRep].push(id);
    } else {
      reps.push(id);
      members[id] = [];
    }
  }
  return { reps, members };
}

function _cnvSvCompatibleSegments(a, b) {
  if (!a || !b) return false;
  if (!["DEL", "DUP"].includes(String(a.sv_type || "").toUpperCase())) return false;
  if ((a.source || "") !== (b.source || "")) return false;
  if ((a.CHROM || "") !== (b.CHROM || "")) return false;
  return String(a.sv_type || "").toUpperCase() === String(b.sv_type || "").toUpperCase();
}

function _cnvSvAdjacentMergeGroups(ids) {
  const sorted = ids.map(_cnvSvBaseVariantById).filter(Boolean)
    .sort((a, b) => String(a.source).localeCompare(String(b.source))
      || String(a.CHROM).localeCompare(String(b.CHROM))
      || String(a.sv_type).localeCompare(String(b.sv_type))
      || Number(a.POS) - Number(b.POS));
  const groups = [];
  let current = [];
  sorted.forEach(v => {
    const prev = current[current.length - 1];
    const gap = prev ? Number(v.POS) - Number(prev.END) : Infinity;
    if (prev && _cnvSvCompatibleSegments(prev, v) && gap >= 0 && gap <= CNV_SV_MERGE_GAP_THRESHOLD) {
      current.push(v);
    } else {
      if (current.length >= 2) groups.push(current);
      current = [v];
    }
  });
  if (current.length >= 2) groups.push(current);
  return groups;
}

function _automaticCnvSvMerges() {
  const variants = [
    ...Object.values(state.data?.cnv_variants || {}),
    ...Object.values(state.data?.sv_variants || {}),
  ];
  return _cnvSvAdjacentMergeGroups(variants.map(v => v.id)).map(segments => {
    const first = segments[0], last = segments[segments.length - 1];
    const source = first.source || "cnv";
    const chrom = first.CHROM || "";
    const start = Number(first.POS), end = Number(last.END);
    const svType = first.sv_type || "";
    const id = _cnvSvMergeId(source, chrom, start, end, svType);
    const memberIds = segments.map(v => v.id);
    return { id, source, member_ids: memberIds };
  });
}

function _effectiveCnvSvMerges() {
  return _cnvSvMergeView().merges;
}

function renderCnvSvTabBar() {
  const bar = document.getElementById("cnv-sv-tab-bar");
  if (!bar) return;
  if (!CNV_SV_TIER_ORDER.includes(activeCnvSvTab)) activeCnvSvTab = null;
  if (!activeCnvSvTab) activeCnvSvTab = "CNV-1A";

  bar.innerHTML = CNV_SV_TIER_ORDER.map(t => {
    const active = t === activeCnvSvTab ? " active" : "";
    const ids = _cnvSvIdsForTier(t);
    const loading = t.startsWith("CNV-") ? !!state.cnvPending : !!state.svPending;
    return `<button type="button" class="tier-tab ${CNV_SV_TIER_CLASS[t]}${active}" data-tier="${t}">
              <span class="tier-tab-title">${escapeHtml(CNV_SV_TITLES[t])}</span>
              <span class="tier-tab-count">${loading ? "…" : "Total " + ids.length}</span>
            </button>`;
  }).join("");

  // Each tier panel wraps its cards in a .block-body so the SNV
  // tier-panel padding rule (`.tier-panel > .block-body { padding-top: 8px }`)
  // gives the same inset-card look. Without the wrapper the cards
  // would butt right up against the colored panel edge.
  CNV_SV_TIER_ORDER.forEach(tier => {
    const panel = document.querySelector(`#cnv-sv-tab-panels .tier-panel[data-tier="${tier}"]`);
    if (!panel) return;
    const ids = _cnvSvIdsForTier(tier);
    const loading = tier.startsWith("CNV-") ? !!state.cnvPending : !!state.svPending;
    const isClinical = tier.endsWith("-1A") || tier.endsWith("-2A");
    panel.innerHTML = "";
    if (loading) {
      const wrap = document.createElement("div");
      wrap.className = "block-body";
      wrap.innerHTML = `<div class="analysis-card-empty">載入中…</div>`;
      panel.appendChild(wrap);
      return;
    }
    if (!ids.length) {
      const empty = document.createElement("div");
      empty.className = "block-body";
      empty.innerHTML = (isClinical && !state.data?.has_phenotype)
        ? `<div class="analysis-card-empty">請先設定 phenotype（HPO / panel），才會有 Clinical 結果。</div>`
        : `<div class="analysis-card-empty">（無資料）</div>`;
      panel.appendChild(empty);
      return;
    }
    const body = document.createElement("div");
    body.className = "block-body open";
    const { reps, members } = _cnvSvClusterIds(ids);
    reps.forEach((repId, i) => {
      const v = _cnvSvVariantById(repId);
      if (!v) return;
      body.appendChild(renderCnvSvCard(v, repId, { tier, index: i + 1 }));
      if (v.is_merged_parent) {
        const det = document.createElement("details");
        det.className = "cnv-sv-cluster-alts";
        det.innerHTML = `<summary>顯示自動整合前 ${v.merged_segment_ids.length} 個原始片段</summary>`;
        const altBody = document.createElement("div");
        altBody.className = "cnv-sv-cluster-alt-body";
        v.merged_segment_ids.forEach(mid => {
          const mv = _cnvSvBaseVariantById(mid);
          if (!mv) return;
          const altCard = renderCnvSvCard(mv, mid, { tier });
          altCard.classList.add("cnv-sv-cluster-alt");
          altBody.appendChild(altCard);
        });
        det.appendChild(altBody);
        body.appendChild(det);
      }
      const memberIds = members[repId] || [];
      if (memberIds.length) {
        const det = document.createElement("details");
        det.className = "cnv-sv-cluster-alts";
        const summary = document.createElement("summary");
        summary.textContent = `⤴ 顯示同位點 ${memberIds.length} 個近似 SV（不同 caller / breakpoint）`;
        det.appendChild(summary);
        const altBody = document.createElement("div");
        altBody.className = "cnv-sv-cluster-alt-body";
        memberIds.forEach(mid => {
          const mv = _cnvSvVariantById(mid);
          if (!mv) return;
          const altCard = renderCnvSvCard(mv, mid, { tier });
          altCard.classList.add("cnv-sv-cluster-alt");
          altBody.appendChild(altCard);
        });
        det.appendChild(altBody);
        body.appendChild(det);
      }
    });
    panel.appendChild(body);
  });
  applyCnvSvTabActive();
}

function applyCnvSvTabActive() {
  const panels = document.getElementById("cnv-sv-tab-panels");
  if (!panels) return;
  panels.querySelectorAll(".tier-panel").forEach(p => {
    p.classList.toggle("active", p.dataset.tier === activeCnvSvTab);
  });
}

// ---------- Mitochondria tier tabs --------------------------------
// MITO-1 = ClinVar P/LP or MITOMAP confirmed/pathogenic, MITO-2 =
// rare/not-observed or reported non-benign mtDNA variants, MITO-3 =
// other PASS mtDNA variants. MITOMAP is used for tiering but not shown
// on the card.
const MITO_TIER_ORDER = ["MITO-1", "MITO-2", "MITO-3"];
const MITO_TITLES = {
  "MITO-1": "1 Pathogenic",
  "MITO-2": "2 Rare / reported mtDNA variant",
  "MITO-3": "3 Other variant",
};
let activeMitoTab = null;

function _mitoIdsForTier(tier) {
  return (state.data?.mito_categories && state.data.mito_categories[tier]) || [];
}

function renderMitoTabBar() {
  const bar = document.getElementById("mito-tab-bar");
  if (!bar) return;
  if (!MITO_TIER_ORDER.includes(activeMitoTab)) activeMitoTab = null;
  if (!activeMitoTab) activeMitoTab = "MITO-1";

  const loading = !!state.mitoPending;
  bar.innerHTML = MITO_TIER_ORDER.map(t => {
    const active = t === activeMitoTab ? " active" : "";
    const ids = _mitoIdsForTier(t);
    return `<button type="button" class="tier-tab tier-mito${active}" data-tier="${t}">
              <span class="tier-tab-title">${escapeHtml(MITO_TITLES[t])}</span>
              <span class="tier-tab-count">${loading ? "…" : "Total " + ids.length}</span>
            </button>`;
  }).join("");

  MITO_TIER_ORDER.forEach(tier => {
    const panel = document.querySelector(`#mito-tab-panels .tier-panel[data-tier="${tier}"]`);
    if (!panel) return;
    const ids = _mitoIdsForTier(tier);
    panel.innerHTML = "";
    if (loading) {
      const wrap = document.createElement("div");
      wrap.className = "block-body";
      wrap.innerHTML = `<div class="analysis-card-empty">載入中…</div>`;
      panel.appendChild(wrap);
      return;
    }
    if (!ids.length) {
      const empty = document.createElement("div");
      empty.className = "block-body";
      empty.innerHTML = `<div class="analysis-card-empty">${
        tier === "MITO-1" ? "（無 pathogenic mtDNA 變異）"
        : tier === "MITO-2" ? "（無 rare / reported mtDNA 變異）"
        : "（無其他 mtDNA 變異）"
      }</div>`;
      panel.appendChild(empty);
      return;
    }
    const body = document.createElement("div");
    body.className = "block-body open";
    ids.forEach((id, i) => {
      const v = state.data?.mito_variants?.[id];
      if (!v) return;
      body.appendChild(renderMitoCard(v, id, { index: i + 1 }));
    });
    panel.appendChild(body);
  });
  applyMitoTabActive();
}

function applyMitoTabActive() {
  const panels = document.getElementById("mito-tab-panels");
  if (!panels) return;
  panels.querySelectorAll(".tier-panel").forEach(p => {
    p.classList.toggle("active", p.dataset.tier === activeMitoTab);
  });
}

// ---------- STR tier tabs -----------------------------------------
const STR_TIER_ORDER = ["STR-P", "STR-I", "STR-N"];
const STR_TITLES = {
  "STR-P": "Pathogenic",
  "STR-I": "Intermediate / Borderline",
  "STR-N": "Normal / No threshold",
};
let activeStrTab = null;

function _strIdsForTier(tier) {
  return (state.data?.str_categories && state.data.str_categories[tier]) || [];
}

function renderStrTabBar() {
  const bar = document.getElementById("str-tab-bar");
  if (!bar) return;
  if (!STR_TIER_ORDER.includes(activeStrTab)) activeStrTab = null;
  if (!activeStrTab) activeStrTab = "STR-P";
  const loading = !!state.strPending;
  bar.innerHTML = STR_TIER_ORDER.map(t => {
    const active = t === activeStrTab ? " active" : "";
    const ids = _strIdsForTier(t);
    return `<button type="button" class="tier-tab tier-str${active}" data-tier="${t}">
              <span class="tier-tab-title">${escapeHtml(STR_TITLES[t])}</span>
              <span class="tier-tab-count">${loading ? "…" : "Total " + ids.length}</span>
            </button>`;
  }).join("");

  STR_TIER_ORDER.forEach(tier => {
    const panel = document.querySelector(`#str-tab-panels .tier-panel[data-tier="${tier}"]`);
    if (!panel) return;
    panel.innerHTML = "";
    if (loading) {
      const wrap = document.createElement("div");
      wrap.className = "block-body";
      wrap.innerHTML = `<div class="analysis-card-empty">載入中…</div>`;
      panel.appendChild(wrap);
      return;
    }
    const ids = _strIdsForTier(tier);
    if (!ids.length) {
      const empty = document.createElement("div");
      empty.className = "block-body";
      empty.innerHTML = `<div class="analysis-card-empty">${
        tier === "STR-P" ? "（無 pathogenic STR）"
        : tier === "STR-I" ? "（無 intermediate / borderline STR）"
        : "（無 normal / no-threshold STR 資料）"
      }</div>`;
      panel.appendChild(empty);
      return;
    }
    const body = document.createElement("div");
    body.className = "block-body open";
    ids.forEach((id, i) => {
      const v = state.data?.str_variants?.[id];
      if (!v) return;
      body.appendChild(renderStrCard(v, id, { index: i + 1 }));
    });
    panel.appendChild(body);
  });
  applyStrTabActive();
}

function applyStrTabActive() {
  const panels = document.getElementById("str-tab-panels");
  if (!panels) return;
  panels.querySelectorAll(".tier-panel").forEach(p => {
    p.classList.toggle("active", p.dataset.tier === activeStrTab);
  });
}

function _strClassificationClass(v) {
  const c = String(v?.CLASSIFICATION || "").toLowerCase();
  if (c === "pathogenic") return "str-class-pathogenic";
  if (c === "intermediate" || c === "borderline") return "str-class-intermediate";
  if (c === "no_threshold") return "str-class-no-threshold";
  return "str-class-normal";
}

function _strThresholdParts(v) {
  const parts = [];
  if (v.BENIGN_MIN || v.BENIGN_MAX) {
    parts.push(`Benign ${escapeHtml(v.BENIGN_MIN || "—")}–${escapeHtml(v.BENIGN_MAX || "—")}`);
  }
  if (v.INTERMEDIATE_MIN || v.INTERMEDIATE_MAX) {
    parts.push(`Intermediate ${escapeHtml(v.INTERMEDIATE_MIN || "—")}–${escapeHtml(v.INTERMEDIATE_MAX || "—")}`);
  }
  if (v.PATHOGENIC_MIN || v.PATHOGENIC_MAX) {
    parts.push(`Pathogenic ${escapeHtml(v.PATHOGENIC_MIN || "—")}–${escapeHtml(v.PATHOGENIC_MAX || "—")}`);
  }
  return parts.join(" · ");
}

function renderStrCard(v, id, opts = {}) {
  const card = document.createElement("div");
  card.className = `variant-card str-card ${_strClassificationClass(v)}`;
  card.dataset.variantId = id;
  const coord = [v.CHROM, v.POS && `${v.POS}-${v.END || v.POS}`].filter(Boolean).join(":");
  const threshold = _strThresholdParts(v) || "Threshold not available";
  card.innerHTML = `
    <div class="variant-card-header str-card-header">
      <div class="variant-title">
        ${opts.index ? `<span class="rank-badge">${opts.index}</span>` : ""}
        <span class="gene">${escapeHtml(v.GENE || "STR")}</span>
        <span class="muted">${escapeHtml(v.STR_ID || "")}</span>
      </div>
      <span class="str-class-badge">${escapeHtml(v.CLASSIFICATION || "normal")}</span>
    </div>
    <div class="str-subtitle">${escapeHtml(v.DISEASE || "—")} · ${escapeHtml(v.INHERITANCE || "—")} · ${escapeHtml(v.TYPE || "—")} · ${escapeHtml(v.MOTIF || "—")}</div>
    <div class="str-repeat-grid">
      <div><span class="label">Allele 1</span><strong>${escapeHtml(v.REPCN_A1 || "—")}</strong></div>
      <div><span class="label">Allele 2</span><strong>${escapeHtml(v.REPCN_A2 || "—")}</strong></div>
      <div><span class="label">CI</span><strong>${escapeHtml(v.REPCI || "—")}</strong></div>
      <div><span class="label">DP</span><strong>${escapeHtml(v.DP || "—")}</strong></div>
    </div>
    <div class="str-threshold">${threshold}</div>
    <details class="str-detail">
      <summary>詳細</summary>
      <div class="str-detail-grid">
        <div><span class="label">Locus structure</span><span>${escapeHtml(v.LOCUS_STRUCTURE || "—")}</span></div>
        <div><span class="label">Coordinate</span><span>${escapeHtml(coord || "—")}</span></div>
        <div><span class="label">Pipeline</span><span>${escapeHtml(v.PIPELINE || "—")}</span></div>
      </div>
    </details>`;
  return card;
}

// Click dispatch for all tier-tab groups (SNV, CNV/SV, Mito, STR).
// The tab's data-tier tells us which group, so the active-class toggle
// stays scoped to that group's bar.
document.addEventListener("click", ev => {
  const tab = ev.target.closest(".tier-tab");
  if (!tab) return;
  const tier = tab.dataset.tier;
  if (!tier) return;
  let barId, current, setActive, applyActive;
  if (CNV_SV_TIER_ORDER.includes(tier)) {
    barId = "cnv-sv-tab-bar"; current = activeCnvSvTab;
    setActive = t => { activeCnvSvTab = t; }; applyActive = applyCnvSvTabActive;
  } else if (MITO_TIER_ORDER.includes(tier)) {
    barId = "mito-tab-bar"; current = activeMitoTab;
    setActive = t => { activeMitoTab = t; }; applyActive = applyMitoTabActive;
  } else if (STR_TIER_ORDER.includes(tier)) {
    barId = "str-tab-bar"; current = activeStrTab;
    setActive = t => { activeStrTab = t; }; applyActive = applyStrTabActive;
  } else {
    barId = "tier-tab-bar"; current = activeTierTab;
    setActive = t => { activeTierTab = t; };
    applyActive = renderActiveSnvTier;
  }
  if (tier === current) return;
  setActive(tier);
  document.querySelectorAll(`#${barId} .tier-tab`).forEach(b => {
    b.classList.toggle("active", b.dataset.tier === tier);
  });
  applyActive();
});

// ---------- Mitochondria variant card -----------------------------
const MITO_LOCUS_LABELS = {
  protein: "protein-coding", tRNA: "tRNA", rRNA: "rRNA",
  control: "control region", intergenic: "intergenic", unknown: "—",
};

function _mitoExternalLinks(v) {
  // gnomAD v3 has a dedicated mtDNA dataset; M-<pos>-<ref>-<alt> works
  // for SNVs (gnomAD normalises).
  const links = [];
  if (v.REF && v.ALT) {
    links.push({
      label: "gnomAD-MT",
      href: `https://gnomad.broadinstitute.org/variant/M-${v.POS}-${escapeAttr(v.REF)}-${escapeAttr(v.ALT)}?dataset=gnomad_r3`,
    });
  }
  if (v.clinvar_variation_id) {
    links.push({
      label: "ClinVar",
      href: `https://www.ncbi.nlm.nih.gov/clinvar/variation/${encodeURIComponent(v.clinvar_variation_id)}/`,
    });
  }
  return links;
}

function _mitoReviewerStatusSel(id) {
  const status = (state.reports?.status?.[id]) || "";
  return _renderStatusRadio(id, status, statusOptions("candidate"));
}

function _mitoHeteroplasmy(v) {
  const h = v.heteroplasmy;
  if (h == null || !Number.isFinite(Number(h))) return "—";
  return `${(Number(h) * 100).toFixed(1)}%`;
}

function _mitoConsequenceLabel(v) {
  return (v.consequence || "—");
}

function _formatMitoAf(v) {
  const parts = [];
  const add = (label, val) => {
    if (val == null || val === "" || val === ".") return;
    const n = Number(val);
    parts.push(Number.isFinite(n) ? `${label} ${(n * 100).toPrecision(3)}%` : `${label} ${val}`);
  };
  add("hom", v.gnomad_mito_af_hom);
  add("het", v.gnomad_mito_af_het);
  add("AF", v.gnomad_mito_af);
  if (v.gnomad_mito_an) parts.push(`AN ${v.gnomad_mito_an}`);
  return parts.join(" · ") || "—";
}

function _mitoClinvarDiseaseList(v, id) {
  const diseases = Array.isArray(v.clinvar_diseases) ? v.clinvar_diseases : [];
  if (!diseases.length) return "";
  const picked = getEdit(id, "report_diseases_clinvar") || {};
  const rows = diseases.map((d, i) => {
    const key = String(i);
    const checked = picked[key] ? "checked" : "";
    return `<label class="disease-row mito-clinvar-disease">
      <input type="checkbox" class="mito-disease-pick" data-id="${escapeAttr(id)}" data-idx="${escapeAttr(key)}" ${checked} title="報告要發這個 ClinVar disease" />
      <span class="disease-summary-text">${escapeHtml(d)}</span>
    </label>`;
  });
  return `<div class="cnv-sv-section">
    <div class="cnv-sv-section-title">ClinVar disease</div>
    <div class="disease-list">${rows.join("")}</div>
  </div>`;
}

// Mutect2-mito FILTER flag → plain-Chinese gloss (shown as a tooltip).
const MITO_FILTER_GLOSS = {
  PASS:             "通過所有過濾",
  weak_evidence:    "變異訊號弱（likelihood 未達門檻）— 常見於低 heteroplasmy 雜訊位點",
  base_qual:        "alt allele 的中位 base quality 偏低",
  blacklisted_site: "落在 mtDNA 已知問題區（poly-C tract、NUMT 高相似區等黑名單）",
  possible_numt:    "疑似來自核基因組的 mtDNA 偽基因片段（NUMT）",
  contamination:    "疑似樣本污染",
  strand_bias:      "證據只來自單一 read 方向",
  strict_strand:    "alt allele 在兩個 read 方向都沒被代表到",
  slippage:         "STR 區域的 polymerase slippage",
  map_qual:         "mapping quality 異常（ref 與 alt 差異大）",
  position:         "alt 變異離 read 末端太近",
  clustered_events: "附近 somatic events 過多",
  haplotype:        "靠近同一 haplotype 上被過濾掉的變異",
  multiallelic:     "此位點 alt allele 過多",
  fragment:         "ref/alt 的中位 fragment length 差異過大",
};
function _mitoFilterTitle(filt) {
  const parts = (filt || "").split(";").map(s => s.trim()).filter(Boolean);
  if (!parts.length) return MITO_FILTER_GLOSS.PASS;
  return parts.map(p => `${p}：${MITO_FILTER_GLOSS[p] || "（未知旗標）"}`).join("\n");
}
const _MITO_TLOD_TITLE = "Mutect2 tumor LOD：log10(變異存在 / 不存在) 的 likelihood ratio。越高 = 越確定是真變異（非測序錯誤）；一般 >6 算可靠，1-2 多為雜訊。";

function _renderMitoDetailBox(v, id) {
  const filt = v.filter && v.filter !== "." ? v.filter : "";
  const tlod = (v.TLOD != null) ? Number(v.TLOD).toFixed(2) : "—";
  const ad   = v.AD || "—";
  const dp   = (v.depth != null) ? v.depth : "—";
  const consL = _mitoConsequenceLabel(v);
  // Mito ACMG override — the mito pipeline does not provide an ACMG class, so
  // the reviewer picks one manually here. Persisted to
  // state.reports.edits[id].ACMG_classification_mito and consumed by
  // the docx report (`_acmg_label` → mito column + note-2 wording).
  const mitoAcmg = (getEdit(id, "ACMG_classification_mito") || "").trim();
  const mitoAcmgSig = classifySignificance(mitoAcmg) || "";
  const acmgOpts = ["", "Pathogenic", "Likely pathogenic", "Uncertain significance",
                    "Likely benign", "Benign"];
  const acmgSelect = `
    <select class="mito-acmg-select ${mitoAcmgSig}" data-id="${escapeAttr(id)}">
      ${acmgOpts.map(o => `<option value="${escapeAttr(o)}" ${o===mitoAcmg?"selected":""}>${o || "—"}</option>`).join("")}
    </select>`;
  return `<div class="cnv-sv-detail-box">
    <div class="cnv-sv-detail-row">
      <span><strong>ACMG:</strong> ${acmgSelect}</span>
      <span><strong>變化:</strong> ${escapeHtml(v.REF || "?")}→${escapeHtml(v.ALT || "?")}</span>
      <span><strong>類型:</strong> ${escapeHtml(MITO_LOCUS_LABELS[v.locus_type] || v.locus_type || "—")}</span>
      <span><strong>Heteroplasmy:</strong> ${_mitoHeteroplasmy(v)} <span class="muted">(AD ${escapeHtml(ad)} · DP ${dp})</span></span>
      ${v.genotype ? `<span><strong>GT:</strong> ${escapeHtml(v.genotype)}</span>` : ""}
      ${filt ? `<span data-tip="${escapeAttr(_mitoFilterTitle(filt))}"><strong>Filter:</strong> ${escapeHtml(filt)} <span class="muted" style="cursor:help">ⓘ</span></span>` : ""}
    </div>
    <div class="cnv-sv-detail-row">
      <span><strong>Consequence:</strong> ${escapeHtml(consL)}</span>
      ${v.impact ? `<span><strong>Impact:</strong> ${escapeHtml(v.impact)}</span>` : ""}
      ${v.biotype ? `<span><strong>Biotype:</strong> ${escapeHtml(v.biotype)}</span>` : ""}
      ${v.aa_change ? `<span><strong>Protein change:</strong> ${escapeHtml(v.aa_change)}</span>` : ""}
      ${(() => {
        const sig = (v.CLNSIG || "").trim();
        if (!sig) return `<span><strong>ClinVar:</strong> —</span>`;
        const cls = classifySignificance(sig) || "";
        const stars = v.clinvar_stars != null && v.clinvar_stars !== "" && Number(v.clinvar_stars) > 0
          ? ` ${"★".repeat(Number(v.clinvar_stars))}` : "";
        return `<span><strong>ClinVar:</strong> <span class="acmg-class ${cls}">${escapeHtml(sig.replace(/_/g," "))}${escapeHtml(stars)}</span></span>`;
      })()}
      <span><strong>gnomAD mito:</strong> ${escapeHtml(_formatMitoAf(v))}</span>
      ${v.TLOD != null ? `<span data-tip="${escapeAttr(_MITO_TLOD_TITLE)}"><strong>TLOD:</strong> ${tlod} <span class="muted" style="cursor:help">ⓘ</span></span>` : ""}
    </div>
  </div>`;
}

function renderMitoCard(v, id, opts = {}) {
  const card = document.createElement("div");
  card.className = "variant-card mito-card";
  card.dataset.id = id;
  const idxTxt = opts.index ? `<span class="card-idx">#${opts.index}</span>` : "";
  const locusCls = `mito-locus-${(v.locus_type || "unknown")}`;
  const locusPill = `<span class="mito-locus-pill ${locusCls}">${escapeHtml(MITO_LOCUS_LABELS[v.locus_type] || v.locus_type || "—")}</span>`;
  const hgvs = v.HGVS_M || `m.${v.POS}${v.REF}>${v.ALT}`;
  const links = _mitoExternalLinks(v).map(l =>
    `<a href="${escapeAttr(l.href)}" target="_blank" rel="noopener">${escapeHtml(l.label)}</a>`
  ).join("");
  const igvLink = `<a href="#" class="btn-igv" data-id="${escapeAttr(id)}" title="在 IGV 內檢視">IGV</a>`;
  const comment = (getEdit(id, "comment") || "");
  card.innerHTML = `
    <div class="variant-head">
      ${idxTxt}
      ${_mitoReviewerStatusSel(id)}
      ${locusPill}
      <span class="cnv-sv-pos">${v.gene_symbol ? `${escapeHtml(v.gene_symbol)} ` : ""}${escapeHtml(hgvs)}<button class="btn-copy" data-copy="${escapeAttr((v.gene_symbol ? v.gene_symbol + " " : "") + hgvs)}" title="複製">${COPY_ICON_SVG}</button>
      </span>
      <span class="mito-het-badge" title="heteroplasmy fraction">${_mitoHeteroplasmy(v)}</span>
      <span style="flex:1"></span>
      <span class="ext-links">${igvLink}${links}</span>
    </div>
    ${_renderMitoDetailBox(v, id)}
    ${_mitoClinvarDiseaseList(v, id)}
    <div class="cnv-sv-section cnv-sv-comment">
      <div class="cnv-sv-section-title">Comment</div>
      <textarea class="cnv-sv-comment-text" data-id="${escapeAttr(id)}" rows="2" placeholder="備註">${escapeHtml(comment)}</textarea>
    </div>
  `;
  return card;
}

// ---------- CNV / SV variant card rendering ----------------------
//
// One card per AnnotSV record. Layout mirrors the SNV variant card's
// visual language (header pills + collapsible body + status dropdown
// + comment + disease list) but the fields are SV-specific:
//   • position / type / cytoband / length / copy-number
//   • AnnotSV's own ACMG class (1-5) and ranking score
//   • per-gene table built from split rows (Tx / Location / OMIM)
//   • pathogenic-region overlap (P_loss / P_gain) + benign AF
//   • AnnotSV reasoning text (collapsed)
//   • disease list synthesised from gene OMIM_phenotype lines
//
// Edits (status / comment / report-disease checkbox) reuse the same
// state.reports.{status,edits} dicts as SNV cards — AnnotSV_IDs and
// SNV chr-pos-ref-alt ids never collide so one flat namespace works.

// ACMG_class 1..5 → human label (drives the per-card dropdown). The
// dropdown writes the integer back; the label only shows in the UI.
const SV_ACMG_LABELS = {
  1: "Benign",
  2: "Likely benign",
  3: "VUS",
  4: "Likely pathogenic",
  5: "Pathogenic",
};
// Mirror of SNV's classifySignificance(): map AnnotSV's 1..5 numeric
// scale onto the SNV sig-* colour classes so the dropdown reads the
// same way visually.
const SV_ACMG_SIG_CLASS = {
  5: "sig-p",
  4: "sig-lp",
  3: "sig-vus",
  2: "sig-lb",
  1: "sig-b",
};

function _fmtPos(n) {
  if (n == null) return "?";
  return Number(n).toLocaleString();
}

function _fmtBp(n) {
  if (n == null) return "";
  const a = Math.abs(Number(n));
  if (a >= 1e6) return `${(a / 1e6).toFixed(2)}Mb`;
  if (a >= 1e3) return `${(a / 1e3).toFixed(1)}kb`;
  return `${a}bp`;
}

function _normalizeChrom(c) {
  if (!c) return "";
  return String(c).startsWith("chr") ? c : `chr${c}`;
}

function _chromNumber(c) {
  // "chr12" → "12"; "12" → "12"; used for cytoband prefix.
  return String(c || "").replace(/^chr/, "");
}

function _cnvSvExternalLinks(v) {
  const build = state.data?.genome_build === "hg19" ? "hg19" : "hg38";
  const chrom = _normalizeChrom(v.CHROM);
  const region = `${chrom}:${v.POS}-${v.END}`;
  const ucscDb = build === "hg19" ? "hg19" : "hg38";
  const links = [
    { label: "UCSC",     href: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=${ucscDb}&position=${region}` },
    { label: "DECIPHER", href: `https://www.deciphergenomics.org/search/patients/results?q=${encodeURIComponent(region)}` },
    { label: "dbVar",    href: `https://www.ncbi.nlm.nih.gov/dbvar/?term=${encodeURIComponent(region)}` },
  ];
  if (v.gene_symbol) {
    links.push({ label: "GeneCards", href: `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${encodeURIComponent(v.gene_symbol)}` });
  }
  return links;
}

function _cnvSvAcmgClassValue(id, v) {
  // Reviewer-edited override (numeric 1-5) takes precedence; falls
  // back to AnnotSV's own ACMG_class.
  const edited = getEdit(id, "ACMG_class_sv");
  if (edited != null && edited !== "") return Number(edited);
  return (v.acmg_class != null) ? Number(v.acmg_class) : null;
}

function _renderCnvSvHeader(v, id, opts) {
  const sourceLabel = v.source === "cnv" ? "CNV" : "SV";
  const typeChip = `<span class="sv-type-pill sv-type-${escapeAttr(v.sv_type || "")}">${escapeHtml(v.sv_type || "?")}</span>`;
  const chrom = _normalizeChrom(v.CHROM);
  const cytoBoth = v.cytoband ? `${_chromNumber(v.CHROM)}${v.cytoband}` : "";
  const lengthPart = v.length != null ? _fmtBp(v.length) : "";
  const region = `${chrom}:${_fmtPos(v.POS)}-${_fmtPos(v.END)}`;
  const regionRaw = `${chrom}:${v.POS}-${v.END}`;
  // SNV-card-style status dropdown (1/2/C/0/X). Reuses the same
  // state.reports.status dict + the same options as SNV — picking C
  // routes the variant into the Candidate variants report section.
  const status = (state.reports?.status?.[id]) || "";
  const statusSel = _renderStatusRadio(id, status, statusOptions("candidate"));
  const idxTxt = opts.index ? `<span class="card-idx">#${opts.index}</span>` : "";
  const extLinks = _cnvSvExternalLinks(v).map(l =>
    `<a href="${escapeAttr(l.href)}" target="_blank" rel="noopener">${escapeHtml(l.label)}</a>`
  ).join("");
  const igvLink = `<a href="#" class="btn-igv" data-id="${escapeAttr(id)}" title="在 IGV 內檢視">IGV</a>`;

  return `<div class="variant-head">
    ${idxTxt}
    ${statusSel}
    <span class="cnv-sv-source-tag">${sourceLabel}</span>
    ${typeChip}
    <span class="cnv-sv-pos">${escapeHtml(region)}<button class="btn-copy" data-copy="${escapeAttr(regionRaw)}" title="複製座標">${COPY_ICON_SVG}</button>
      ${lengthPart ? ` <span class="muted" style="font-size:11px">· ${escapeHtml(lengthPart)}</span>` : ""}
      ${cytoBoth ? ` <span class="muted" style="font-size:11px">· ${escapeHtml(cytoBoth)}</span>` : ""}
    </span>
    <span class="ext-links">${igvLink}${extLinks}</span>
  </div>`;
}

function _renderCnvSvDetailBox(v, id) {
  const cn = (v.copy_number != null) ? ` · CN ${v.copy_number}` : "";
  const filter = v.filter && v.filter !== "." ? v.filter : "PASS";
  const qual = (v.qual != null) ? Number(v.qual).toFixed(2) : "—";
  const zyg = v.zygosity || "—";
  // ACMG dropdown borrows SNV's sig-* colour scale (sig-p…sig-b) so
  // the field colour matches the rest of the app — Pathogenic red,
  // VUS yellow, etc. AnnotSV's numeric class is the default; the
  // reviewer's override lives on state.reports.edits[id].ACMG_class_sv
  // (separate field from SNV's `ACMG_classification` so they don't
  // collide).
  const acmgVal = _cnvSvAcmgClassValue(id, v);
  const sigClass = SV_ACMG_SIG_CLASS[acmgVal] || "";
  const acmgSelect = `
    <select class="cnv-sv-acmg-select ${sigClass}" data-id="${escapeAttr(id)}">
      <option value="" ${acmgVal==null ? "selected" : ""}>—</option>
      ${[5,4,3,2,1].map(n =>
        `<option value="${n}" ${acmgVal===n?"selected":""}>${escapeHtml(SV_ACMG_LABELS[n])}</option>`
      ).join("")}
    </select>`;
  const score = (v.ranking_score != null) ? Number(v.ranking_score).toFixed(2) : "—";
  const reasoning = v.ranking_criteria
    ? (() => {
        const items = v.ranking_criteria.split(";").map(s => s.trim()).filter(Boolean);
        return `<details class="cnv-sv-reasoning">
          <summary>AnnotSV 評分依據 <span class="cnv-sv-reasoning-score"><strong>Score:</strong> ${escapeHtml(score)}</span></summary>
          <ul class="cnv-sv-reasoning-list">${
            items.map(s => `<li><code>${escapeHtml(s)}</code></li>`).join("")
          }</ul>
        </details>`;
      })()
    : `<div class="cnv-sv-reasoning"><strong>Score:</strong> ${escapeHtml(score)}</div>`;

  // Disease-related gene count = genes with a non-empty OMIM
  // phenotype text. OMIM_morbid would be the strictest signal but we
  // don't always parse it; OMIM_phenotype non-empty is a reliable
  // proxy and is the same field we already render in the gene table.
  const totalGenes = (v.gene_count != null) ? v.gene_count : (v.genes || []).length;
  const diseaseGenes = (v.genes || []).filter(g => g.omim_phenotype).length;
  const geneCountText = totalGenes != null
    ? `${totalGenes}（疾病相關：${diseaseGenes}）`
    : "—";

  return `<div class="cnv-sv-detail-box">
    <div class="cnv-sv-detail-row">
      <span><strong>ACMG:</strong> ${acmgSelect}</span>
      <span><strong>涵蓋基因數:</strong> ${escapeHtml(geneCountText)}</span>
      <span><strong>基因型:</strong> ${escapeHtml(v.GT || "—")} (${escapeHtml(zyg)})${cn}</span>
      <span><strong>Filter:</strong> ${escapeHtml(filter)}</span>
      <span><strong>Qual:</strong> ${qual}</span>
    </div>
    ${reasoning}
  </div>`;
}

function _renderCnvSvGeneTable(v, id) {
  // Backend trims `genes` to the visible-table set (≤10 rows + any
  // in-panel overflow), and ships the long tail in `genes_compact`
  // with only the chip-display fields. `genes_total` is the original
  // gene_count so the section header reads correctly.
  const genes = v.genes || [];
  const genesCompact = v.genes_compact || [];
  const total = (v.genes_total != null) ? v.genes_total : genes.length + genesCompact.length;
  if (!total) return "";

  const picked = getEdit(id, "report_genes") || {};
  const _fmtW = w => (w % 1 === 0) ? String(w | 0) : Number(w).toFixed(1);
  const _firstLine = s => (s || "").split("\n")[0] || "";
  const rowHtml = (g) => {
    const checked = picked[g.gene] ? "checked" : "";
    const triggerMark = g.in_panel ? `<span class="pheno-star" title="HPO/panel match">★</span>` : "";
    const omimCell = g.omim_id
      ? `<a href="https://www.omim.org/entry/${escapeAttr(g.omim_id)}" target="_blank" rel="noopener">${escapeHtml(g.omim_id)}</a>`
      : "—";
    // AnnotSV emits Overlapped_CDS_percent as 0..100 already (saw
    // 100 → "10000%" pre-fix). Treat the value as the percent itself,
    // no extra ×100.
    const cdsPct = (g.overlap_cds_pct != null)
      ? `${Math.round(Number(g.overlap_cds_pct))}%` : "—";
    // Pheno reads as `matched/total` so the reviewer sees how many
    // input HPO/panel weights implicate this gene. Falls back to "—"
    // when phenotype isn't configured (denominator 0).
    const pheno = (g.pheno_total && g.pheno_total > 0)
      ? `${_fmtW(g.pheno_matched || 0)}/${_fmtW(g.pheno_total)}`
      : "—";
    const inh     = g.omim_inheritance || "";
    const phenAll = g.omim_phenotype   || "";
    const sameGeneCell = g.gene
      ? `<button class="same-gene-btn" data-gene="${escapeAttr(g.gene)}" type="button" title="列出 ${escapeAttr(g.gene)} 的所有 SNV/Indel + CNV/SV 變異">搜尋同基因</button>`
      : "";
    return `<tr class="${g.in_panel ? "gene-row-in-panel" : ""}" data-gene="${escapeAttr(g.gene || "")}">
      <td class="gene-pick-cell">
        <input type="checkbox" class="gene-pick" data-id="${escapeAttr(id)}" data-gene="${escapeAttr(g.gene || "")}" ${checked} title="勾選=放進報告" />
      </td>
      <td><strong>${escapeHtml(g.gene || "?")}</strong>${triggerMark}</td>
      <td>${escapeHtml(g.tx || "")}</td>
      <td>${escapeHtml(g.location || "")}</td>
      <td>${cdsPct}</td>
      <td class="gene-clip-cell" data-full="${escapeAttr(inh)}" title="點此展開">${escapeHtml(inh) || "—"}</td>
      <td>${omimCell}</td>
      <td class="gene-clip-cell" data-full="${escapeAttr(phenAll)}" title="點此展開">${escapeHtml(_firstLine(phenAll)) || "—"}</td>
      <td>${pheno}</td>
      <td class="gene-search-cell">${sameGeneCell}</td>
    </tr>`;
  };

  const tableHead = `<thead><tr>
    <th></th><th>Gene</th><th>Tx</th><th>Location</th><th>CDS%</th>
    <th>Inheritance</th><th>OMIM</th><th>Phenotype</th><th>Pheno</th><th></th>
  </tr></thead>`;
  const relevantGenes = genes.filter(g => g.in_panel);
  const hiddenFullGenes = genes.filter(g => !g.in_panel);
  const visibleRows = relevantGenes.map(rowHtml).join("");

  // Overflow body is rendered lazily on first <details> open. For
  // SVs that span 1500+ genes, eagerly building the chip DOM was
  // adding ~100 ms per card even though the panel was hidden.
  const genesOverflow = v.genes_overflow || [];
  const overflowCount = hiddenFullGenes.length + genesOverflow.length + genesCompact.length;
  let overflowHtml = "";
  if (overflowCount) {
    overflowHtml = `<details class="cnv-sv-gene-overflow" data-id="${escapeAttr(id)}" data-rendered="0">
      <summary class="muted">展開其餘 ${overflowCount} 個基因…</summary>
      <div class="gene-overflow-body"></div>
    </details>`;
  }

  return `<div class="cnv-sv-section">
    <div class="cnv-sv-section-title">基因 (${total})</div>
    <table class="cnv-sv-gene-table">${tableHead}<tbody>${visibleRows}</tbody></table>
    ${overflowHtml}
  </div>`;
}

function _renderCnvSvOverlap(v) {
  // Type-specific filter: a deletion only meaningfully overlaps loss
  // pathogenic regions; a duplication only gain regions; everything
  // else (INV / INS / TRA) shows all three so reviewers can pick.
  const allowed = new Set();
  if (v.sv_type === "DEL") allowed.add("p_loss");
  else if (v.sv_type === "DUP") allowed.add("p_gain");
  else { allowed.add("p_loss"); allowed.add("p_gain"); allowed.add("p_ins"); }

  // Each block clamps to 2 visible lines via CSS line-clamp (the
  // `\n`-split approach broke when AnnotSV puts the entire phen text
  // on one wrapped line). A toggle button below each block flips a
  // `.expanded` class to reveal the rest.
  const groups = [];
  for (const [key, label] of [["p_loss", "P_loss"], ["p_gain", "P_gain"], ["p_ins", "P_ins"]]) {
    if (!allowed.has(key)) continue;
    const p = v[key];
    if (!p || (!p.phens && !(p.sources || []).length)) continue;
    const phenLine = p.phens ? `<div class="cnv-sv-overlap-phen">${escapeHtml(p.phens)}</div>` : "";
    const sources = p.sources || [];
    const sourcesHtml = sources.length
      ? `<div class="muted cnv-sv-overlap-sources">${sources.map(escapeHtml).join("； ")}</div>`
      : "";
    groups.push(`<div class="cnv-sv-overlap-row">
      <div class="cnv-sv-overlap-head"><strong>${label}:</strong></div>
      <div class="cnv-sv-overlap-content">${phenLine}${sourcesHtml}</div>
      <button type="button" class="cnv-sv-overlap-toggle">▸ 展開全部</button>
    </div>`);
  }
  if (!groups.length) {
    return `<div class="cnv-sv-section">
      <div class="cnv-sv-section-title">已知致病區域重疊</div>
      <div class="cnv-sv-overlap-empty muted">無已知致病區域重疊</div>
    </div>`;
  }
  return `<div class="cnv-sv-section">
    <div class="cnv-sv-section-title">已知致病區域重疊</div>
    ${groups.join("")}
  </div>`;
}

function _renderCnvSvBenign(v) {
  // Type-specific filter mirrors the pathogenic-overlap one:
  // DEL → B_loss, DUP → B_gain, INV / INS / TRA → all four blocks.
  const allowed = new Set();
  if (v.sv_type === "DEL") allowed.add("b_loss");
  else if (v.sv_type === "DUP") allowed.add("b_gain");
  else { allowed.add("b_loss"); allowed.add("b_gain"); allowed.add("b_ins"); allowed.add("b_inv"); }

  const groups = [];
  for (const [key, label] of [["b_loss","B_loss"], ["b_gain","B_gain"], ["b_ins","B_ins"], ["b_inv","B_inv"]]) {
    if (!allowed.has(key)) continue;
    const b = v[key];
    if (!b || (!b.sources?.length && !b.coords?.length)) continue;
    const afHead = (b.af_max != null)
      ? ` <span class="muted" style="font-size:11px">max AF ${Number(b.af_max).toFixed(4)}</span>`
      : "";
    // Pair source+coord+AF when we can; fall back to source-only line
    // if the lengths don't agree (defensive — AnnotSV usually emits
    // them in lock-step).
    const lines = [];
    const n = Math.max(b.sources.length, b.coords.length, (b.afs || []).length);
    for (let i = 0; i < n; i++) {
      const src   = b.sources[i] || "";
      const coord = b.coords[i] || "";
      const af    = (b.afs && b.afs[i]) || "";
      const segs = [];
      if (src)   segs.push(escapeHtml(src));
      if (coord) segs.push(`<span class="muted">${escapeHtml(coord)}</span>`);
      if (af)    segs.push(`<span class="muted">AF ${escapeHtml(af)}</span>`);
      if (segs.length) lines.push(segs.join(" · "));
    }
    const sourcesHtml = lines.length
      ? `<div class="muted cnv-sv-overlap-sources">${lines.join("； ")}</div>`
      : "";
    groups.push(`<div class="cnv-sv-overlap-row cnv-sv-benign-row">
      <div class="cnv-sv-overlap-head"><strong>${label}:</strong>${afHead}</div>
      <div class="cnv-sv-overlap-content">${sourcesHtml}</div>
      <button type="button" class="cnv-sv-overlap-toggle">▸ 展開全部</button>
    </div>`);
  }
  if (!groups.length) {
    return `<div class="cnv-sv-section">
      <div class="cnv-sv-section-title">已知良性區域重疊</div>
      <div class="cnv-sv-overlap-empty muted">無已知良性區域重疊</div>
    </div>`;
  }
  return `<div class="cnv-sv-section">
    <div class="cnv-sv-section-title">已知良性區域重疊</div>
    ${groups.join("")}
  </div>`;
}

function _renderCnvSvComment(v, id) {
  const comment = (getEdit(id, "comment") || "");
  return `<div class="cnv-sv-section cnv-sv-comment">
    <div class="cnv-sv-section-title">Comment</div>
    <textarea class="cnv-sv-comment-text" data-id="${escapeAttr(id)}" rows="2" placeholder="備註">${escapeHtml(comment)}</textarea>
  </div>`;
}

function _renderCnvSvDisease(v, id) {
  const disease = (getEdit(id, "disease") || "");
  return `<div class="cnv-sv-section cnv-sv-disease">
    <div class="cnv-sv-section-title">Disease</div>
    <textarea class="cnv-sv-disease-text" data-id="${escapeAttr(id)}" rows="2" placeholder="疾病名稱">${escapeHtml(disease)}</textarea>
  </div>`;
}

function renderCnvSvCard(v, id, opts = {}) {
  const card = document.createElement("div");
  card.className = "variant-card cnv-sv-card";
  card.dataset.id = id;
  card.innerHTML = `
    ${_renderCnvSvHeader(v, id, opts)}
    ${_renderCnvSvDetailBox(v, id)}
    ${_renderCnvSvGeneTable(v, id)}
    ${_renderCnvSvOverlap(v)}
    ${_renderCnvSvBenign(v)}
    ${_renderCnvSvDisease(v, id)}
    ${_renderCnvSvComment(v, id)}
  `;
  return card;
}

// CNV/SV-specific edit hooks. These piggy-back on the existing
// state.reports.{status, edits} dicts the SNV cards use; AnnotSV_IDs
// and chr-pos-ref-alt SNV ids never collide so one flat namespace is
// fine. Selectors are scoped to .cnv-sv-card so the SNV handlers in
// renderVariantCard's setup don't double-fire. The status dropdown
// itself shares the .status-select class with SNV — its existing
// document-level handler updates state.reports.status keyed by id,
// which works for either kind of variant.
document.addEventListener("change", ev => {
  const t = ev.target;
  // Both CNV/SV cards and Mito cards use this listener block — match
  // either so changes on a Mito card fire correctly.
  const card = t.closest?.(".cnv-sv-card, .mito-card");
  if (!card) return;
  const id = card.dataset.id;
  if (!id) return;
  if (t.matches(".gene-pick")) {
    const picked = { ...(getEdit(id, "report_genes") || {}) };
    const gene = t.dataset.gene;
    if (t.checked) picked[gene] = true; else delete picked[gene];
    setEdit(id, "report_genes", picked);
    updateSaveHint();
  } else if (t.matches(".cnv-sv-acmg-select")) {
    setEdit(id, "ACMG_class_sv", t.value);
    t.classList.remove("sig-p","sig-lp","sig-vus","sig-lb","sig-b");
    const next = SV_ACMG_SIG_CLASS[Number(t.value)];
    if (next) t.classList.add(next);
    updateSaveHint();
  } else if (t.matches(".mito-acmg-select")) {
    setEdit(id, "ACMG_classification_mito", t.value);
    t.classList.remove("sig-p","sig-lp","sig-vus","sig-lb","sig-b");
    const sig = classifySignificance(t.value);
    if (sig) t.classList.add(sig);
    try { renderReportSections(); } catch (_e) {}
    updateSaveHint();
  }
});

document.addEventListener("input", ev => {
  const t = ev.target;
  if (!t.matches?.(".cnv-sv-comment-text, .cnv-sv-disease-text")) return;
  const id = t.dataset.id;
  if (!id) return;
  setEdit(id, t.matches(".cnv-sv-disease-text") ? "disease" : "comment", t.value);
  updateSaveHint();
});

// Click on a truncated cell (Inheritance / Phenotype) → expand it
// to show the full text. Click again to collapse. Each cell tracks
// its own state via .gene-clip-expanded so phen and inh expand
// independently.
document.addEventListener("click", ev => {
  const cell = ev.target.closest?.(".cnv-sv-gene-table .gene-clip-cell");
  if (!cell) return;
  if (cell.classList.contains("gene-clip-expanded")) {
    cell.classList.remove("gene-clip-expanded");
    cell.textContent = (cell.dataset.full || "").split("\n")[0] || "—";
  } else {
    cell.classList.add("gene-clip-expanded");
    cell.textContent = cell.dataset.full || "";
  }
});

// Lazy-render the gene-overflow chip body on first <details> open.
// SVs that span thousands of genes shipped ~600 KB of chip DOM up
// front; deferring it keeps card render fast and only pays the cost
// when the reviewer actually expands the section. Toggle events
// don't bubble, so the listener attaches in capture phase.
document.addEventListener("toggle", ev => {
  const det = ev.target;
  if (!det || !det.classList?.contains("cnv-sv-gene-overflow")) return;
  if (!det.open) return;
  if (det.dataset.rendered === "1") return;
  const id = det.dataset.id;
  const v = _cnvSvVariantById(id);
  if (!v) return;
  const overflowFull = [
    ...(v.genes || []).filter(g => !g.in_panel),
    ...(v.genes_overflow || []),
  ];
  const compact      = v.genes_compact  || [];
  const body = det.querySelector(".gene-overflow-body");
  if (!body) return;
  // In-panel rows beyond the visible cap stay in full table format
  // (so reviewers can still see Tx / Location / Phenotype for them).
  // Non-in-panel rows collapse to compact chips since they were
  // shipped without those fields.
  const picked = getEdit(id, "report_genes") || {};
  const _fmtW = w => (w % 1 === 0) ? String(w | 0) : Number(w).toFixed(1);
  const _firstLine = s => (s || "").split("\n")[0] || "";
  const fullRowHtml = (g) => {
    const checked = picked[g.gene] ? "checked" : "";
    const triggerMark = g.in_panel ? `<span class="pheno-star" title="HPO/panel match">★</span>` : "";
    const omimCell = g.omim_id
      ? `<a href="https://www.omim.org/entry/${escapeAttr(g.omim_id)}" target="_blank" rel="noopener">${escapeHtml(g.omim_id)}</a>`
      : "—";
    const cdsPct = (g.overlap_cds_pct != null)
      ? `${Math.round(Number(g.overlap_cds_pct))}%` : "—";
    const pheno = (g.pheno_total && g.pheno_total > 0)
      ? `${_fmtW(g.pheno_matched || 0)}/${_fmtW(g.pheno_total)}`
      : "—";
    const inh     = g.omim_inheritance || "";
    const phenAll = g.omim_phenotype   || "";
    const sameGeneCell = g.gene
      ? `<button class="same-gene-btn" data-gene="${escapeAttr(g.gene)}" type="button" title="列出 ${escapeAttr(g.gene)} 的所有 SNV/Indel + CNV/SV 變異">搜尋同基因</button>`
      : "";
    return `<tr class="${g.in_panel ? "gene-row-in-panel" : ""}" data-gene="${escapeAttr(g.gene || "")}">
      <td class="gene-pick-cell"><input type="checkbox" class="gene-pick" data-id="${escapeAttr(id)}" data-gene="${escapeAttr(g.gene || "")}" ${checked} title="勾選=放進報告" /></td>
      <td><strong>${escapeHtml(g.gene || "?")}</strong>${triggerMark}</td>
      <td>${escapeHtml(g.tx || "")}</td>
      <td>${escapeHtml(g.location || "")}</td>
      <td>${cdsPct}</td>
      <td class="gene-clip-cell" data-full="${escapeAttr(inh)}" title="點此展開">${escapeHtml(inh) || "—"}</td>
      <td>${omimCell}</td>
      <td class="gene-clip-cell" data-full="${escapeAttr(phenAll)}" title="點此展開">${escapeHtml(_firstLine(phenAll)) || "—"}</td>
      <td>${pheno}</td>
      <td class="gene-search-cell">${sameGeneCell}</td>
    </tr>`;
  };
  const tableHead = `<thead><tr>
    <th></th><th>Gene</th><th>Tx</th><th>Location</th><th>CDS%</th>
    <th>Inheritance</th><th>OMIM</th><th>Phenotype</th><th>Pheno</th><th></th>
  </tr></thead>`;
  const fullTable = overflowFull.length
    ? `<table class="cnv-sv-gene-table">${tableHead}<tbody>${overflowFull.map(fullRowHtml).join("")}</tbody></table>`
    : "";
  const chipBlock = compact.length
    ? `<div class="gene-overflow-chips">${compact.map(g =>
        `<span class="gene-overflow-chip${g.in_panel ? " gene-row-in-panel" : ""}">${escapeHtml(g.gene || "?")}${
          g.omim_id ? ` <a href="https://www.omim.org/entry/${escapeAttr(g.omim_id)}" target="_blank" rel="noopener" class="muted">${escapeHtml(g.omim_id)}</a>` : ""
        }</span>`
      ).join("")}</div>`
    : "";
  body.innerHTML = fullTable + chipBlock;
  det.dataset.rendered = "1";
}, true);

// 致病區域重疊 expand/collapse: each row has its content wrapped in
// a CSS-line-clamped div; this toggle flips the expanded class and
// updates the button label.
document.addEventListener("click", ev => {
  const btn = ev.target.closest?.(".cnv-sv-overlap-toggle");
  if (!btn) return;
  const row = btn.closest(".cnv-sv-overlap-row");
  if (!row) return;
  const expanded = row.classList.toggle("expanded");
  btn.textContent = expanded ? "▾ 收合" : "▸ 展開全部";
});

// Sidebar nav: clicking a button with data-target scrolls the matching
// card into view. Cards declare an id (scroll-margin-top keeps the
// landing position below the topbar). On narrow viewports we also
// auto-collapse the sidebar after the click so it doesn't sit on top
// of the freshly-revealed content.
document.addEventListener("click", ev => {
  const link = ev.target.closest(".sidebar-link[data-target]");
  if (!link) return;
  const target = document.getElementById(link.dataset.target);
  if (!target) return;
  target.scrollIntoView({ behavior: "smooth", block: "start" });
  if (window.matchMedia("(max-width: 768px)").matches) {
    document.body.classList.add("sidebar-collapsed");
    _setSidebarToggleAria(false);
  }
});

// Sidebar open/close. Default: expanded on desktop, collapsed on
// mobile. localStorage remembers an explicit toggle so the choice
// persists across reloads.
function _setSidebarToggleAria(open) {
  const btn = document.getElementById("btn-sidebar-toggle");
  if (btn) btn.setAttribute("aria-expanded", open ? "true" : "false");
}
(function initSidebar() {
  const stored = localStorage.getItem("ngs-sidebar");
  const isMobile = window.matchMedia("(max-width: 768px)").matches;
  const open = stored === null ? !isMobile : stored === "open";
  document.body.classList.toggle("sidebar-collapsed", !open);
  // aria-expanded is set on DOMContentLoaded since the button isn't
  // guaranteed to be in the DOM when this IIFE runs (module script).
  document.addEventListener("DOMContentLoaded", () => _setSidebarToggleAria(open), { once: true });
})();
// Secondary findings header toggle. Plain triangle-button toggling
// .open on the sibling body div — no animation, matches the
// lightweight visual style the reviewer asked for.
document.querySelector(".secondary-findings-toggle")?.addEventListener("click", (ev) => {
  const btn = ev.currentTarget;
  const body = btn.nextElementSibling;
  const expand = btn.getAttribute("aria-expanded") !== "true";
  btn.setAttribute("aria-expanded", expand ? "true" : "false");
  body?.classList.toggle("open", expand);
});

document.getElementById("btn-sidebar-toggle")?.addEventListener("click", () => {
  const collapsed = document.body.classList.toggle("sidebar-collapsed");
  localStorage.setItem("ngs-sidebar", collapsed ? "collapsed" : "open");
  _setSidebarToggleAria(!collapsed);
});

// Tally how many currently-loaded SNV variants pass each high-level gene
// scope flag so reviewers can tell whether the filters are doing work.
function updateInPanelCount() {
  const variants = state.data?.variants || {};
  const total = Object.keys(variants).length;
  const inPanel = Object.values(variants).filter(v => v.in_panel).length;
  const diseaseAssociated = Object.values(variants).filter(v => v.disease_associated).length;
  const el = document.getElementById("in-panel-count");
  if (el) el.textContent = total ? `(${inPanel} / ${total})` : "";
  const daEl = document.getElementById("disease-associated-count");
  if (daEl) daEl.textContent = total ? `(${diseaseAssociated} / ${total})` : "";
}

function renderPharmcatBlock(hostId) {
  const host = document.getElementById(hostId);
  host.innerHTML = "";
  host.style.display = "";

  const pc = state.data?.pgx || state.data?.pharmcat || {};
  const reportView = pc.report_view || {};
  const reportGenes = Array.isArray(reportView.report_genes) ? reportView.report_genes : [];
  const loading = !!state.pgxPending;
  const isReportPreview = hostId === "sec-pharmcat";

  // Mark the block as data-bearing so the CSS only paints the gray header
  // background when PharmCAT actually returned something.
  host.classList.toggle("has-data", reportGenes.length > 0);

  const wasOpen = toggledBlocks.has(hostId)
    ? host.dataset.wasOpen === "1"
    : currentSampleTestType() === "TITAN-WGS" || hostId === "cat-pharmcat-c";
  host.dataset.wasOpen = wasOpen ? "1" : "0";

  const header = document.createElement("div");
  header.className = "block-header" + (wasOpen ? " open" : "");
  header.innerHTML = `
    <span><span class="arrow"></span><span class="title">PGx / PharmCAT</span></span>
    <span class="count">${loading ? "…" : `${reportGenes.length} genes`}</span>`;
  host.appendChild(header);

  const body = document.createElement("div");
  body.className = "block-body" + (wasOpen ? " open" : "");
  if (loading) {
    body.innerHTML = `<div class="muted">PGx 載入中…</div>`;
  } else if (!reportGenes.length) {
    body.innerHTML = `<div class="muted">尚無 PGx / PharmCAT 結果。</div>`;
  } else {
    body.innerHTML = isReportPreview
      ? renderPharmcatReportBody(pc)
      : renderPharmcatAnalysisBody(pc);
  }
  host.appendChild(body);

  header.addEventListener("click", () => {
    const open = body.classList.toggle("open");
    header.classList.toggle("open", open);
    host.dataset.wasOpen = open ? "1" : "0";
    toggledBlocks.add(hostId);
  });
}

function _pgxMetaHtml(pc) {
  const meta = [
    pc.timestamp ? `Generated ${escapeHtml(pc.timestamp)}` : "",
    pc.pharmcat_version ? `PharmCAT ${escapeHtml(pc.pharmcat_version)}` : "",
    pc.data_version ? `Data ${escapeHtml(pc.data_version)}` : "",
  ].filter(Boolean).join(" · ");
  return meta ? `<div class="muted pharmcat-ts">${meta}</div>` : "";
}

function _pgxSourceWarning(pc) {
  return pc.pharmcat_available
    ? ""
    : `<div class="pgx-source-warning">缺少 PharmCAT JSON；除 MT-RNR1 可由 TSV 補值外，其餘基因不會使用 TSV 代替。</div>`;
}

function _pgxTable(headers, rows, className = "") {
  const head = headers.map(value => `<th>${escapeHtml(value)}</th>`).join("");
  const body = rows.map(row => `<tr>${row.map(value => `<td>${escapeHtml(value ?? "")}</td>`).join("")}</tr>`).join("");
  return `<div class="pgx-table-wrap"><table class="pgx-report-table ${escapeAttr(className)}"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>`;
}

function renderPharmcatReportBody(pc) {
  const view = pc.report_view || {};
  const summary = view.summary || {};
  const actionableGenes = Array.isArray(summary.actionable_genes) ? summary.actionable_genes : [];
  const actionCategories = Array.isArray(view.action_categories) ? view.action_categories : [];
  const summaryRows = Array.isArray(view.summary_rows) ? view.summary_rows : [];
  const genotypeRows = Array.isArray(view.genotype_rows) ? view.genotype_rows : [];
  const overview = actionCategories.length
    ? actionCategories.map(category => `
        <div class="pgx-overview-row">
          <strong>${escapeHtml(category.action || "")}</strong>
          <span>${escapeHtml((category.drugs || []).join("、"))}</span>
        </div>`).join("")
    : `<div class="muted">本次檢測未發現符合目前回報規則之明確臨床可應用藥物基因體結果。</div>`;
  const summaryTable = summaryRows.length
    ? _pgxTable(
        ["藥物", "基因", "建議處置", "建議依據及等級"],
        summaryRows.map(row => [row.drug, row.gene, row.action, row.source_level]),
        "pgx-summary-table",
      )
    : `<div class="muted">（無符合摘要回報規則的藥物）</div>`;
  const genotypeTable = _pgxTable(
    ["基因", "基因型", "表型"],
    genotypeRows.map(row => [row.gene, row.genotype, row.phenotype]),
    "pgx-genotype-table",
  );
  const drugCount = Number(summary.actionable_drugs || 0);
  const summaryDrugCount = Number(summary.summary_drugs || 0);
  const appendixOnly = Number(summary.appendix_only_drugs || 0);
  const intro = drugCount
    ? `本次檢測發現 ${drugCount} 項具臨床用藥參考價值的藥物結果${actionableGenes.length ? `，涉及 ${actionableGenes.join("、")} 基因` : ""}；摘要表列出 ${summaryDrugCount} 項${appendixOnly ? `，另 ${appendixOnly} 項僅列於分析區完整用藥建議` : ""}。`
    : "此處依健檢報告回報規則顯示 PGx 結果。";
  return `${_pgxSourceWarning(pc)}${_pgxMetaHtml(pc)}
    <div class="pgx-report-intro">${escapeHtml(intro)}</div>
    <div class="pgx-section">
      <h4 class="pgx-heading">用藥建議概覽</h4>
      <div class="pgx-overview">${overview}</div>
      <div class="pgx-residual-note">其餘未列之藥物，未發現符合本報告回報規則之明確處方調整建議。</div>
    </div>
    <div class="pgx-section">
      <h4 class="pgx-heading">藥物建議摘要</h4>
      ${summaryTable}
    </div>
    <div class="pgx-section">
      <h4 class="pgx-heading">基因型與表現型</h4>
      ${genotypeTable}
    </div>
  `;
}

function _hasPharmcatGeneDetails(g) {
  const details = g?.details || {};
  const meaningful = value => {
    const text = typeof value === "string" ? value.trim().toLowerCase() : "";
    return !!text && !["unknown", "no result", "n/a", "na", ".", "-", "—"].includes(text);
  };
  return !!(
    meaningful(details.allele1_name) || meaningful(details.allele1_function)
    || meaningful(details.allele2_name) || meaningful(details.allele2_function)
    || (details.variants || []).length
    || (details.uncalled || []).length
    || (details.messages || []).some(message => typeof message === "string" && message.trim())
  );
}

function renderPharmcatAnalysisGene(pc, row) {
  const gene = row.gene || "";
  const payload = pc.genes?.[gene] || { gene, details: {} };
  const mtRnr1Supplement = gene === "MT-RNR1"
    ? [payload.mtrn1_risk, payload.outside_caller, payload.notes].filter(Boolean).join(" · ")
    : "";
  const detail = renderPharmcatGeneDetails(payload, row);
  const attentionClass = row.requires_attention ? " pgx-gene-head-attention" : "";
  return `<div class="pgx-analysis-gene">
    <div class="pgx-gene-head${attentionClass}">
      <strong>${escapeHtml(gene)}</strong>
      <span>${escapeHtml(row.genotype || "—")}</span>
      <span class="pgx-gene-phenotype">${escapeHtml(row.phenotype || "No phenotype assigned")}</span>
      ${mtRnr1Supplement ? `<span class="pgx-evidence-badge">${escapeHtml(mtRnr1Supplement)}</span>` : ""}
    </div>
    ${detail}
  </div>`;
}

function renderPharmcatAnalysisBody(pc) {
  const view = pc.report_view || {};
  const fullRows = Array.isArray(view.full_recommendations) ? view.full_recommendations : [];
  const genotypeRows = Array.isArray(view.analysis_genes)
    ? view.analysis_genes
    : (Array.isArray(view.genotype_rows) ? view.genotype_rows : []);
  const actionCategories = Array.isArray(view.analysis_action_categories)
    ? view.analysis_action_categories
    : (Array.isArray(view.action_categories) ? view.action_categories : []);
  const overview = actionCategories.map(category => {
    const drugs = Array.isArray(category.drugs) ? category.drugs.filter(Boolean) : [];
    return `<div class="pgx-overview-row">
      <strong>${escapeHtml(category.action || "")}</strong>
      <span>${escapeHtml(drugs.length ? drugs.join("、") : "無")}</span>
    </div>`;
  }).join("") || `<div class="muted">（無 PGx 用藥分類資料）</div>`;
  const recommendationTable = fullRows.length
    ? _pgxTable(
        ["藥物", "基因與表型", "CPIC/FDA 建議"],
        fullRows.map(row => [row.drug, row.gene_phenotype, row.recommendation]),
        "pgx-full-recommendation-table",
      )
    : `<div class="muted">（無符合完整回報規則的用藥建議）</div>`;
  const genes = genotypeRows.map(row => renderPharmcatAnalysisGene(pc, row)).join("");
  const globalMessages = (pc.messages || []).filter(message => typeof message === "string" && message.trim());
  const messages = globalMessages.length
    ? `<div class="pgx-global-messages"><strong>PharmCAT messages</strong><ul>${globalMessages.map(message => `<li>${escapeHtml(message)}</li>`).join("")}</ul></div>`
    : "";
  return `${_pgxSourceWarning(pc)}${_pgxMetaHtml(pc)}
    <div class="pgx-section">
      <h4 class="pgx-heading">用藥建議概覽</h4>
      <div class="pgx-overview">${overview}</div>
    </div>
    <div class="pgx-section">
      <h4 class="pgx-heading">完整用藥建議</h4>
      ${recommendationTable}
    </div>
    <div class="pgx-section">
      <h4 class="pgx-heading">CPIC Level A 基因詳細結果</h4>
      <div class="pgx-analysis-genes">${genes}</div>
    </div>
    ${messages}`;
}

function _pharmcatMeaningfulDetail(value) {
  if (typeof value !== "string") return "";
  const text = value.trim();
  return ["", "unknown", "no result", "n/a", "na", ".", "-", "—"].includes(text.toLowerCase())
    ? ""
    : text;
}

function _pharmcatStarAllelesCell(value) {
  const alleles = String(value || "").split(",").map(item => item.trim()).filter(Boolean);
  if (!alleles.length) return "—";
  if (alleles.length <= 4) return escapeHtml(alleles.join(", "));
  return `<details class="pgx-star-alleles"><summary>${alleles.length} 個可能的 star allele 定義</summary><div>${escapeHtml(alleles.join(", "))}</div></details>`;
}

function renderPharmcatGeneDetails(g, analysisRow = {}) {
  const fnLine = (name, fn) => {
    name = _pharmcatMeaningfulDetail(name);
    fn = _pharmcatMeaningfulDetail(fn);
    if (!name && !fn) return "";
    const left  = name ? escapeHtml(name) : "—";
    const right = fn   ? ` <span class="muted">(${escapeHtml(fn)})</span>` : "";
    return `<div>${left}${right}</div>`;
  };
  const details = g.details || {};
  const vrows = (details.variants || []).map(v => `
    <tr>
      <td>${escapeHtml(v.rsid || "—")}</td>
      <td>${escapeHtml(v.chr || "")}:${escapeHtml(String(v.pos ?? ""))}</td>
      <td>${escapeHtml(v.call || "")}</td>
      <td>${_pharmcatStarAllelesCell(v.alleles)}</td>
    </tr>`).join("");
  const uncalledValues = (details.uncalled || []).filter(value => typeof value === "string" && value.trim());
  const uncalled = uncalledValues.length
    ? `<div class="muted">Uncalled: ${uncalledValues.map(escapeHtml).join(", ")}</div>`
    : "";
  const messageValues = (details.messages || []).filter(value => typeof value === "string" && value.trim());
  const messages = messageValues.length
    ? `<div class="muted">Messages: ${messageValues.map(escapeHtml).join(", ")}</div>`
    : "";
  const variants = vrows
    ? `<div class="pgx-star-alleles-note">Call 是病人實際檢出的基因型；「可對應之 star alleles」列出此位點可參與定義的 PharmCAT star alleles，不代表病人同時具有全部 alleles。</div>
       <table class="pharmcat-variants"><thead><tr><th>rsID</th><th>Position</th><th>Call</th><th>可對應之 star alleles</th></tr></thead><tbody>${vrows}</tbody></table>`
    : "";
  const recommendationRows = Array.isArray(analysisRow.recommendation_rows)
    ? analysisRow.recommendation_rows
    : [];
  const relatedDrugs = Array.isArray(analysisRow.related_drugs)
    ? analysisRow.related_drugs.filter(Boolean)
    : [];
  const recommendations = recommendationRows.length
    ? _pgxTable(
        ["藥物", "基因與表型", "CPIC/FDA 建議"],
        recommendationRows.map(row => [row.drug, row.gene_phenotype, row.recommendation]),
        "pgx-full-recommendation-table pgx-gene-recommendation-table",
      )
    : `<div class="muted pgx-related-drugs">相關藥物：${escapeHtml(relatedDrugs.length ? relatedDrugs.join("、") : "無")}</div>`;
  const alleleEvidence = _hasPharmcatGeneDetails(g)
    ? `<div class="pharmcat-alleles">
         ${fnLine(details.allele1_name, details.allele1_function)}
         ${fnLine(details.allele2_name, details.allele2_function)}
       </div>${uncalled}${messages}${variants}`
    : `<div class="muted pgx-allele-placeholder">（無其他 PharmCAT allele 證據）</div>`;
  return `<details class="pgx-detail"><summary>詳細</summary>
            <div class="pgx-gene-recommendations">${recommendations}</div>
            ${alleleEvidence}
          </details>`;
}

function renderAll() {
  if (!state.data) return;
  updateWelcomeVisibility();
  renderSampleMeta();
  renderGeneticCounseling();
  renderClinicalDescription();
  renderPhenotype();
  renderVersionPicker();
  renderDeadZoneCard();
  renderComment();
  renderReportSections();
  renderCandidateSections();
  renderDiseaseAssociatedReportWarning();
  updateSaveHint();
}

// ---------- Welcome / version notes --------------------------------

function updateWelcomeVisibility() {
  document.getElementById("welcome-card")?.classList.toggle("hidden", !!state.data);
}

async function loadWelcomeVersion() {
  const host = document.getElementById("welcome-version-content");
  if (!host) return;
  try {
    const resp = await fetch("./VERSION.md", { cache: "no-store" });
    if (!resp.ok) throw new Error(`${resp.status} ${resp.statusText}`);
    const text = await resp.text();
    host.innerHTML = renderSimpleMarkdown(text);
  } catch (e) {
    host.innerHTML = `<div class="muted">版本資訊載入失敗：${escapeHtml(e.message || String(e))}</div>`;
  }
}

function renderSimpleMarkdown(text) {
  const lines = String(text || "").replace(/\r\n/g, "\n").split("\n");
  const out = [];
  let listOpen = false;
  const closeList = () => {
    if (listOpen) {
      out.push("</ul>");
      listOpen = false;
    }
  };
  for (const raw of lines) {
    const line = raw.trim();
    if (!line) {
      closeList();
      continue;
    }
    const h = line.match(/^(#{1,3})\s+(.+)$/);
    if (h) {
      closeList();
      const level = h[1].length;
      out.push(`<h${level}>${renderMarkdownInline(h[2])}</h${level}>`);
      continue;
    }
    const item = line.match(/^-\s+(.+)$/);
    if (item) {
      if (!listOpen) {
        out.push("<ul>");
        listOpen = true;
      }
      out.push(`<li>${renderMarkdownInline(item[1])}</li>`);
      continue;
    }
    closeList();
    out.push(`<p>${renderMarkdownInline(line)}</p>`);
  }
  closeList();
  return out.join("");
}

function renderMarkdownInline(text) {
  return escapeHtml(text).replace(/`([^`]+)`/g, "<code>$1</code>");
}

// Auto-save state. Every dirty edit schedules a debounced background
// save through saveChanges({silent:true}); the hint line under the
// 💾 buttons reflects the current phase (dirty / saving / saved /
// error). beforeunload below catches the rare case where the user
// closes the tab during the debounce window.
let _autoSaveTimer = null;
let _saveInflight = false;
let _saveError = "";
let _lastSavedAt = null;   // Date of the most recent successful save, for the hint label

function scheduleAutoSave(delayMs = 1500) {
  if (!state.currentLIS) return;
  clearTimeout(_autoSaveTimer);
  _autoSaveTimer = setTimeout(_doAutoSave, delayMs);
}

async function _doAutoSave() {
  if (!state.dirty || !state.currentLIS || _saveInflight) return;
  await saveChanges({ silent: true });
}

async function flushPendingSave(timeoutMs = 10000) {
  clearTimeout(_autoSaveTimer);
  const deadline = Date.now() + timeoutMs;
  while (_saveInflight && Date.now() < deadline) {
    await new Promise(resolve => setTimeout(resolve, 50));
  }
  if (state.dirty && !_saveInflight) {
    await saveChanges({ silent: true });
  }
  while (_saveInflight && Date.now() < deadline) {
    await new Promise(resolve => setTimeout(resolve, 50));
  }
  return !state.dirty && !_saveInflight;
}

function updateSaveHint() {
  let msg;
  if (_saveError) {
    msg = "⚠ 儲存失敗：" + _saveError;
  } else if (_saveInflight) {
    msg = "儲存中…";
  } else if (state.dirty) {
    msg = "有變更（自動儲存中…）";
    scheduleAutoSave();
  } else {
    // Append the timestamp of the last successful save so the
    // reviewer knows whether the auto-save fired recently. Locale
    // toLocaleTimeString output is e.g. "下午2:30:15" / "14:30:15".
    msg = _lastSavedAt
      ? `已儲存（${_lastSavedAt.toLocaleTimeString()}）`
      : "已儲存";
  }
  // Update every save-hint span on the page (top / mid / bottom).
  const txt = msg;
  document.querySelectorAll(".js-save-hint").forEach(el => { el.textContent = txt; });
}

// Native browser confirmation when the user tries to leave with
// unsaved edits or a save still in flight. Auto-save normally fires
// after 1.5 s of inactivity, but a closed tab between keystrokes can
// still drop the most-recent change.
window.addEventListener("beforeunload", (ev) => {
  if (state.dirty || _saveInflight) {
    ev.preventDefault();
    ev.returnValue = "";
  }
});

// ---------- Event wiring -------------------------------------------

document.addEventListener("click", ev => {
  const t = ev.target;
  if (t.matches(".js-btn-save")) {
    saveChanges();
  } else if (t.matches("#btn-close-bottom")) {
    collapseCandidateSections();
  } else if (t.matches(".btn-export-clinical")) {
    exportDiagnosticDocx();
  } else if (t.matches(".btn-export-health")) {
    exportHealthDocx();
  } else if (t.matches(".btn-print-report")) {
    printReportCards();
  } else if (t.matches(".btn-export-html, .btn-export-screening")) {
    // Other export targets (analysis HTML / screening PDF) are not yet
    // ported from the legacy GitHub-Pages tool.
    alert("此匯出格式尚未實作。");
  } else if (t.closest(".clinical-header")) {
    toggleCollapsibleCard(t.closest(".clinical-header"));
  } else if (t.closest(".btn-copy")) {
    ev.stopPropagation();
    copyToClipboard(t.closest(".btn-copy"));
  } else if (t.matches(".consequence-toggle")) {
    ev.stopPropagation();
    const wrap = t.closest(".consequence-multi");
    const rest = wrap?.querySelector(".consequence-rest");
    if (rest) {
      const willHide = !rest.classList.contains("hidden");
      rest.classList.toggle("hidden", willHide);
      t.textContent = willHide ? "▾" : "▴";
    }
  } else if (t.matches(".btn-more")) {
    ev.stopPropagation();
    toggleVariantExtras(t);
  } else if (t.matches(".tag-remove")) {
    ev.stopPropagation();
    removeTag(t.dataset.tag);
  } else if (t.matches(".btn-add-manual")) {
    ev.stopPropagation();
    addManualVariant(t.dataset.status);
  } else if (t.matches(".btn-remove-manual")) {
    ev.stopPropagation();
    removeManualVariant(t.dataset.mid);
  } else if (t.matches(".disease-pick, .mito-disease-pick")) {
    // Don't let clicking the checkbox also toggle its <details> container.
    ev.stopPropagation();
  } else if (t.matches(".disease-collapse")) {
    // "▴ 收合" at the bottom of the expanded yellow detail box —
    // closes the parent <details>.
    ev.stopPropagation();
    const det = t.closest("details.disease-row");
    if (det) det.open = false;
  }
});

function copyToClipboard(btn) {
  const text = btn.dataset.copy || "";
  if (!text) return;
  // Save the SVG icon as innerHTML so we can restore it after the
  // brief ✓ / ✗ flash; textContent would strip the child <svg>.
  const orig = btn.innerHTML;
  const flash = (mark, ms) => {
    btn.textContent = mark;
    setTimeout(() => { btn.innerHTML = orig; }, ms);
  };
  // navigator.clipboard requires a secure context (HTTPS or localhost).
  // The hospital intranet serves this app over plain HTTP, so we fall
  // back to the legacy textarea + execCommand approach when the modern
  // API is unavailable or rejects.
  const legacyCopy = () => {
    try {
      const ta = document.createElement("textarea");
      ta.value = text;
      ta.setAttribute("readonly", "");
      ta.style.position = "fixed";
      ta.style.top = "0";
      ta.style.left = "0";
      ta.style.opacity = "0";
      document.body.appendChild(ta);
      ta.focus();
      ta.select();
      const ok = document.execCommand("copy");
      document.body.removeChild(ta);
      flash(ok ? "✓" : "✗", ok ? 900 : 1200);
    } catch {
      flash("✗", 1200);
    }
  };
  if (navigator.clipboard && window.isSecureContext) {
    navigator.clipboard.writeText(text)
      .then(() => flash("✓", 900))
      .catch(legacyCopy);
  } else {
    legacyCopy();
  }
}

function toggleVariantExtras(btn) {
  const card = btn.closest(".variant-card");
  if (!card) return;
  // A card may have multiple .more-extras blocks (e.g. one under the
  // Total/Variant/Pheno column for Exomiser+LIRICAL, another under the
  // AlphaMissense column for LoGoFunc / MaxEntScan / PDIVAS). The single
  // More button toggles them in lockstep.
  const blocks = card.querySelectorAll(".more-extras");
  if (!blocks.length) return;
  const willHide = !blocks[0].classList.contains("hidden");
  blocks.forEach(b => b.classList.toggle("hidden", willHide));
  btn.textContent = willHide ? "▾ More" : "▴ Less";
}

// Legacy panel-status radio inputs cannot normally be unchecked. Modern
// secondary ✓ chips are checkboxes; this remains only for old stragglers.
document.addEventListener("pointerdown", ev => {
  const chip = ev.target.closest?.(".status-radio-chip");
  const input = chip?.querySelector('input[type="radio"]');
  if (input) input.dataset.wasChecked = input.checked ? "1" : "0";
});

document.addEventListener("click", ev => {
  const chip = ev.target.closest?.(".status-radio-chip");
  const input = chip?.querySelector('input[type="radio"]');
  if (!input || input.dataset.wasChecked !== "1") return;
  ev.preventDefault();
  const wrap = input.closest(".status-radio");
  if (!wrap) return;
  const panel = wrap.dataset.panel;
  if (panel) setPanelStatus(wrap.dataset.id, panel, "");
  else       setStatus(wrap.dataset.id, "");
});

document.addEventListener("change", ev => {
  const t = ev.target;
  // Status chips replace the old <select.status-select>. Main chips
  // use checkboxes (C + 0 may coexist); secondary ✓ chips are checkboxes too.
  if (t.matches('.status-radio input[type="radio"], .status-radio input[type="checkbox"]')) {
    const wrap = t.closest(".status-radio");
    if (!wrap) return;
    const panel = wrap.dataset.panel;
    if (panel) setPanelStatus(wrap.dataset.id, panel, t.checked ? t.value : "");
    else       toggleStatus(wrap.dataset.id, t.value, t.checked);
  } else if (t.matches(".status-select")) {
    // Legacy fallback for any straggling <select> instance.
    const panel = t.dataset.panel;
    if (panel) setPanelStatus(t.dataset.id, panel, t.value);
    else       setStatus(t.dataset.id, t.value);
  } else if (t.matches("#m-sry-confirmed")) {
    state.reports.sry_confirmed = t.checked;
    state.dirty = true;
    updateSaveHint();
  } else if (t.matches("#m-category")) {
    state.reports.category = t.value || null;
    state.dirty = true;
    updateSaveHint();
  } else if (t.matches(".acmg-class")) {
    setEdit(t.dataset.id, "ACMG_classification", t.value);
    _refreshAcmgSourceHints(t.dataset.id);
    updateSaveHint();
    // Re-apply significance color to match the edited value
    t.classList.remove(...SIG_CLASSES);
    const cls = classifySignificance(t.value);
    if (cls) t.classList.add(cls);
    renderReportSections();
    renderCandidateSections();
  } else if (t.matches(".acmg-score")) {
    setEdit(t.dataset.id, "ACMG_score", t.value);
    _refreshAcmgSourceHints(t.dataset.id);
    updateSaveHint();
  } else if (t.matches(".acmg-crit")) {
    setEdit(t.dataset.id, "ACMG_criteria", t.value);
    _refreshAcmgSourceHints(t.dataset.id);
    updateSaveHint();
  } else if (t.matches(".variant-comment")) {
    setEdit(t.dataset.id, "comment", t.value);
    updateSaveHint();
  } else if (t.matches(".transcript-select")) {
    setEdit(t.dataset.id, "selected_transcript_key", t.value);
    renderAll();
    updateSaveHint();
  } else if (t.matches(".disease-pick")) {
    const picked = { ...(getEdit(t.dataset.id, "report_diseases") || {}) };
    const idx = t.dataset.idx;
    if (t.checked) picked[idx] = true;
    else           delete picked[idx];
    setEdit(t.dataset.id, "report_diseases", picked);
    _syncVariantCheckboxes(".disease-pick", t.dataset.id, idx, t.checked, t);
    updateSaveHint();
  } else if (t.matches(".mito-disease-pick")) {
    const picked = { ...(getEdit(t.dataset.id, "report_diseases_clinvar") || {}) };
    const idx = t.dataset.idx;
    if (t.checked) picked[idx] = true;
    else           delete picked[idx];
    setEdit(t.dataset.id, "report_diseases_clinvar", picked);
    _syncVariantCheckboxes(".mito-disease-pick", t.dataset.id, idx, t.checked, t);
    updateSaveHint();
  }
});

function autoGrow(ta) {
  if (!ta) return;
  ta.style.height = "auto";
  ta.style.height = ta.scrollHeight + "px";
}

document.addEventListener("input", ev => {
  const t = ev.target;
  if (t.matches("#clinical-text")) {
    state.reports.clinical_description = t.value;
    state.dirty = true;
    updateSaveHint();
    autoGrow(t);
  } else if (t.matches("#counseling-text")) {
    state.reports.genetic_counseling = t.value;
    state.dirty = true;
    updateSaveHint();
    autoGrow(t);
  } else if (t.matches("#comment-text")) {
    state.reports.comment = t.value;
    state.dirty = true;
    updateSaveHint();
  } else if (t.matches(".variant-comment")) {
    setEdit(t.dataset.id, "comment", t.value);
    _syncEditControls(t.dataset.id, "comment", t.value, t);
    updateSaveHint();
  } else if (t.matches(".acmg-score")) {
    setEdit(t.dataset.id, "ACMG_score", t.value);
    _syncEditControls(t.dataset.id, "ACMG_score", t.value, t);
    _refreshAcmgSourceHints(t.dataset.id);
    updateSaveHint();
  } else if (t.matches(".acmg-crit")) {
    setEdit(t.dataset.id, "ACMG_criteria", t.value);
    _syncEditControls(t.dataset.id, "ACMG_criteria", t.value, t);
    _refreshAcmgSourceHints(t.dataset.id);
    updateSaveHint();
  } else if (t.matches(".manual-position")) {
    updateManualVariant(t.dataset.mid, "position", t.value);
    // Keep the adjacent 📋 button copying the latest position string.
    const btn = t.parentElement?.querySelector(".btn-copy");
    if (btn) btn.dataset.copy = t.value;
    updateSaveHint();
  } else if (t.matches(".manual-comment")) {
    updateManualVariant(t.dataset.mid, "comment", t.value);
    updateSaveHint();
  } else if (t.matches(".manual-disease")) {
    updateManualVariant(t.dataset.mid, "disease", t.value);
    updateSaveHint();
  }
});

// Tag input: Enter or comma commits the typed value and clears the field;
// picking a datalist suggestion fires 'change' which we commit too.
document.addEventListener("keydown", ev => {
  if (!ev.target.matches(".tag-input")) return;
  if (ev.key === "Enter" || ev.key === ",") {
    ev.preventDefault();
    const v = ev.target.value;
    ev.target.value = "";
    addTag(v);
    setTimeout(() => {
      const fresh = document.querySelector(".tag-input");
      if (fresh) fresh.focus();
    }, 0);
  } else if (ev.key === "Escape") {
    ev.target.value = "";
  }
});
document.addEventListener("change", ev => {
  if (ev.target.matches(".tag-input")) {
    const v = ev.target.value;
    ev.target.value = "";
    addTag(v);
  }
});

// Toggle any collapsible header (Clinical presentation, Comment, …).
// Body element is always the next sibling of the header by convention,
// and the .card section wraps both — that's where wasOpen lives.
function toggleCollapsibleCard(header) {
  const card = header.closest(".card");
  const body = header.nextElementSibling;
  if (!card || !body) return;
  const open = !(card.dataset.wasOpen === "1");
  card.dataset.wasOpen = open ? "1" : "0";
  header.classList.toggle("open", open);
  body.classList.toggle("open", open);
  toggledBlocks.add(card.id);
  // The clinical textarea auto-grows; resize it to fit existing
  // content the moment the body becomes visible (scrollHeight is 0
  // while display:none).
  if (open) {
    const ta = body.querySelector("#clinical-text, #counseling-text");
    if (ta) requestAnimationFrame(() => autoGrow(ta));
  }
}

function collapseCandidateSections() {
  const host = document.getElementById("category-sections");
  host.querySelectorAll(".cat-block:not(.tier-panel)").forEach(block => {
    block.dataset.wasOpen = "0";
    toggledBlocks.add(block.id);
    block.querySelector(".block-header")?.classList.remove("open");
    block.querySelector(".block-body")?.classList.remove("open");
  });
  // Also close the currently visible tier panels (SNV + CNV/SV) by
  // deselecting their tabs. The tab strips stay put so the user can
  // reopen any tier with one click.
  activeTierTab = null;
  activeCnvSvTab = null;
  document.querySelectorAll(".tier-tab").forEach(b => b.classList.remove("active"));
  document.querySelectorAll(".tier-panel").forEach(p => p.classList.remove("active"));
}

// ---------- Save: write reports JSON to GitHub via Contents API ----

async function ghGetSha(path) {
  const token = getToken();
  if (!token) throw new Error("No GitHub token");
  const url = `${API_BASE}/contents/${encodePath(path)}?ref=${encodeURIComponent(BRANCH)}`;
  const resp = await fetch(url, {
    headers: { Authorization: `token ${token}`, Accept: "application/vnd.github+json" },
    cache: "no-store",
  });
  if (resp.status === 404) return null;
  if (!resp.ok) throw new Error(`${resp.status} ${resp.statusText} on GET ${path}`);
  const j = await resp.json();
  return j.sha || null;
}

function toBase64Utf8(s) {
  // Chunked to avoid blowing the call-stack on large HTML payloads.
  const bytes = new TextEncoder().encode(s);
  let binary = "";
  const CHUNK = 0x8000;
  for (let i = 0; i < bytes.length; i += CHUNK) {
    binary += String.fromCharCode.apply(null, bytes.subarray(i, i + CHUNK));
  }
  return btoa(binary);
}

async function ghPutContent(path, text, message) {
  const token = getToken();
  if (!token) throw new Error("No GitHub token");
  const b64 = toBase64Utf8(text);

  const put = async (sha) => {
    const body = { message, content: b64, branch: BRANCH };
    if (sha) body.sha = sha;
    return fetch(`${API_BASE}/contents/${encodePath(path)}`, {
      method: "PUT",
      headers: {
        Authorization: `token ${token}`,
        Accept: "application/vnd.github+json",
        "Content-Type": "application/json",
      },
      body: JSON.stringify(body),
    });
  };

  let sha = await ghGetSha(path);
  let resp = await put(sha);
  // Conflict: refresh sha once and retry
  if (resp.status === 409 || resp.status === 422) {
    sha  = await ghGetSha(path);
    resp = await put(sha);
  }
  if (!resp.ok) {
    const txt = await resp.text();
    throw new Error(`${resp.status} ${resp.statusText}: ${txt}`);
  }
  return resp.json();
}

// Same Contents API push but accepts already-base64'd content (so binary
// payloads like a PDF can go up without going through toBase64Utf8).
async function ghPutBinary(path, base64Content, message) {
  const token = getToken();
  if (!token) throw new Error("No GitHub token");
  const put = async (sha) => {
    const body = { message, content: base64Content, branch: BRANCH };
    if (sha) body.sha = sha;
    return fetch(`${API_BASE}/contents/${encodePath(path)}`, {
      method: "PUT",
      headers: {
        Authorization: `token ${token}`,
        Accept: "application/vnd.github+json",
        "Content-Type": "application/json",
      },
      body: JSON.stringify(body),
    });
  };
  let sha = await ghGetSha(path);
  let resp = await put(sha);
  if (resp.status === 409 || resp.status === 422) {
    sha  = await ghGetSha(path);
    resp = await put(sha);
  }
  if (!resp.ok) {
    const txt = await resp.text();
    throw new Error(`${resp.status} ${resp.statusText}: ${txt}`);
  }
  return resp.json();
}

async function ghPutJSON(path, obj, message) {
  return ghPutContent(path, JSON.stringify(obj, null, 2), message);
}

function encodePath(p) {
  return p.split("/").map(encodeURIComponent).join("/");
}

async function saveChanges(opts = {}) {
  // opts.silent: don't pop an alert on failure (auto-save calls this
  // with silent=true; manual 💾 clicks pass nothing → loud failure).
  if (!state.currentLIS || _saveInflight) return;
  if (!state.dirty) return;
  _saveInflight = true;
  _saveError = "";
  // Cancel any pending debounced save; this call is the save.
  clearTimeout(_autoSaveTimer);

  const saveBtns = document.querySelectorAll(".js-btn-save");
  const setBusy = b => saveBtns.forEach(btn => { btn.disabled = b; });

  setBusy(true);
  updateSaveHint();
  try {
    const statuses = state.reports.status || {};
    const hasAutoCausative = Object.values(statuses).some(v => v === "1");
    const hasManualCausative = (state.reports.manual_variants || []).some(
      m => m.status === "1" && (m.position || "").trim() !== ""
    );
    state.reports.yield = (hasAutoCausative || hasManualCausative) ? 1 : 0;
    const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
    const sid = row?.sample_id || state.currentLIS;
    await apiPut(`/samples/${encodeURIComponent(sid)}/report`, state.reports);
    state.dirty = false;
    _lastSavedAt = new Date();
  } catch (e) {
    _saveError = e.message || "未知錯誤";
    if (!opts.silent) alert("儲存失敗：" + e.message);
  } finally {
    _saveInflight = false;
    setBusy(false);
    updateSaveHint();
  }
}

// ---------- Export: static analysis HTML → GitHub ------------------

async function exportAnalysisHTML() {
  if (!state.currentLIS) return;
  const btns = document.querySelectorAll(".btn-export-html");
  const setBusy = b => btns.forEach(x => { x.disabled = b; });
  const hint = msg => {
    document.querySelectorAll(".js-save-hint").forEach(el => { el.textContent = msg; });
  };

  setBusy(true);
  hint("產生匯出檔…");
  try {
    // 1) Inline CSS (fetch current stylesheet)
    let css = "";
    try {
      const cssResp = await fetch("./style.css", { cache: "no-store" });
      if (cssResp.ok) css = await cssResp.text();
    } catch {/* ignore — export still works, just unstyled */}

    // 2) Clone the main app area and freeze interactive widgets
    const orig  = document.getElementById("app");
    const clone = orig.cloneNode(true);
    freezeForExport(orig, clone);

    // 3) Wrap in a self-contained HTML document
    const meta = state.data?.meta || {};
    const title = `VCF Analysis — ${meta.LIS_ID || state.currentLIS}`;
    const html = buildExportHTML({ css, bodyInner: clone.outerHTML, title, meta });

    // 4) Trigger a local download (always-on per the user's preference)
    //    and push the same content to GitHub. Local download is fired
    //    first so the user sees immediate feedback even if the network
    //    push later fails.
    const fname = `${state.currentLIS}.html`;
    downloadBlob(new Blob([html], { type: "text/html;charset=utf-8" }), fname);
    const path = `output/analysis_html/${fname}`;
    await ghPutContent(path, html, `export: analysis HTML for ${state.currentLIS}`);
    hint(`已下載 + 匯出 → ${path}`);
  } catch (e) {
    hint("");
    alert("匯出失敗：" + e.message);
  } finally {
    setBusy(false);
  }
}

function freezeForExport(orig, clone) {
  // cloneNode() copies markup, not live input state — snapshot values from
  // the ORIGINAL DOM in document order, then apply by index to the clone.
  const origControls  = Array.from(orig.querySelectorAll("input, select, textarea"));
  const cloneControls = Array.from(clone.querySelectorAll("input, select, textarea"));
  const values = origControls.map(el => {
    if (el.tagName === "INPUT" && el.type === "checkbox") {
      return { kind: "checkbox", checked: el.checked };
    }
    if (el.tagName === "SELECT") {
      const opt = el.options[el.selectedIndex];
      return { kind: "select", text: opt ? (opt.textContent || opt.value || "") : "" };
    }
    return { kind: "text", value: el.value };
  });

  // Mirror the live page's default expansion: only Causative / Other are open.
  // The rest stay collapsed but remain clickable in the export thanks to the
  // small toggle script inlined by buildExportHTML().
  ["sec-causative", "sec-other"].forEach(id => {
    const host = clone.querySelector("#" + id);
    if (!host) return;
    host.querySelector(".block-header")?.classList.add("open");
    host.querySelector(".block-body")?.classList.add("open");
  });

  // Remove the search UI and all save-rows — no interactive controls in export.
  clone.querySelectorAll(".save-row").forEach(el => el.remove());
  clone.querySelectorAll("#q-lis, #q-lis-dropdown, #search-status").forEach(el => el.remove());
  // If the first card is the (id-less) search card, drop it.
  const firstCard = clone.querySelector(".card");
  if (firstCard && !firstCard.id) firstCard.remove();

  // Freeze form controls to static spans/divs, using the snapshotted values.
  cloneControls.forEach((el, i) => {
    const v   = values[i];
    const tag = el.tagName;
    let replacement;
    if (v.kind === "checkbox") {
      replacement = document.createElement("span");
      replacement.textContent = v.checked ? "☑" : "☐";
    } else if (v.kind === "select") {
      replacement = document.createElement("span");
      replacement.textContent = (v.text || "").trim() || "—";
    } else if (tag === "TEXTAREA") {
      replacement = document.createElement("div");
      replacement.textContent = v.value || "";
    } else {
      replacement = document.createElement("span");
      replacement.textContent = v.value || "—";
    }
    replacement.className = `export-static ${el.className || ""}`.trim();
    if (tag === "TEXTAREA") replacement.classList.add("export-multiline");
    el.replaceWith(replacement);
  });
}

function buildExportHTML({ css, bodyInner, title, meta }) {
  const exportCss = `
/* --- Export-only overrides --- */
body { background: #fff; }
.topbar { background: #1f2328; color: #fff; padding: 12px 20px; }
.topbar h1 { margin: 0; font-size: 18px; font-weight: 600; }
.export-static {
  display: inline-block;
  border: none !important;
  background: transparent !important;
  padding: 0 !important;
  font: inherit;
  color: inherit;
  min-width: 0;
  vertical-align: baseline;
}
.export-static.export-multiline {
  display: block;
  white-space: pre-wrap;
  overflow-wrap: anywhere;
}
.export-static.status-select { font-weight: 600; }
.export-static.acmg-crit { white-space: pre-wrap; display: block; }
.export-static.variant-comment,
.export-static.acmg-class,
.export-static.acmg-score { font-family: ui-monospace, monospace; }
.export-banner {
  text-align: center;
  color: #6a737d;
  font-size: 12px;
  padding: 12px 0;
  border-top: 1px solid #d0d7de;
  margin-top: 24px;
}
`;
  const exportedAt = new Date().toISOString();
  const headerLine = [meta.LIS_ID, meta.Name, meta.MRN].filter(Boolean).map(escapeHtml).join(" — ");
  return `<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>${escapeHtml(title)}</title>
<style>${css}
${exportCss}
</style>
</head>
<body>
<header class="topbar"><h1>VCF Analysis${headerLine ? " — " + headerLine : ""}</h1></header>
${bodyInner}
<div class="export-banner">Exported ${escapeHtml(exportedAt)} · static snapshot</div>
<script>
// Click-to-toggle for collapsed sections (no other interactivity in the export).
document.addEventListener("click", function (e) {
  var h = e.target.closest(".block-header, .clinical-header");
  if (!h) return;
  var body = h.nextElementSibling;
  h.classList.toggle("open");
  if (body) body.classList.toggle("open");
});
</script>
</body>
</html>`;
}

// ---------- Export: clinical (TXT) + screening (PDF) reports -------
//
// Two more buttons in the save-row produce per-sample reports:
//   - Clinical TXT (the "診斷報告"), pushed to output/clinical_reports/
//   - Screening PDF (the "健檢報告"), pushed to output/screening_reports/
// Both follow the same "fetch state → render → ghPutContent" shape as
// exportAnalysisHTML above. Helpers are defined below; the main entry
// points are exportClinicalReport() and exportScreeningReport().

// VEP Consequence → 中文 + 報告解釋句. Fallback below for unknown terms.
const CONSEQUENCE_CH = {
  stop_gained:                     { label: "無義突變 (Nonsense)",
    explain: "此變異會形成過早的終止密碼子，使蛋白質轉譯提前終止，通常導致蛋白質功能喪失。" },
  missense_variant:                { label: "誤義突變 (Missense)",
    explain: "此變異使密碼子改變為不同的胺基酸，可能影響蛋白質的結構穩定性、功能或相互作用。" },
  frameshift_variant:              { label: "移碼突變 (Frameshift)",
    explain: "此變異改變了開放閱讀框，使下游胺基酸序列完全錯亂並通常產生過早的終止密碼子，常導致蛋白質功能喪失。" },
  stop_lost:                       { label: "終止密碼子失去 (Stop loss)",
    explain: "此變異使原本的終止密碼子消失，導致蛋白質轉譯延長，可能影響蛋白質功能。" },
  start_lost:                      { label: "起始密碼子失去 (Start loss)",
    explain: "此變異使蛋白質的起始密碼子消失，可能導致蛋白質無法正常轉譯。" },
  inframe_insertion:               { label: "框內插入 (Inframe insertion)",
    explain: "此變異在不改變閱讀框的前提下插入若干胺基酸，可能影響蛋白質結構或功能。" },
  inframe_deletion:                { label: "框內缺失 (Inframe deletion)",
    explain: "此變異在不改變閱讀框的前提下缺失若干胺基酸，可能影響蛋白質結構或功能。" },
  splice_donor_variant:            { label: "剪接供體位點變異 (Splice donor)",
    explain: "此變異位於 intron 的 5' 剪接位點，可能導致 intron 無法正確剪除，影響蛋白質序列。" },
  splice_acceptor_variant:         { label: "剪接受體位點變異 (Splice acceptor)",
    explain: "此變異位於 intron 的 3' 剪接位點，可能導致 intron 無法正確剪除，影響蛋白質序列。" },
  splice_region_variant:           { label: "剪接區變異 (Splice region)",
    explain: "此變異位於剪接位點附近，可能影響 intron 剪除的準確度。" },
  splice_donor_5th_base_variant:   { label: "剪接供體第 5 鹼基變異",
    explain: "此變異位於剪接供體下游第 5 個位置，可能干擾正常的剪接過程。" },
  splice_donor_region_variant:     { label: "剪接供體區變異",
    explain: "此變異位於剪接供體位點附近，可能影響剪接準確度。" },
  splice_polypyrimidine_tract_variant: { label: "聚嘧啶區變異 (Polypyrimidine tract)",
    explain: "此變異位於剪接受體上游的聚嘧啶區，可能影響剪接效率。" },
  protein_altering_variant:        { label: "蛋白質改變變異",
    explain: "此變異會造成蛋白質序列改變（非單純胺基酸取代），可能影響蛋白質功能。" },
  synonymous_variant:              { label: "同義變異 (Synonymous)",
    explain: "此變異不改變胺基酸序列，通常不影響蛋白質功能。" },
  intron_variant:                  { label: "內含子變異 (Intronic)",
    explain: "此變異位於內含子區，通常不影響蛋白質序列，但若靠近剪接位點仍可能干擾剪接。" },
  "5_prime_UTR_variant":           { label: "5' 非轉譯區變異",
    explain: "此變異位於基因 5' 非轉譯區，可能影響轉錄或轉譯效率。" },
  "3_prime_UTR_variant":           { label: "3' 非轉譯區變異",
    explain: "此變異位於基因 3' 非轉譯區，可能影響 mRNA 穩定性或轉譯調控。" },
  upstream_gene_variant:           { label: "上游基因區變異",
    explain: "此變異位於基因上游，可能影響轉錄調控。" },
  downstream_gene_variant:         { label: "下游基因區變異",
    explain: "此變異位於基因下游，臨床意義通常不明確。" },
  intergenic_variant:              { label: "基因間變異 (Intergenic)",
    explain: "此變異位於基因間區，臨床意義通常不明確。" },
  non_coding_transcript_exon_variant: { label: "非編碼轉錄本外顯子變異",
    explain: "此變異位於非編碼轉錄本的外顯子區，臨床意義通常不明確。" },
  mature_miRNA_variant:            { label: "成熟 miRNA 變異",
    explain: "此變異位於成熟 miRNA 序列，可能影響其調控功能。" },
  coding_sequence_variant:         { label: "編碼序列變異",
    explain: "此變異位於編碼序列，但具體影響需個別評估。" },
  TF_binding_site_variant:         { label: "轉錄因子結合位點變異",
    explain: "此變異位於轉錄因子結合位點，可能影響基因表達調控。" },
  regulatory_region_variant:       { label: "調控區域變異",
    explain: "此變異位於調控區域，可能影響基因表達。" },
};

// ACMG / ClinVar 5-tier classifier → 中文 (used in the "此為...之變異位點" sentence).
const ACMG_CH = {
  "Pathogenic":             "致病性",
  "Likely pathogenic":      "疑似致病性",
  "Uncertain significance": "不確定意義",
  "Likely benign":          "可能良性",
  "Benign":                 "良性",
  "Conflicting":            "意義分歧",
};

// OMIM inheritance code → 中文. Codes appear inside parens in Disease text.
const INHERITANCE_LABELS = {
  AD:  "體染色體顯性遺傳",
  AR:  "體染色體隱性遺傳",
  XLD: "性染色體顯性遺傳",
  XLR: "性染色體隱性遺傳",
  XL:  "性聯遺傳",
  YL:  "Y 染色體遺傳",
  MT:  "粒線體遺傳",
  DD:  "雙等位基因顯性遺傳",
  IC:  "細胞質遺傳",
};

const ZYG_CH = { het: "Heterozygous", hom: "Homozygous", hemi: "Hemizygous" };

// East-Asian-aware visual width — Chinese / Japanese / Korean glyphs occupy
// 2 monospace cells each, ASCII / Latin take 1.
function visualWidth(s) {
  let w = 0;
  for (const ch of String(s || "")) {
    const c = ch.codePointAt(0);
    if (c >= 0x1100 && (
      c <= 0x115F ||
      (c >= 0x2E80 && c <= 0x9FFF) ||
      (c >= 0xA000 && c <= 0xA4CF) ||
      (c >= 0xAC00 && c <= 0xD7A3) ||
      (c >= 0xF900 && c <= 0xFAFF) ||
      (c >= 0xFE30 && c <= 0xFE4F) ||
      (c >= 0xFF00 && c <= 0xFF60) ||
      (c >= 0xFFE0 && c <= 0xFFE6)
    )) w += 2;
    else w += 1;
  }
  return w;
}
function padToWidth(s, target) {
  return s + " ".repeat(Math.max(0, target - visualWidth(s)));
}
function chunkByWidth(s, target) {
  const out = [];
  let line = "", w = 0;
  for (const ch of String(s || "")) {
    const cw = visualWidth(ch);
    if (w + cw > target && line) { out.push(line); line = ""; w = 0; }
    line += ch; w += cw;
  }
  if (line) out.push(line);
  if (!out.length) out.push("");
  return out;
}

// Build a fixed-width ASCII table with === separators. Long cells wrap onto
// successive lines within the same column; rows are joined with a single
// " " between columns and one leading space.
function formatVariantTable(rows) {
  const header = ["基因", "結構", "核苷酸", "基因型", "ClinVar", "ACMG&AMP指引"];
  const widths = [8, 8, 30, 14, 14, 14];
  const totalW = widths.reduce((a, b) => a + b, 0) + widths.length; // +N joiner spaces
  const sep = "=".repeat(totalW);
  const renderLine = cells => " " + cells.map((c, i) => padToWidth(c, widths[i])).join(" ");
  const out = [sep, renderLine(header), sep];
  for (const row of rows) {
    const cellLines = row.map((c, i) => chunkByWidth(String(c == null ? "" : c), widths[i]));
    const maxL = Math.max(...cellLines.map(c => c.length || 1));
    for (let li = 0; li < maxL; li++) {
      out.push(renderLine(cellLines.map(c => c[li] || "")));
    }
  }
  out.push(sep);
  return out.join("\n");
}

// Disease helpers ---------------------------------------------------

// Which Disease{i} did the user tick on this variant card? First ticked,
// or fall back to Disease1 if nothing ticked.
function pickedDiseaseSlot(id, v) {
  const picked = (state.reports?.edits?.[id]?.report_diseases) || {};
  const idxs = Object.keys(picked).filter(k => picked[k]).map(Number)
    .filter(n => Number.isFinite(n)).sort((a, b) => a - b);
  for (const i of idxs) {
    const d = v[`Disease${i}`];
    if (d && d !== "NA") return { idx: i, text: d };
  }
  for (let i = 1; i <= 5; i++) {
    const d = v[`Disease${i}`];
    if (d && d !== "NA") return { idx: i, text: d };
  }
  return { idx: 1, text: "" };
}

// Disease text format: "<Name> (INH) [: description]" or "<Name>, somatic"
// — extract a clean disease name + the inheritance code(s). Allows
// comma + whitespace inside the parens so "(AR, DD)" parses as
// inheritance="AR, DD" instead of falling through to a verbatim suffix.
function diseaseInfo(text) {
  if (!text) return { name: "", inheritance: "" };
  const firstLine = String(text).split("\n")[0].trim();
  const inhMatch = firstLine.match(/\(([A-Z][A-Z?\/,\s]*)\)/);
  const inh = inhMatch ? inhMatch[1].trim() : "";
  let name = firstLine;
  if (inhMatch) name = name.slice(0, firstLine.indexOf(inhMatch[0]));
  name = name.replace(/[:,;]+\s*$/, "").trim();
  return { name, inheritance: inh };
}
function inheritanceCH(code) {
  if (!code) return "遺傳模式未明確";
  // Compound codes like "AD/AR" or "AR, DD" — translate each, join with "或".
  const parts = code.split(/[\/,]/).map(s => s.trim()).filter(Boolean);
  const labels = parts.map(p => INHERITANCE_LABELS[p] || p);
  return labels.join("或");
}

// HGVS = "<gene>:<transcript>:<cdna>[:<protein>]" — split into pieces.
function parseHGVS(hgvs) {
  const parts = String(hgvs || "").split(":");
  return {
    gene:       parts[0] || "",
    transcript: parts[1] || "",
    cdna:       parts[2] || "",
    protein:    parts[3] || "",
  };
}
function hgvsCellText(v) {
  const h = parseHGVS(v.HGVS);
  return h.protein ? `${h.cdna}(${h.protein})` : h.cdna;
}

// ACMG / ClinVar text helpers ---------------------------------------

function acmgClassCH(cls) {
  if (!cls) return "";
  return ACMG_CH[cls] || cls;
}
function consequenceEntry(consequence) {
  // VEP can emit multiple terms joined with "&" or ","; pick the first known.
  const terms = String(consequence || "").split(/[&,]/).map(s => s.trim());
  for (const t of terms) {
    if (CONSEQUENCE_CH[t]) return { term: t, ...CONSEQUENCE_CH[t] };
  }
  const first = terms[0] || "";
  return {
    term: first,
    label: first ? `${first} 變異` : "未分類變異",
    explain: "此變異的功能影響需個別評估。",
  };
}

// gnomAD AF — pick the first numeric available across the field-name
// variants emitted by hg38 / hg19 pipelines.
function variantGnomadAF(v) {
  for (const k of ["gnomad41_genome_AF", "gnomad41_exome_AF", "AF"]) {
    const x = parseFloat(v[k]);
    if (!isNaN(x)) return x;
  }
  return null;
}
function gnomadAFText(v) {
  const af = variantGnomadAF(v);
  if (af == null || af === 0) {
    return "該變異位點在族群資料庫 gnomAD 中未報導過發生率，顯示其為罕見變異位點。";
  }
  // Plain-decimal percent rendering. Pick precision from the magnitude of the
  // value (rather than fixed 2-sig-figs) so very small AFs don't fall into
  // scientific notation — e.g. 7e-7 → "0.00007%" not "7.0e-5%".
  const pct = af * 100;
  let pretty;
  if (pct >= 1) {
    pretty = pct.toFixed(2).replace(/\.?0+$/, "");
  } else {
    const expo = Math.floor(Math.log10(pct));   // pct=7e-5 → expo=-5
    const decimals = Math.max(2, -expo + 1);    // expo=-5 → 6 decimals
    pretty = pct.toFixed(decimals).replace(/0+$/, "").replace(/\.$/, "");
  }
  if (af < 0.001) {
    return `該變異位點在族群資料庫 gnomAD 中報導過發生率為 ${pretty}%，顯示其為罕見變異位點。`;
  }
  return `該變異位點在族群資料庫 gnomAD 中報導過發生率為 ${pretty}%。`;
}

function variantClinSig(v) {
  // ClinVar exports use underscores in CLNSIG values (e.g. "Likely_pathogenic",
  // "Pathogenic/Likely_pathogenic"). Render with spaces in the report.
  const raw = v.CLNSIG || v.CLINSIG || v.CLNSIGn || "";
  const sig = String(raw).trim();
  return sig === "." ? "" : sig.replace(/_/g, " ");
}
function clinvarText(v) {
  const sig = variantClinSig(v);
  if (!sig) return "在疾病資料庫 (ClinVar) 中未被報導過。";
  return `在疾病資料庫 (ClinVar) 中此變異位點被報導為「${sig}」。`;
}
function acmgGuidelineText(v) {
  const cls = v.ACMG_classification;
  if (!cls) return "目前無 ACMG 評測。";
  return `根據美國醫學遺傳學暨基因體學學會 (American College of Medical Genetics and Genomics) 與分子病理學學會 (Association for Molecular Pathology) 於 2015 年發表之準則，評測此變異位點為「${cls}」。`;
}

// ymd -> "2025年05月04日"
function ymdToCnDate(ymd) {
  const s = String(ymd || "");
  const m = s.match(/^(\d{4})(\d{2})(\d{2})$/);
  if (!m) return s;
  return `${m[1]}年${m[2]}月${m[3]}日`;
}
function todayYmd() {
  const d = new Date();
  const pad = n => String(n).padStart(2, "0");
  return `${d.getFullYear()}${pad(d.getMonth() + 1)}${pad(d.getDate())}`;
}

// Patient phenotype list (二、檢驗套組) — render HPO rows as
// "Name (HP:id)" and panel rows as just the panel name.
function phenotypeSummaryCH(list) {
  const out = [];
  for (const r of list || []) {
    const ph    = (r.phenotype || "").trim();
    const label = (r.label || r.hpo_name || "").trim() || ph;
    if (ph.startsWith("HP:")) {
      out.push(label !== ph ? `${label} (${ph})` : ph);
    } else if (label) {
      out.push(label);
    }
  }
  return out.join("、") || "—";
}

// Word-wrap a paragraph at a target visual width with a fixed indent.
// Visual-width based (not whitespace-tokenised) so Chinese sentences wrap
// at the right column even when they contain no spaces.
function wrapText(text, target, indent) {
  const pad = " ".repeat(indent);
  const lines = [];
  let cur = "", w = 0;
  for (const ch of String(text || "")) {
    const cw = visualWidth(ch);
    if (w + cw > target && cur) {
      lines.push(pad + cur);
      cur = ""; w = 0;
    }
    cur += ch; w += cw;
  }
  if (cur) lines.push(pad + cur);
  return lines;
}

// Render the per-variant block: title line, ASCII table, the two numbered
// remarks, then the descriptive paragraph. `kind` is "causative" or "other"
// (it only changes the second remark sentence).
function renderVariantBlock(vid, v, kind) {
  const out = [];
  const h = parseHGVS(v.HGVS);
  const gene = v.gene_symbol || h.gene || "";
  const transcript = h.transcript || "—";

  // Title
  out.push("    " + (transcript ? `${gene} (${transcript})` : gene));

  // Table
  const tableText = formatVariantTable([[
    gene,
    v.exon_or_intron || "—",
    hgvsCellText(v),
    ZYG_CH[v.zygosity] || v.zygosity || "—",
    variantClinSig(v) || "—",
    v.ACMG_classification || "—",
  ]]);
  for (const ln of tableText.split("\n")) out.push("    " + ln);

  // Numbered remarks (use the user-picked Disease for inheritance + name)
  const dis = pickedDiseaseSlot(vid, v);
  const info = diseaseInfo(dis.text);
  const inhTxt = inheritanceCH(info.inheritance);
  const acmgTxt = acmgClassCH(v.ACMG_classification) || "—";
  const tail = kind === "causative"
    ? "與臨床症狀相關"
    : "無法完全解釋受檢者全部之臨床症狀，其臨床意義須由醫師配合其他相關資料進行最佳綜合判斷";
  out.push(`    1. ${gene}為${info.name || "—"}的致病基因之一，其遺傳模式屬於${inhTxt}。`);
  out.push(`    2. 此為${acmgTxt}之變異位點，${tail}。`);

  // Descriptive paragraph
  const cons = consequenceEntry(v.Consequence);
  const cdna = h.cdna || "";
  const prot = h.protein ? ` (${h.protein})` : "";
  const paragraph = [
    `在個案之檢體中，檢測到 1 個位於基因 ${gene} 的變異位點。`,
    `變異位點 ${cdna}${prot} 為${cons.label}，${cons.explain}`,
    gnomadAFText(v),
    clinvarText(v),
    acmgGuidelineText(v),
    "此報告僅供參考，臨床判斷仍應以病患的實際狀況為主。建議比對臨床表徵並進行父母親與家族成員之變異位點檢測，以釐清上述變異致病之可能性；根據家族成員變異位點檢測報告或相關資料庫更新，可能影響變異位點 ACMG 判讀結果。",
  ].join("");
  out.push("");
  out.push(...wrapText(paragraph, 76, 4));
  return out;
}

// Whole-document builder. Returns a TXT string. Pulls everything from
// state.data + state.reports (already loaded for the current sample).
function buildClinicalTXT() {
  const data = state.data || {};
  const meta = data.meta || {};
  const isWGS = isWgsTestType(meta.Test);
  const build = data.genome_build || "hg38";
  const clinvarDate = ymdToCnDate(data.clinvar_date) || "—";

  const statusMap = state.reports?.status || {};
  const causIds  = Object.keys(statusMap).filter(id => statusMap[id] === "1");
  const otherIds = Object.keys(statusMap).filter(id => statusMap[id] === "2");

  const lines = [];
  lines.push(`一、檢驗項目: 次世代定序${isWGS ? "全基因組" : "全外顯子"}定序檢測`);
  lines.push("");
  lines.push(`二、檢驗套組: ${phenotypeSummaryCH(data.patient_phenotype)}`);
  lines.push("");
  lines.push("三、檢測結果");
  lines.push("  檢體說明:");
  lines.push("    檢體類別：血液");
  lines.push("  綜合說明:");
  lines.push("");
  lines.push("    第一類：與臨床症狀相關基因之已知致病性變異位點");
  if (!causIds.length) {
    lines.push("    未找到與臨床症狀相關基因之已知致病性變異位點。");
  } else {
    for (const vid of causIds) {
      const v = data.variants?.[vid];
      if (!v) continue;
      lines.push("");
      lines.push(...renderVariantBlock(vid, v, "causative"));
    }
  }
  lines.push("");
  lines.push("    第二類：其他變異位點");
  if (!otherIds.length) {
    lines.push("    未找到其他變異位點。");
  } else {
    for (const vid of otherIds) {
      const v = data.variants?.[vid];
      if (!v) continue;
      lines.push("");
      lines.push(...renderVariantBlock(vid, v, "other"));
    }
  }

  lines.push("");
  lines.push("四、檢測方法說明");
  lines.push(`  1. 本次檢測使用次世代定序儀分析 (Illumina ${isWGS ? "NovaSeq X Plus" : "NextSeq 2000"})。`);
  lines.push("  2. 本次檢測變異位點的錯誤率 ≦ 0.1% (Phred-scaled Q score ≧ 30)。");
  lines.push(`  3. 本次檢測平均定序深度 ≧ ${isWGS ? "27.5X" : "50X"}。`);
  lines.push(...wrapText(
    "4. 本檢測僅能檢測出基因內單一核苷酸 (single nucleotide)、小片段的缺失或插入 (small indel)、大片段缺失 (deletion) 及擴增 (duplication)，無法檢測出轉位 (translocation)、倒轉 (inversion) 或其他複雜性結構變異 (complex structural variation)、組織特異性的鑲嵌 (tissue-specific mosaicism)、串聯重複 (tandem repeat) 以及未定序區域 (例如 promoter、intron)。",
    74, 2));
  lines.push(...wrapText("5. 本檢測報告僅供醫療專業人員參考，需配合其他相關臨床資料與家族成員之相關檢驗。", 74, 2));
  lines.push("  6. 目前次世代定序分子遺傳診斷皆屬研究性質。");

  lines.push("");
  lines.push("五、檢測結果注釋");
  lines.push(`  1. 本檢測結果比對參考序列為人類 ${build} 版本。`);
  lines.push(...wrapText(
    `2. ClinVar 及 ACMG&AMP 指引: 引用 ClinVar 資料庫截至 ${clinvarDate} 更新的註解，及美國醫學遺傳學暨基因體學學會 (ACMG) 與分子病理學學會 (AMP) 2015 年頒佈的指引，並且主要列入致病 (Pathogenic) 及疑似致病 (Likely pathogenic) 變異；其他類別變異經醫師判斷認為與疾病相關時亦可列入。`,
    74, 2));
  lines.push("  3. 參考資料:");
  lines.push("     a. 疾病資料庫: OMIM、ClinVar");
  lines.push(`     b. 族群資料庫: gnomAD (v4.1${isWGS ? " genome" : " exome"})`);
  lines.push("     c. 序列資料庫: RefSeqGene");
  lines.push("  4. 本次檢測基因包括:");
  const phenoGenes = (data.pheno_genes || []).slice().sort();
  if (phenoGenes.length) {
    lines.push(...wrapText(phenoGenes.join(", "), 74, 5));
  } else {
    lines.push("     —");
  }

  return lines.join("\n");
}

async function exportClinicalReport() {
  if (!state.currentLIS) return;
  const btns = document.querySelectorAll(".btn-export-clinical");
  const setBusy = b => btns.forEach(x => { x.disabled = b; });
  const hint = msg => {
    document.querySelectorAll(".js-save-hint").forEach(el => { el.textContent = msg; });
  };

  setBusy(true);
  hint("產生診斷報告…");
  try {
    const txt = buildClinicalTXT();
    const meta = state.data?.meta || {};
    const lis = state.currentLIS;
    const mrn = meta.MRN || "MRN";
    const fname = `${lis}_${mrn}_clinical_${todayYmd()}.txt`;
    downloadBlob(new Blob([txt], { type: "text/plain;charset=utf-8" }), fname);
    const path = `output/clinical_reports/${fname}`;
    await ghPutContent(path, txt, `export: clinical TXT for ${lis}`);
    hint(`已下載 + 匯出 → ${path}`);
  } catch (e) {
    hint("");
    alert("匯出失敗：" + e.message);
  } finally {
    setBusy(false);
  }
}

// Chinese-font loader for jsPDF. Fetches both Regular and Bold variants of
// Noto Sans TC, base64-encodes them in chunks (avoids the call-stack blowup
// when spreading a 5 MB byte array into String.fromCharCode), caches once
// per page lifetime. Bold is used for headings and KV labels; Regular for
// body / paragraphs. CDN paths get reorganised over time, so each chain
// has multiple fallbacks ordered by likelihood-of-success.
const PDF_FONT_NAME = "NotoSansTC";
const PDF_FONT_REGULAR_SOURCES = [
  "./fonts/NotoSansTC-Regular.ttf",
  "https://cdn.jsdelivr.net/gh/google/fonts@main/ofl/notosanstc/static/NotoSansTC-Regular.ttf",
  "https://cdn.jsdelivr.net/gh/google/fonts@main/ofl/notosanstc/NotoSansTC%5Bwght%5D.ttf",
  "https://raw.githubusercontent.com/google/fonts/main/ofl/notosanstc/static/NotoSansTC-Regular.ttf",
  "https://cdn.jsdelivr.net/npm/@fontsource/noto-sans-tc@4.5.13/files/noto-sans-tc-traditional-400-normal.ttf",
  "https://cdn.jsdelivr.net/gh/notofonts/noto-cjk@main/Sans/SubsetOTF/TC/NotoSansCJKtc-Regular.otf",
];
const PDF_FONT_BOLD_SOURCES = [
  "./fonts/NotoSansTC-Bold.ttf",
  "https://cdn.jsdelivr.net/gh/google/fonts@main/ofl/notosanstc/static/NotoSansTC-Bold.ttf",
  "https://cdn.jsdelivr.net/gh/google/fonts@main/ofl/notosanstc/NotoSansTC%5Bwght%5D.ttf",
  "https://raw.githubusercontent.com/google/fonts/main/ofl/notosanstc/static/NotoSansTC-Bold.ttf",
  "https://cdn.jsdelivr.net/npm/@fontsource/noto-sans-tc@4.5.13/files/noto-sans-tc-traditional-700-normal.ttf",
  "https://cdn.jsdelivr.net/gh/notofonts/noto-cjk@main/Sans/SubsetOTF/TC/NotoSansCJKtc-Bold.otf",
];
let _pdfFontRegular = null, _pdfFontBold = null;
async function _fetchFontB64(sources) {
  const tried = [];
  for (const url of sources) {
    try {
      const r = await fetch(url, { cache: "force-cache" });
      if (!r.ok) {
        tried.push(`${url} → ${r.status}`);
        console.warn(`[font] ${url} → ${r.status} ${r.statusText}`);
        continue;
      }
      const buf = new Uint8Array(await r.arrayBuffer());
      let bin = "";
      const CHUNK = 0x8000;
      for (let i = 0; i < buf.length; i += CHUNK) {
        bin += String.fromCharCode.apply(null, buf.subarray(i, i + CHUNK));
      }
      console.info(`[font] loaded ${url} (${(buf.length / 1024 / 1024).toFixed(1)} MB)`);
      return btoa(bin);
    } catch (e) {
      tried.push(`${url} → ${e.message}`);
      console.warn(`[font] ${url} → ${e.message}`);
    }
  }
  throw new Error("無法下載中文字型，所有來源都失敗：\n" + tried.join("\n"));
}
async function loadPdfFonts() {
  if (!_pdfFontRegular) _pdfFontRegular = await _fetchFontB64(PDF_FONT_REGULAR_SOURCES);
  if (!_pdfFontBold)    _pdfFontBold    = await _fetchFontB64(PDF_FONT_BOLD_SOURCES);
  return { regular: _pdfFontRegular, bold: _pdfFontBold };
}

// Gene-panel files (uploaded to docs/ alongside this app). One gene per
// line; rendered as a comma-joined list at the end of the PDF in small
// type. Cached after the first fetch.
const PANEL_FILES = [
  { title: "重大疾病風險篩檢基因清單", url: "./ACMG_SF_v3.3.txt" },
  { title: "帶因者篩檢基因清單",       url: "./carrier_mackenzie_1300+.txt" },
];
let _panelCache = null;
async function loadGenePanels() {
  if (_panelCache) return _panelCache;
  const out = [];
  for (const p of PANEL_FILES) {
    try {
      const r = await fetch(p.url, { cache: "force-cache" });
      if (!r.ok) { out.push({ title: p.title, genes: [] }); continue; }
      const txt = await r.text();
      const genes = txt.split(/\r?\n/).map(s => s.trim()).filter(Boolean);
      out.push({ title: p.title, genes });
    } catch {
      out.push({ title: p.title, genes: [] });
    }
  }
  _panelCache = out;
  return out;
}

// Char-by-char text wrapper. jsPDF's built-in splitTextToSize does word-wrap,
// which on Chinese-with-embedded-English text breaks early at the CJK→Latin
// boundary (e.g. "剪接供體位點變異" then linefeed then "(Splice donor)"
// instead of breaking inside the line). Wrapping per glyph eliminates that
// awkwardness and also preserves whitespace (jsPDF was eating spaces around
// CJK boundaries).
function pdfWrapText(doc, text, maxWidth) {
  const lines = [];
  let cur = "";
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (ch === "\n") { lines.push(cur); cur = ""; continue; }
    const tentative = cur + ch;
    if (doc.getTextWidth(tentative) > maxWidth) {
      if (cur) { lines.push(cur); cur = ch; }
      else     { lines.push(ch);  cur = "";  }
    } else {
      cur = tentative;
    }
  }
  if (cur) lines.push(cur);
  return lines;
}

// Minimal layout helper around a jsPDF doc. Tracks a y cursor, auto-
// page-breaks. Bold/Regular variants of NotoSansTC are pre-registered;
// helpers default to Bold so the body text reads "thicker" (the user
// specifically asked for less-thin output).
function makePdfWriter(doc) {
  const pageW = doc.internal.pageSize.getWidth();
  const pageH = doc.internal.pageSize.getHeight();
  const margin = 16;
  const w = { doc, y: margin, pageW, pageH, margin };
  const ensureSpace = need => {
    if (w.y + need > pageH - margin) { doc.addPage(); w.y = margin; }
  };

  w.heading = (text, level = 1) => {
    const sz = level === 1 ? 20 : level === 2 ? 14 : 12;
    if (level !== 1 && w.y > margin + 4) w.y += 4;
    ensureSpace(sz + 4);
    doc.setFont(PDF_FONT_NAME, "bold");
    doc.setFontSize(sz);
    doc.setTextColor(level === 1 ? 30 : 70, level === 1 ? 30 : 70, level === 1 ? 30 : 70);
    doc.text(text, margin, w.y + sz * 0.7);
    w.y += sz * 0.95 + 1;
    if (level === 2) {
      doc.setDrawColor(160, 160, 160);
      doc.setLineWidth(0.4);
      doc.line(margin, w.y, pageW - margin, w.y);
      w.y += 4;
    } else if (level === 1) {
      w.y += 1;
      doc.setDrawColor(60, 60, 60);
      doc.setLineWidth(0.7);
      doc.line(margin, w.y, pageW - margin, w.y);
      w.y += 6;
    }
    doc.setTextColor(20, 20, 20);
  };

  // Body paragraph. Defaults to Bold body — user asked for thicker text.
  // opts.weight = "regular" | "bold"
  w.para = (text, opts = {}) => {
    const sz = opts.size || 10.5;
    const indent = opts.indent || 0;
    const weight = opts.weight === "regular" ? "normal" : "bold";
    doc.setFont(PDF_FONT_NAME, weight);
    doc.setFontSize(sz);
    doc.setTextColor(opts.color || 20, opts.color || 20, opts.color || 20);
    const lineH = sz * 0.55 + 1.8;
    const maxW = pageW - 2 * margin - indent;
    const lines = pdfWrapText(doc, String(text || ""), maxW);
    for (const ln of lines) {
      ensureSpace(lineH + 1);
      doc.text(ln, margin + indent, w.y + lineH * 0.7);
      w.y += lineH;
    }
  };

  w.gap = (h = 4) => { w.y += h; };

  // Two-column key/value row. Label in subtle gray, value bold.
  w.kv = (key, val, opts = {}) => {
    const sz = opts.size || 10.5;
    doc.setFontSize(sz);
    const lineH = sz * 0.55 + 2.2;
    const labelW = opts.labelWidth || 26;
    const maxValW = pageW - 2 * margin - labelW;
    doc.setFont(PDF_FONT_NAME, "bold");
    const lines = pdfWrapText(doc, String(val || ""), maxValW);
    ensureSpace(lineH * lines.length);
    doc.setTextColor(110, 110, 110);
    doc.text(key, margin, w.y + lineH * 0.7);
    doc.setTextColor(20, 20, 20);
    doc.text(lines, margin + labelW, w.y + lineH * 0.7);
    w.y += lineH * lines.length;
  };

  // Subheading inside a section (e.g. variant gene + transcript title).
  w.subheading = (text, opts = {}) => {
    const sz = opts.size || 12;
    ensureSpace(sz + 4);
    doc.setFont(PDF_FONT_NAME, "bold");
    doc.setFontSize(sz);
    doc.setTextColor(30, 30, 30);
    doc.text(text, margin, w.y + sz * 0.7);
    w.y += sz * 0.85 + 2;
  };

  return w;
}

// Per-variant entry in the screening PDF — mirrors the clinical TXT
// report's structure: gene + transcript title, KV block with the table
// columns, two numbered remarks, then the descriptive paragraph.
function pdfWriteVariant(w, vid, v, kind) {
  const h = parseHGVS(v.HGVS);
  const gene = v.gene_symbol || h.gene || "";
  const transcript = h.transcript || "";
  const dis = pickedDiseaseSlot(vid, v);
  const info = diseaseInfo(dis.text);
  const inhTxt = inheritanceCH(info.inheritance);
  const acmgTxt = acmgClassCH(v.ACMG_classification) || "—";

  // Variant title — gene + transcript + HGVS (cdna+protein) all on the
  // subheading row, so the table doesn't need a separate "HGVS" line.
  const titleParts = [gene];
  if (transcript) titleParts.push(transcript);
  const hgvsTxt = hgvsCellText(v);
  if (hgvsTxt) titleParts.push(hgvsTxt);
  w.subheading(titleParts.join(" "));

  // KV info block (the diagnostic-report table cells, just rendered as
  // labelled rows since proportional Chinese kills ASCII alignment). HGVS
  // already lives in the title above, so it's deliberately not repeated.
  w.kv("結構",   v.exon_or_intron || "—");
  w.kv("基因型", ZYG_CH[v.zygosity] || v.zygosity || "—");
  w.kv("ClinVar", variantClinSig(v) || "—");
  w.kv("ACMG",   v.ACMG_classification || "—");
  w.gap(2);

  // The two clinical-style remarks. Second sentence mirrors the
  // 致病性 / 不確定意義 / etc. mapping; recommendation tail is the
  // standard "建議比對臨床表徵" line.
  const tail = kind === "causative"
    ? "與臨床症狀相關"
    : "建議比對臨床表徵";
  w.para(`1. ${gene}為${info.name || "—"}的致病基因之一，其遺傳模式屬於${inhTxt}。`, { indent: 4 });
  w.para(`2. 此為${acmgTxt}之變異位點，${tail}。`, { indent: 4 });
  w.gap(2);

  // Descriptive paragraph — same composition as the clinical TXT block.
  const cons = consequenceEntry(v.Consequence);
  const cdna = h.cdna || "";
  const prot = h.protein ? ` (${h.protein})` : "";
  const para = [
    `在個案之檢體中，檢測到 1 個位於基因 ${gene} 的變異位點。`,
    `變異位點 ${cdna}${prot} 為${cons.label}，${cons.explain}`,
    gnomadAFText(v),
    clinvarText(v),
    acmgGuidelineText(v),
    "此報告僅供參考，臨床判斷仍應以病患的實際狀況為主。建議比對臨床表徵並進行父母親與家族成員之變異位點檢測，以釐清上述變異致病之可能性；根據家族成員變異位點檢測報告或相關資料庫更新，可能影響變異位點 ACMG 判讀結果。",
  ].join("");
  w.para(para, { indent: 4 });
  w.gap(8);
}

// Only the variants the user has marked V on the candidate card (which
// promotes them into the Report area's panel section) end up in the PDF.
// pickedDiseaseSlot() inside pdfWriteVariant already falls back to Disease1
// when no Disease checkbox is ticked.
function pdfWriteSection(w, title, ids, dataVariants, panelKey) {
  w.heading(title, 2);
  const filtered = (ids || []).filter(id =>
    dataVariants?.[id] && _isSecondaryEligible(id) && isSecondarySelected(id, panelKey)
  );
  if (!filtered.length) {
    w.para("（未偵測到致病性之變異位點）", { indent: 4, weight: "regular" });
    w.gap(4);
    return;
  }
  for (const id of filtered) {
    pdfWriteVariant(w, id, dataVariants[id]);
  }
}

// Pharmacogenomics block — paste only the Actionable / Routine summaries,
// no per-gene details (mirrors the toggle the user keeps closed in the UI).
function pdfWritePharmacogenomics(w, pc) {
  w.heading("藥物基因體學", 2);
  if (!pc || !pc.genes || !Object.keys(pc.genes).length) {
    w.para("（無 PharmCAT 結果）", { indent: 4 });
    w.gap(4);
    return;
  }
  const actionable = (pc.actionable || []).map(gene => pc.genes?.[gene]).filter(Boolean);
  const routine = (pc.routine || []).map(gene => pc.genes?.[gene]).filter(Boolean);

  // Actionable
  w.doc.setFont(PDF_FONT_NAME, "bold");
  w.doc.setFontSize(12);
  w.doc.setTextColor(180, 50, 50);
  if (w.y + 14 > w.pageH - w.margin) { w.doc.addPage(); w.y = w.margin; }
  w.doc.text("與用藥相關", w.margin, w.y + 12 * 0.7);
  w.y += 12 * 0.7 + 4;
  w.doc.setTextColor(20, 20, 20);

  for (const g of actionable) {
    w.para(`${g.gene} — ${_pgxGeneSubtitle(g)}`, { size: 10 });
    const drugs = g.drugs || [];
    if (!drugs.length) {
      w.para("（PharmCAT 報告中尚無藥物層級建議）", { size: 9, indent: 6 });
    } else {
      for (const drug of drugs) {
        w.doc.setFont(PDF_FONT_NAME, "normal");
        w.doc.setFontSize(10);
        const dh = 10 * 0.55 + 2;
        if (w.y + dh > w.pageH - w.margin) { w.doc.addPage(); w.y = w.margin; }
        w.doc.setTextColor(70, 70, 70);
        w.doc.text("• " + (drug.drug || "general"), w.margin + 6, w.y + dh * 0.7);
        w.doc.setTextColor(20, 20, 20);
        w.y += dh;
        for (const rec of (drug.recommendations || [])) {
          const txt = String(rec.recommendation || "").replace(/\s+/g, " ").trim();
          if (!txt) continue;
          const label = [rec.source, rec.evidence].filter(Boolean).join(" / ");
          w.para(`(${label || "PGx"}) ${txt}`, { size: 9, indent: 14 });
        }
      }
    }
    w.gap(2);
  }
  if (!actionable.length) w.para("（無）", { size: 9, indent: 6 });
  w.gap(4);

  // Routine.
  w.doc.setFont(PDF_FONT_NAME, "bold");
  w.doc.setFontSize(12);
  w.doc.setTextColor(60, 90, 60);
  if (w.y + 14 > w.pageH - w.margin) { w.doc.addPage(); w.y = w.margin; }
  w.doc.text("標準處方", w.margin, w.y + 12 * 0.7);
  w.y += 12 * 0.7 + 4;
  w.doc.setTextColor(20, 20, 20);

  const groups = new Map();
  for (const g of routine) {
    const sym = String(g.gene || "");
    const label = g.phenotype || "Other";
    const key = label.toLowerCase();
    if (!groups.has(key)) groups.set(key, { label, items: [] });
    groups.get(key).items.push({ gene: sym, subtitle: _pgxGeneSubtitle(g) });
  }
  const sortedKeys = [...groups.keys()].sort();
  for (const key of sortedKeys) {
    const { label, items } = groups.get(key);
    const list = items.map(it => it.subtitle ? `${it.gene} (${it.subtitle})` : it.gene).join(", ");
    w.para(`${label}: ${list}`, { size: 10, indent: 4 });
  }
  if (!groups.size) w.para("（無）", { size: 9, indent: 6 });
  w.gap(4);
}

async function buildScreeningPDF() {
  const fonts  = await loadPdfFonts();
  const panels = await loadGenePanels();

  if (!window.jspdf || !window.jspdf.jsPDF) {
    throw new Error("jsPDF library not loaded");
  }
  const { jsPDF } = window.jspdf;
  const doc = new jsPDF({ unit: "mm", format: "a4" });
  doc.addFileToVFS(`${PDF_FONT_NAME}-Regular.ttf`, fonts.regular);
  doc.addFont(`${PDF_FONT_NAME}-Regular.ttf`, PDF_FONT_NAME, "normal");
  doc.addFileToVFS(`${PDF_FONT_NAME}-Bold.ttf`, fonts.bold);
  doc.addFont(`${PDF_FONT_NAME}-Bold.ttf`, PDF_FONT_NAME, "bold");
  doc.setFont(PDF_FONT_NAME, "normal");

  const data = state.data || {};
  const meta = data.meta || {};
  const cats = data.categories || {};
  const w = makePdfWriter(doc);

  // Title band
  w.heading("全基因組基因篩檢報告", 1);
  w.kv("檢體編號", meta.LIS_ID || state.currentLIS || "");
  w.kv("姓名",     meta.Name || "");
  w.kv("病歷號",   meta.MRN || "");
  w.kv("檢驗項目", meta.Test ? `次世代定序${isWgsTestType(meta.Test) ? "全基因組" : "全外顯子"}定序檢測` : "");
  w.kv("產生日期", todayYmd().replace(/^(\d{4})(\d{2})(\d{2})$/, "$1-$2-$3"));
  w.gap(6);

  pdfWriteSection(w, "重大疾病風險篩檢（美國遺傳醫學會 ACMG 次要發現基因）", cats.acmg_sf,   data.variants, "acmg_sf");
  pdfWriteSection(w, "中風相關基因",                                        cats.stroke, data.variants, "stroke");
  pdfWriteSection(w, "帶因者篩檢",                                          cats.carrier,   data.variants, "carrier");
  pdfWritePharmacogenomics(w, data.pgx || data.pharmcat);

  // 檢測方法說明 — fixed wording for the screening report (assumes WGS,
  // since that's what's used for screening). Numbered list rendered with
  // a hanging indent so wrapped lines align under the text, not the
  // number.
  w.heading("檢測方法說明", 2);
  const methodLines = [
    "1. 本次檢測使用次世代定序儀分析 (Illumina NovaSeq X Plus)。",
    "2. 本次檢測變異位點的錯誤率 ≦ 0.1% (Phred-scaled Q score ≧ 30)。",
    "3. 本次檢測平均定序深度 ≧ 27.5X。",
    "4. 本檢測僅能檢測出基因內單一核苷酸 (single nucleotide)、小片段的缺失或插入 (small indel)，無法檢測出拷貝數變異 (copy number variants)、轉位 (translocation)、倒轉 (inversion) 或其他複雜性結構變異 (complex structural variation)、組織特異性的鑲嵌 (tissue-specific mosaicism)、串聯重複 (tandem repeat) 以及未定序區域 (例如 promoter、intron)。",
    "5. 本檢測報告僅供醫療專業人員參考，需配合其他相關臨床資料與家族成員之相關檢驗。",
    "6. 目前次世代定序分子遺傳診斷皆屬研究性質。",
  ];
  for (const ln of methodLines) w.para(ln, { indent: 4 });
  w.gap(4);

  // 檢測結果注釋 — ClinVar date is dynamic (state.data.clinvar_date),
  // formatted via the same helper the clinical TXT uses; everything else
  // is fixed for the screening report.
  const clinvarDate = ymdToCnDate(data.clinvar_date) || "—";
  w.heading("檢測結果注釋", 2);
  const noteLines = [
    "1. 本檢測結果比對參考序列為人類 hg38 版本。",
    `2. ClinVar 及 ACMG&AMP 指引: 引用 ClinVar 資料庫截至 ${clinvarDate} 更新的註解，及美國醫學遺傳學暨基因體學學會 (ACMG) 與分子病理學學會 (AMP) 2015 年頒佈的指引，並且主要列入致病 (Pathogenic) 及疑似致病 (Likely pathogenic) 變異；其他類別變異經醫師判斷認為與疾病相關時亦可列入。`,
    "3. 參考資料:",
    "     a. 疾病資料庫: OMIM、ClinVar",
    "     b. 族群資料庫: gnomAD (v4.1 genome)",
    "     c. 序列資料庫: RefSeqGene",
  ];
  for (const ln of noteLines) w.para(ln, { indent: 4 });
  w.gap(4);

  // Gene panel listings appended at the end, deliberately small font so
  // they don't dominate the report. Each panel = heading + comma-joined
  // gene list wrapped to page width.
  if (panels.length) {
    w.heading("檢測基因清單", 2);
    for (const p of panels) {
      doc.setFont(PDF_FONT_NAME, "bold");
      doc.setFontSize(10);
      doc.setTextColor(60, 60, 60);
      if (w.y + 12 > w.pageH - w.margin) { doc.addPage(); w.y = w.margin; }
      doc.text(`${p.title}（${p.genes.length} 個基因）`, w.margin, w.y + 7);
      w.y += 9;
      doc.setTextColor(20, 20, 20);
      w.para(p.genes.length ? p.genes.join(", ") : "—",
             { size: 7.5, indent: 2, weight: "regular" });
      w.gap(3);
    }
  }

  return doc;
}

async function exportScreeningReport() {
  if (!state.currentLIS) return;
  const btns = document.querySelectorAll(".btn-export-screening");
  const setBusy = b => btns.forEach(x => { x.disabled = b; });
  const hint = msg => {
    document.querySelectorAll(".js-save-hint").forEach(el => { el.textContent = msg; });
  };
  setBusy(true);
  hint("產生健檢報告 PDF…");
  try {
    const doc = await buildScreeningPDF();
    const meta = state.data?.meta || {};
    const lis = state.currentLIS;
    const mrn = meta.MRN || "MRN";
    const fname = `${lis}_${mrn}_screening_${todayYmd()}.pdf`;
    // Local download first (jsPDF gives us a Blob directly), then push the
    // same bytes to GitHub via the base64 content path.
    downloadBlob(doc.output("blob"), fname);
    const dataUri = doc.output("datauristring");
    const b64 = dataUri.split(",", 2)[1] || "";
    const path = `output/screening_reports/${fname}`;
    await ghPutBinary(path, b64, `export: screening PDF for ${lis}`);
    hint(`已下載 + 匯出 → ${path}`);
  } catch (e) {
    hint("");
    alert("匯出失敗：" + e.message);
  } finally {
    setBusy(false);
  }
}

// ---------- Boot ----------------------------------------------------

async function loadByRow(row) {
  const st = document.getElementById("search-status");
  try {
    st.textContent = `載入中: ${row.LIS_ID} / ${row.Name || "?"}`;
    await loadSample(row.LIS_ID);
    st.textContent = `已載入: ${row.LIS_ID}`;
    renderAll();
    // Samples with > 1 analysis versions get the picker; on confirm
    // the picker reloads + re-renders, so the initial render above
    // is the placeholder until the user picks.
    maybeShowVersionPicker(() => renderAll());
  } catch (e) {
    st.textContent = "錯誤: " + e.message;
  }
}

// Selecting a row from the typeahead drops here via setupCombobox.pick().
// Pressing Enter without a highlighted row tries to resolve whatever
// the user typed against the index; the old [Load] button was redundant
// once dropdown picks already auto-loaded.
async function loadByQuery(q) {
  const st = document.getElementById("search-status");
  st.textContent = "查詢中...";
  let row;
  try { row = resolveLIS(q); }
  catch (e) { st.textContent = e.message; return; }
  if (!row) { st.textContent = "找不到對應樣本（請從下拉選單選擇）"; return; }
  await loadByRow(row);
}

// ---------- LIS_ID combobox typeahead ------------------------------

function setupCombobox() {
  const input = document.getElementById("q-lis");
  const list  = document.getElementById("q-lis-dropdown");
  let activeIdx = -1;

  function renderOptions(rows) {
    list.innerHTML = "";
    activeIdx = -1;
    if (!rows.length) {
      list.classList.add("hidden");
      return;
    }
    rows.forEach((r, i) => {
      const li = document.createElement("li");
      li.className = "combobox-option";
      li.dataset.idx = i;
      li.innerHTML = `
        <span class="opt-lis">${escapeHtml(r.LIS_ID || "")}</span>
        <span class="opt-name">${escapeHtml(maskName(r.Name || ""))}</span>
        <span class="opt-mrn">${escapeHtml(maskMrn(r.MRN || ""))}</span>`;
      li.addEventListener("mousedown", ev => {
        ev.preventDefault();
        pick(r);
      });
      list.appendChild(li);
    });
    list.classList.remove("hidden");
  }

  function currentRows() {
    return matchSamples(input.value);
  }

  async function pick(row) {
    input.value = row.LIS_ID || "";
    list.classList.add("hidden");
    await loadByRow(row);
  }

  input.addEventListener("focus", () => renderOptions(currentRows()));
  input.addEventListener("input", () => renderOptions(currentRows()));
  input.addEventListener("blur", () => {
    // Slight delay to let click land first
    setTimeout(() => list.classList.add("hidden"), 120);
  });
  input.addEventListener("keydown", ev => {
    const opts = Array.from(list.querySelectorAll(".combobox-option"));
    if (ev.key === "ArrowDown") {
      ev.preventDefault();
      activeIdx = Math.min(opts.length - 1, activeIdx + 1);
      opts.forEach((el, i) => el.classList.toggle("active", i === activeIdx));
    } else if (ev.key === "ArrowUp") {
      ev.preventDefault();
      activeIdx = Math.max(0, activeIdx - 1);
      opts.forEach((el, i) => el.classList.toggle("active", i === activeIdx));
    } else if (ev.key === "Enter") {
      ev.preventDefault();
      const rows = currentRows();
      if (activeIdx >= 0 && rows[activeIdx]) {
        pick(rows[activeIdx]);
      } else {
        loadByQuery(input.value);
      }
    } else if (ev.key === "Escape") {
      list.classList.add("hidden");
    }
  });
  document.querySelectorAll(".search-combobox .sample-test-chip input").forEach(filter => {
    filter.addEventListener("change", () => {
      if (filter.checked) sampleTestFilters.add(filter.value);
      else sampleTestFilters.delete(filter.value);
      renderOptions(currentRows());
    });
  });
  document.querySelectorAll(".search-combobox .sample-test-only").forEach(button => {
    button.addEventListener("click", () => {
      const selected = button.dataset.testType || "";
      sampleTestFilters.clear();
      if (selected) sampleTestFilters.add(selected);
      document.querySelectorAll(".search-combobox .sample-test-chip input").forEach(filter => {
        filter.checked = filter.value === selected;
      });
      renderOptions(currentRows());
    });
  });
}

// EMR sync: re-fetch from EMR for the current sample and merge into
// sample_metadata.json server-side, then reload to surface the new
// sex / dob / genetic_counseling. Failures (network, no MRN) bubble
// up as alerts since the button is an explicit reviewer action.
document.getElementById("btn-emr-sync")?.addEventListener("click", async () => {
  if (!state.currentLIS) return;
  if (state.dirty) {
    if (!confirm("有未儲存的編輯，EMR 同步會覆蓋部分欄位（sex / 看診紀錄）。先儲存還是覆蓋？\n\n按取消先去儲存，按確定立即同步。")) {
      return;
    }
  }
  const btn = document.getElementById("btn-emr-sync");
  const orig = btn.textContent;
  btn.disabled = true;
  btn.textContent = "同步中…";
  try {
    const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
    const sid = row?.sample_id || state.currentLIS;
    const result = await apiPost(`/samples/${encodeURIComponent(sid)}/sync_emr`, {});
    await loadSample(state.currentLIS);
    renderAll();
    const cnts = result.changes || {};
    const note = Object.keys(cnts).length
      ? `已同步：${Object.keys(cnts).join(" / ")}`
      : "EMR 無新資料";
    document.getElementById("search-status").textContent = note;
  } catch (e) {
    alert("EMR 同步失敗：" + (e.message || e));
  } finally {
    btn.disabled = false;
    btn.textContent = orig;
  }
});

async function bootAfterAuth() {
  showSampleLoading();
  try {
    await loadIndex();
    const n = state.index ? state.index.length : 0;
    document.getElementById("search-status").textContent =
      n ? `索引已載入：${n} 筆樣本` : "索引為空";
  } catch (e) {
    document.getElementById("search-status").textContent = "載入索引失敗: " + e.message;
  } finally {
    hideSampleLoading();
  }
  // Probe whether the EMR client_id is configured server-side. The
  // 🔄 EMR sync button stays hidden when disabled so the UI doesn't
  // dangle a button that would only ever 503.
  try {
    const probe = await apiFetch("/emr/enabled");
    state.emrEnabled = !!(probe && probe.enabled);
  } catch {
    state.emrEnabled = false;
  }
}

(async function boot() {
  setupCombobox();
  setupDiagnosticAnalysisToggle();
  loadWelcomeVersion();
  updateWelcomeVisibility();
  setupSnvDisplayFilters();
  setupOmimFilter();
  setupHpoSearchInput();
  setupPanelSearchInput();
  setupPhenotypeEvents();

  // Wire login form + logout button.
  document.getElementById("login-form")?.addEventListener("submit", handleLogin);
  document.getElementById("btn-logout")?.addEventListener("click", (ev) => {
    if (ev.currentTarget.dataset.loggedIn === "1") handleLogout();
    else showLoginModal();
  });
  setupPatientListUpload();
  setupCaseList();
  setupGeneSearch();

  // Probe /auth/me; show login modal if no session, otherwise boot the
  // sample index. /auth/me bypasses the global 401 handler because we
  // explicitly catch the failure here.
  try {
    const me = await fetch(`${API_BASE}/auth/me`, { credentials: "same-origin" })
      .then(r => r.ok ? r.json() : null);
    if (!me) { setLoggedInUser(""); showLoginModal(); return; }
    setLoggedInUser(me.username);
    await bootAfterAuth();
  } catch (e) {
    showLoginModal(`啟動失敗：${e.message}`);
  }
})();

// Click + change events for the phenotype editor (delegated so chips
// added by re-render still respond).
function setupPhenotypeEvents() {
  const card = document.getElementById("phenotype-card");
  if (!card) return;

  card.addEventListener("click", ev => {
    const btn = ev.target;
    if (btn.matches(".chip-remove[data-idx]")) {
      ev.stopPropagation();
      removeHpo(Number(btn.dataset.idx));
    } else if (btn.matches(".chip-remove[data-panel-idx]")) {
      ev.stopPropagation();
      removePanel(Number(btn.dataset.panelIdx));
    } else if (btn.matches("#btn-start-analysis")) {
      requestAnalysis();
    }
  });

  card.addEventListener("change", ev => {
    if (ev.target.matches(".chip-weight[data-idx]")) {
      setHpoWeight(Number(ev.target.dataset.idx), ev.target.value);
    } else if (ev.target.matches(".chip-weight[data-panel-idx]")) {
      setPanelWeight(Number(ev.target.dataset.panelIdx), ev.target.value);
    }
  });
}

// SNV/Indel display filters are UI-only within the compact main-screen
// payload. The complete source TSV remains searchable through the modal.
// Re-render candidate tiers so card lists and tab counts stay in sync.
function setupSnvDisplayFilters() {
  for (const id of [
    "filter-disease-associated",
    "filter-in-panel-only",
    "filter-vaf",
    "filter-impact-modifier",
  ]) {
    document.getElementById(id)?.addEventListener("change", () => {
      if (state.data) renderCandidateSections();
    });
  }
}

// OMIM-display toggle: unchecked → hide the .disease-list block under
// every SNV variant card via a class on #category-sections.
function setupOmimFilter() {
  const cb = document.getElementById("filter-omim");
  const modalCb = document.getElementById("gene-search-filter-omim");
  const applyMain = () => document.getElementById("category-sections")
    ?.classList.toggle("hide-omim", !cb?.checked);
  const applyModal = () => document.getElementById("gene-search-modal")
    ?.classList.toggle("hide-omim", !modalCb?.checked);
  applyMain();
  applyModal();
  cb?.addEventListener("change", applyMain);
  modalCb?.addEventListener("change", applyModal);
}

// ---------- tiny utils ---------------------------------------------

// PII masking for the sample-picker dropdown — display only, never
// persisted. Search (matchSamples) still runs against the raw values
// so users can type the real name / MRN to find a patient.
//   maskName("張中民") → "張O民"
//   maskName("李華")   → "李O"
//   maskMrn("12345678") → "12X45X78"
// Array.from() iterates by code point so any surrogate-pair characters
// don't end up half-masked.
function maskName(s) {
  if (!s) return "";
  const chars = Array.from(s);
  if (chars.length < 2) return s;
  chars[1] = "O";
  return chars.join("");
}
function maskMrn(s) {
  if (!s) return "";
  const chars = Array.from(s);
  if (chars.length > 2) chars[2] = "X";
  if (chars.length > 4) chars[4] = "X";
  return chars.join("");
}

function escapeHtml(s) {
  return String(s ?? "").replace(/[&<>"']/g, c => ({
    "&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"
  }[c]));
}
function escapeAttr(s) { return escapeHtml(s); }

// Trigger a browser download for a Blob via a programmatic <a download>.
// Always called from a user-gesture chain (button click → async export
// → here), so popup blockers leave it alone. The actual save path is
// whatever the browser is configured to use (typically ~/Downloads);
// we don't get to pick that.
function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.style.display = "none";
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  // Defer revoke so the download has a chance to start.
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

// ============================================================
// Phase 3: load-new-case + version management UI
// ============================================================

// Generic modal show/hide. Modals carry an `.modal` class and toggle
// `.hidden` to enter/leave the page.
function showModal(id) {
  document.getElementById(id)?.classList.remove("hidden");
}
function hideModal(id) {
  document.getElementById(id)?.classList.add("hidden");
  // IGV needs an explicit browser teardown to stop background fetches.
  if (id === "igv-modal" && typeof closeIgvModal === "function") closeIgvModal();
}

// Wire every modal's close button + outside-click + ESC.
document.addEventListener("click", ev => {
  const t = ev.target;
  const closer = t.closest?.("[data-close]");
  if (closer) {
    hideModal(closer.dataset.close);
    return;
  }
  // Click on the dim backdrop closes the modal too. The form card
  // is the only direct child; everything else (the dim) is the modal
  // itself. Route through hideModal() so per-modal teardown runs.
  if (t.matches?.(".modal")) {
    if (t.id === "sample-loading-modal") return;
    hideModal(t.id);
  }
});
document.addEventListener("keydown", ev => {
  if (ev.key === "Escape") {
    document.querySelectorAll(".modal:not(.hidden):not(#sample-loading-modal)").forEach(m => m.classList.add("hidden"));
  }
});

// ---- Load new case ------------------------------------------------

// The two file inputs each have an "upload | path" tab. Clicking a tab
// hides the other input + clears its value so we don't accidentally
// submit both.
document.addEventListener("click", ev => {
  const tab = ev.target.closest?.(".form-source-tabs .form-tab");
  if (!tab) return;
  const tabs = tab.parentElement;
  const target = tabs.dataset.target;       // "tsv" or "phenotype"
  const mode   = tab.dataset.mode;          // "upload" or "path"
  tabs.querySelectorAll(".form-tab").forEach(b => {
    b.classList.toggle("active", b === tab);
  });
  const form = tab.closest("form");
  const fileInp = form.querySelector(`input[name="${target}_file"]`);
  const pathInp = form.querySelector(`input[name="${target}_path"]`);
  if (mode === "upload") {
    fileInp.hidden = false;
    pathInp.hidden = true;  pathInp.value = "";
  } else {
    pathInp.hidden = false;
    fileInp.hidden = true;  fileInp.value = "";
  }
});

// In-memory map: LIS_ID → entry from /samples/unregistered. Used by
// the dropdown change handler so we don't have to re-fetch the
// preview each time the reviewer scrubs the list.
let _unregisteredById = {};
let _unregisteredList = [];
const UNREGISTERED_CACHE_TTL_MS = 24 * 60 * 60 * 1000;
const _unregisteredCache = {
  loadedAt: 0,
  list: null,
};

// Editable HPO/panel state for the load-new-case modal. Mirrors
// phenoEdit on the analysis page but kept separate so the analysis
// page's running session isn't disturbed while the modal is open.
const newCaseEdit = {
  hpo: [],
  panels: [],
  emrPhenotype: null,   // raw EMR phenotype payload (read-only ref)
  source: "",           // 'reviewer-txt' / 'EMR' / 'edited' — for the source line
};

const NEW_CASE_WGS_VCF_SIZE_BYTES = 100 * 1024 * 1024;

function inferNewCaseTestType(entry) {
  if (!entry) return "";
  const sampleIds = [entry.lis_id, entry.source_sample_id].filter(Boolean);
  if (sampleIds.some(sampleId => normalizeSampleTestType("", sampleId) === "TITAN-WGS")) {
    return "TITAN-WGS";
  }
  const pipelineType = String(entry.pipeline_type || "").toLowerCase();
  const sourcePath = String(entry.source_vcf_path || "").toLowerCase();
  const sourceSample = String(entry.source_sample_id || entry.lis_id || "").toLowerCase();
  if (pipelineType === "dragen" || sourcePath.includes("dragen") || sourceSample.endsWith("-dragen")) {
    return "WGS";
  }
  const sourceSize = Number(entry.source_vcf_size || 0);
  if (sourceSize > NEW_CASE_WGS_VCF_SIZE_BYTES) return "WGS";
  return "";
}

function _formatNewCasePendingLabel(entry) {
  if (!entry) return "";
  const ts = entry.mtime ? new Date(entry.mtime * 1000).toLocaleString() : "";
  const size = entry.tsv_size ? `${(entry.tsv_size / 1024).toFixed(0)} KB` : "";
  return [entry.lis_id, ts, size].filter(Boolean).join("  ·  ");
}

function _newCasePendingSearchText(entry) {
  const roster = entry?.roster || {};
  return [
    entry?.lis_id,
    entry?.source_sample_id,
    entry?.source_vcf_path,
    roster.name,
    roster.mrn,
    roster.department,
    roster.test_type,
  ].filter(Boolean).join(" ").toLowerCase();
}

function _renderNewCaseLisDropdown(query = "", { showAll = false } = {}) {
  const drop = document.getElementById("new-case-lis-id-dropdown");
  if (!drop) return;
  const q = String(query || "").trim().toLowerCase();
  let rows = _unregisteredList || [];
  if (q) rows = rows.filter(r => _newCasePendingSearchText(r).includes(q));
  if (!showAll && !q) {
    drop.classList.add("hidden");
    drop.innerHTML = "";
    return;
  }
  rows = rows.slice(0, 80);
  if (!rows.length) {
    drop.innerHTML = `<li class="combobox-option combobox-empty"><span class="opt-name">（沒有符合的未登錄個案）</span></li>`;
    _comboClearActive(drop);
    drop.classList.remove("hidden");
    return;
  }
  drop.innerHTML = rows.map(r => {
    const roster = r.roster || {};
    const runMeta = [
      r.mtime ? new Date(r.mtime * 1000).toLocaleString() : "",
      r.tsv_size ? `${(r.tsv_size / 1024).toFixed(0)} KB` : "",
    ].filter(Boolean).join(" · ");
    const detail = [
      roster.name || "",
      roster.mrn || "",
      roster.department || "",
    ].filter(Boolean).join(" · ");
    return `<li class="combobox-option" data-new-case-lis-id="${escapeAttr(r.lis_id)}">`
      + `<span class="opt-lis">${escapeHtml(r.lis_id || r.source_sample_id || "")}</span>`
      + `<span class="opt-name">${escapeHtml(runMeta)}</span>`
      + `<span class="opt-mrn">${escapeHtml(detail)}</span>`
      + `</li>`;
  }).join("");
  _comboClearActive(drop);
  drop.classList.remove("hidden");
}

function _clearNewCaseLisSelection() {
  const hidden = document.getElementById("new-case-lis-id");
  if (hidden) hidden.value = "";
  newCaseEdit.hpo = [];
  newCaseEdit.panels = [];
  newCaseEdit.source = "";
  renderNewCasePhenoEditor();
}

function _selectNewCaseLisId(lis_id) {
  const entry = _unregisteredById[lis_id];
  const input = document.getElementById("new-case-lis-id-search");
  const hidden = document.getElementById("new-case-lis-id");
  const drop = document.getElementById("new-case-lis-id-dropdown");
  if (!entry || !input || !hidden) return false;
  input.value = _formatNewCasePendingLabel(entry);
  input.title = [
    entry.lis_id,
    entry.source_sample_id,
    entry.source_vcf_path,
  ].filter(Boolean).join("\n");
  hidden.value = entry.lis_id || "";
  drop?.classList.add("hidden");
  if (drop) _comboClearActive(drop);
  _applyNewCaseLisSelection(entry.lis_id);
  return true;
}

function _pickNewCaseLisOption(opt) {
  const lis_id = opt?.dataset?.newCaseLisId || "";
  if (!lis_id) return false;
  return _selectNewCaseLisId(lis_id);
}

function _setUnregisteredList(list) {
  _unregisteredList = Array.isArray(list) ? list : [];
  _unregisteredById = {};
  _unregisteredList.forEach(r => { if (r?.lis_id) _unregisteredById[r.lis_id] = r; });
}

function _unregisteredCacheFresh() {
  return Array.isArray(_unregisteredCache.list)
    && (Date.now() - Number(_unregisteredCache.loadedAt || 0)) < UNREGISTERED_CACHE_TTL_MS;
}

async function _loadUnregisteredSamples({ force = false } = {}) {
  if (!force && _unregisteredCacheFresh()) {
    _setUnregisteredList(_unregisteredCache.list);
    return _unregisteredList;
  }
  const list = await apiFetch("/samples/unregistered") || [];
  _unregisteredCache.loadedAt = Date.now();
  _unregisteredCache.list = list;
  _setUnregisteredList(list);
  return _unregisteredList;
}

function _removeUnregisteredFromCache(lis_id) {
  if (!lis_id || !Array.isArray(_unregisteredCache.list)) return;
  _unregisteredCache.list = _unregisteredCache.list.filter(r => r?.lis_id !== lis_id);
  _setUnregisteredList(_unregisteredCache.list);
}

function _updateNewCaseLisPlaceholder() {
  const input = document.getElementById("new-case-lis-id-search");
  if (!input) return;
  input.placeholder = _unregisteredList.length
    ? "輸入 LIS ID / sample / 姓名 / MRN 搜尋"
    : "（沒有未登錄的個案）";
}

document.getElementById("btn-new-case")?.addEventListener("click", async () => {
  const form = document.getElementById("new-case-form");
  form?.reset();
  document.getElementById("new-case-error")?.classList.add("hidden");
  newCaseEdit.hpo = [];
  newCaseEdit.panels = [];
  newCaseEdit.emrPhenotype = null;
  newCaseEdit.source = "";
  renderNewCasePhenoEditor();
  renderNewCaseEmrRef();

  // EMR sync button only shows when the server has client_id; mirrors
  // the sample-card behaviour. The probe value is cached on
  // bootAfterAuth so this is just a state read.
  const emrBtn = document.getElementById("btn-new-case-emr");
  if (emrBtn) emrBtn.hidden = !state.emrEnabled;

  // Populate the Category dropdown from /api/options so this modal +
  // the sample-card Category select share one source of truth.
  const catSel = document.getElementById("new-case-category");
  if (catSel) {
    const opts = (state.options && state.options.category_options) || [];
    catSel.innerHTML = `<option value="" selected>—</option>` +
      opts.map(o => `<option value="${escapeAttr(o)}">${escapeHtml(o)}</option>`).join("");
  }

  // Populate the LIS_ID typeahead from a one-day front-end cache. The
  // manual refresh button below can force a rescan when a just-finished
  // pipeline output should appear immediately.
  const lisInput = document.getElementById("new-case-lis-id-search");
  const lisHidden = document.getElementById("new-case-lis-id");
  if (lisInput) {
    lisInput.value = "";
    lisInput.title = "";
    lisInput.placeholder = _unregisteredCacheFresh()
      ? "輸入 LIS ID / sample / 姓名 / MRN 搜尋"
      : "未登錄個案清單載入中…";
  }
  if (lisHidden) lisHidden.value = "";
  _initNewCasePanelTabs();
  _resetNewCasePanelTabs();
  renderNewCaseFixedPanelHosts();
  showModal("new-case-modal");
  try {
    await _loadUnregisteredSamples();
    _updateNewCaseLisPlaceholder();
  } catch (e) {
    if (lisInput) lisInput.placeholder = `讀取失敗：${String(e.message || e)}`;
  }
});

document.getElementById("btn-refresh-unregistered")?.addEventListener("click", async () => {
  const btn = document.getElementById("btn-refresh-unregistered");
  const input = document.getElementById("new-case-lis-id-search");
  const hidden = document.getElementById("new-case-lis-id");
  const oldText = btn?.textContent || "";
  if (btn) {
    btn.disabled = true;
    btn.textContent = "更新中…";
  }
  if (input) input.placeholder = "未登錄個案清單更新中…";
  try {
    await _loadUnregisteredSamples({ force: true });
    if (hidden) hidden.value = "";
    if (input) {
      input.value = "";
      input.title = "";
    }
    _updateNewCaseLisPlaceholder();
    _renderNewCaseLisDropdown("", { showAll: true });
  } catch (e) {
    if (input) input.placeholder = `讀取失敗：${String(e.message || e)}`;
  } finally {
    if (btn) {
      btn.disabled = false;
      btn.textContent = oldText || "更新清單";
    }
  }
});

// Picking a sample preloads phenotype from the reviewer txt + auto-
// fills the MRN that was embedded in the phenotype filename.
function _applyNewCaseLisSelection(lis_id) {
  const entry = _unregisteredById[lis_id];
  if (!entry) {
    _clearNewCaseLisSelection();
    return;
  }
  // Auto-fill MRN / 姓名 / Test type from the uploaded clinic-list
  // roster when this LIS_ID is on it. Fall back to the MRN parsed out
  // of the phenotype.txt filename for samples not yet on any list.
  const roster = entry.roster || null;
  const mrnInput  = document.getElementById("new-case-mrn");
  const nameInput = document.getElementById("new-case-name");
  const testSel   = document.querySelector('#new-case-form select[name="test_type"]');
  const fillMrn = (roster && roster.mrn) || (entry.phenotype && entry.phenotype.mrn) || "";
  if (mrnInput && !mrnInput.value && fillMrn) mrnInput.value = fillMrn;
  if (nameInput && !nameInput.value && roster && roster.name) nameInput.value = roster.name;
  if (testSel) {
    const inferredType = inferNewCaseTestType(entry);
    if (inferredType) testSel.value = inferredType;
    else if (roster && roster.test_type) testSel.value = roster.test_type;
  }
  // Show the ordering department as a hint next to Category — the
  // canonical Category list is in English so we can't auto-pick it
  // from the Chinese 科別, but surfacing it helps the reviewer choose.
  const deptHint = document.getElementById("new-case-dept-hint");
  if (deptHint) deptHint.textContent = (roster && roster.department) ? `科別：${roster.department}` : "";

  if (entry.phenotype && (entry.phenotype.hpo?.length || entry.phenotype.panels?.length)) {
    newCaseEdit.hpo = (entry.phenotype.hpo || []).map(h => ({...h}));
    newCaseEdit.panels = (entry.phenotype.panels || []).map(p => ({...p}));
    newCaseEdit.source = "Web phenotype input tool";
  } else {
    newCaseEdit.hpo = [];
    newCaseEdit.panels = [];
    newCaseEdit.source = "未找到 Web phenotype input tool 紀錄";
  }
  renderNewCasePhenoEditor();
}

// EMR sync button on the modal: pull name / sex / dob / phenotype
// from the EMR APIs and merge into the form. Sex overwrites whatever
// the reviewer picked (per spec). HPO chips get REPLACED with the EMR
// list (so the EMR-reference column below shows what's available;
// reviewer can then edit).
document.getElementById("btn-new-case-emr")?.addEventListener("click", async () => {
  const mrnInput  = document.getElementById("new-case-mrn");
  const nameInput = document.getElementById("new-case-name");
  const sexInput  = document.querySelector('#new-case-form select[name="sex"]');
  const errEl     = document.getElementById("new-case-error");
  const mrn = (mrnInput?.value || "").trim();
  if (!mrn) {
    errEl.textContent = "請先填 MRN 才能 EMR 同步";
    errEl.classList.remove("hidden");
    return;
  }
  errEl.classList.add("hidden");
  const btn = document.getElementById("btn-new-case-emr");
  const orig = btn.textContent;
  btn.disabled = true;
  btn.textContent = "同步中…";
  try {
    const data = await apiFetch(`/emr/${encodeURIComponent(mrn)}`);
    if (!data) throw new Error("EMR 無回應");
    const consult = data.consultation || {};
    const pheno   = data.phenotype    || {};
    if (consult.sex && sexInput)            sexInput.value = consult.sex;
    if (consult.records?.[0] && nameInput && !nameInput.value) {
      // The consultation API doesn't carry the patient's name; nothing
      // to fill from there. Left as-is for the reviewer to type.
    }
    if (pheno.hpo && pheno.hpo.length) {
      // txt phenotype is authoritative: if the reviewer-curated txt
      // had any HPO/panel chips, EMR sync only refreshes the read-only
      // reference row below. Reviewer can manually copy into the
      // editable chips. EMR populates the editable chips only when txt
      // was missing.
      const hasTxt = (newCaseEdit.source || "").startsWith("Web phenotype input tool");
      if (!hasTxt) {
        newCaseEdit.hpo = pheno.hpo.map(h => ({...h}));
        newCaseEdit.source = "EMR phenotype API";
      }
    }
    newCaseEdit.emrPhenotype = pheno;
    renderNewCasePhenoEditor();
    renderNewCaseEmrRef();
  } catch (e) {
    errEl.textContent = "EMR 同步失敗：" + (e.message || e);
    errEl.classList.remove("hidden");
  } finally {
    btn.disabled = false;
    btn.textContent = orig;
  }
});

function renderNewCasePhenoEditor() {
  const hpoUl   = document.getElementById("new-case-hpo-chips");
  const panelUl = document.getElementById("new-case-panel-chips");
  const srcEl   = document.getElementById("new-case-pheno-source");
  if (srcEl) srcEl.textContent = newCaseEdit.source ? `來源：${newCaseEdit.source}` : "";
  if (hpoUl) {
    hpoUl.innerHTML = (newCaseEdit.hpo || []).map((h, i) => {
      const w = Number(h.weight ?? 1);
      const opts = [1,2,3,4,5].map(n => `<option value="${n}" ${n===w?"selected":""}>w=${n}</option>`).join("");
      return `<li class="chip chip-hpo">`
        + `<span class="hpo-id">${escapeHtml(h.phenotype || "")}</span>`
        + `<span class="chip-label">${escapeHtml(h.label || "")}</span>`
        + `<select class="chip-weight" data-nc-hpo-idx="${i}" title="Weight">${opts}</select>`
        + `<button type="button" class="chip-remove" data-nc-hpo-idx="${i}" title="移除">×</button>`
        + `</li>`;
    }).join("");
  }
  if (panelUl) {
    panelUl.innerHTML = (newCaseEdit.panels || []).map((p, i) => {
      const w = Number(p.weight ?? 1);
      const opts = [1,2,3,4,5].map(n => `<option value="${n}" ${n===w?"selected":""}>w=${n}</option>`).join("");
      const name = p.name || "";
      const label = _fixedPanelMeta.has(name) ? _fixedPanelDisplayName(name) : name;
      return `<li class="chip chip-panel">`
        + `<span class="chip-label" title="${escapeAttr(name)}">${escapeHtml(label)}</span>`
        + `<select class="chip-weight" data-nc-panel-idx="${i}" title="Weight">${opts}</select>`
        + `<button type="button" class="chip-remove" data-nc-panel-idx="${i}" title="移除">×</button>`
        + `</li>`;
    }).join("");
  }
  syncNewCaseFixedPanelChipState();
}

function renderNewCaseEmrRef() {
  const host = document.getElementById("new-case-emr-pheno");
  if (!host) return;
  const p = newCaseEdit.emrPhenotype;
  if (!p || !p.found) {
    host.innerHTML = `<span class="muted" style="font-size:12px">尚未從 EMR 抓取（或 EMR 無資料）。</span>`;
    return;
  }
  // Show the EMR text exactly as it lives in the EMR — no chip
  // parsing — so the reviewer can read EMR's own wording and decide
  // what to copy into the editable phenotype chips above.
  const raw = p.raw_content || "";
  host.innerHTML = `
    <textarea class="emr-ref-text" readonly rows="6" placeholder="（EMR 無內容）">${escapeHtml(raw)}</textarea>
    <div class="muted" style="font-size:11px;margin-top:4px">EMR date: ${escapeHtml(p.date || "")}</div>
  `;
}

// Chip remove + weight editing for the modal. Document-level so the
// Stamp "（已編輯）" onto the phenotype source label exactly once,
// so the prefix doesn't accumulate "（已編輯）（已編輯）..." every
// time the reviewer adds/removes a chip.
function _markNewCaseEdited() {
  const tag = "（已編輯）";
  const src = newCaseEdit.source || "";
  if (src.includes(tag)) return;
  newCaseEdit.source = src ? src + tag : "已編輯";
}

// chips can be re-rendered without rebinding listeners.
document.addEventListener("click", ev => {
  const btn = ev.target.closest("[data-nc-hpo-idx], [data-nc-panel-idx]");
  if (!btn || !btn.matches(".chip-remove")) return;
  const hpoIdx = btn.getAttribute("data-nc-hpo-idx");
  const pnlIdx = btn.getAttribute("data-nc-panel-idx");
  if (hpoIdx !== null) newCaseEdit.hpo.splice(Number(hpoIdx), 1);
  if (pnlIdx !== null) newCaseEdit.panels.splice(Number(pnlIdx), 1);
  _markNewCaseEdited();
  renderNewCasePhenoEditor();
});
document.addEventListener("change", ev => {
  const sel = ev.target;
  if (!sel.matches(".chip-weight")) return;
  const hpoIdx = sel.getAttribute("data-nc-hpo-idx");
  const pnlIdx = sel.getAttribute("data-nc-panel-idx");
  if (hpoIdx === null && pnlIdx === null) return;
  const w = Number(sel.value);
  if (!Number.isFinite(w)) return;
  if (hpoIdx !== null && newCaseEdit.hpo[Number(hpoIdx)]) newCaseEdit.hpo[Number(hpoIdx)].weight = w;
  if (pnlIdx !== null && newCaseEdit.panels[Number(pnlIdx)]) newCaseEdit.panels[Number(pnlIdx)].weight = w;
});

// Search dropdowns for the modal. Wire to /api/hpo/search and
// /api/panels with the same shapes the analysis page uses, but
// scoped to #new-case-* element ids so the analysis-page handlers
// don't pick these up.
let _ncHpoSearchTimer = null;
let _ncPanelSearchTimer = null;
document.addEventListener("input", ev => {
  if (ev.target.id === "new-case-hpo-search") {
    clearTimeout(_ncHpoSearchTimer);
    _ncHpoSearchTimer = setTimeout(() => _ncRunHpoSearch(ev.target.value), 200);
  } else if (ev.target.id === "new-case-panel-search") {
    clearTimeout(_ncPanelSearchTimer);
    _ncPanelSearchTimer = setTimeout(() => _ncRunPanelSearch(ev.target.value), 200);
  } else if (ev.target.id === "new-case-lis-id-search") {
    const hidden = document.getElementById("new-case-lis-id");
    if (hidden) hidden.value = "";
    _renderNewCaseLisDropdown(ev.target.value, { showAll: false });
  }
});
async function _ncRunHpoSearch(q) {
  const drop = document.getElementById("new-case-hpo-search-dropdown");
  if (!drop) return;
  q = (q || "").trim();
  if (!q) { drop.classList.add("hidden"); drop.innerHTML = ""; return; }
  let rows = [];
  try { rows = await apiFetch(`/hpo/search?q=${encodeURIComponent(q)}&limit=15`) || []; }
  catch { rows = []; }
  if (!rows.length) { drop.classList.add("hidden"); drop.innerHTML = ""; return; }
  drop.innerHTML = rows.map(r =>
    `<li class="combobox-option" data-nc-hpo-pick='${escapeAttr(JSON.stringify(r))}'>`
    + `<span class="opt-lis">${escapeHtml(r.hpo_id || "")}</span>`
    + `<span class="opt-name">${escapeHtml(r.name || "")}</span>`
    + (r.gene_count ? `<span class="opt-mrn">${r.gene_count} genes</span>` : "")
    + `</li>`
  ).join("");
  _comboClearActive(drop);
  drop.classList.remove("hidden");
}
async function _ncRunPanelSearch(q) {
  const drop = document.getElementById("new-case-panel-search-dropdown");
  if (!drop) return;
  q = (q || "").trim();
  let rows = [];
  try {
    rows = await apiFetch("/panels") || [];
    await loadFixedPanelIndex();
  }
  catch { rows = []; }
  const picked = new Set((newCaseEdit.panels || []).map(p => p.name));
  rows = rows.filter(r => !picked.has(r.name) && !_fixedPanelKeys.has(r.name));
  if (q) {
    const ql = q.toLowerCase();
    rows = rows.filter(r => (r.name || "").toLowerCase().includes(ql));
  }
  rows = rows.slice(0, 15);
  if (!rows.length) { drop.classList.add("hidden"); drop.innerHTML = ""; return; }
  drop.innerHTML = rows.map(r =>
    `<li class="combobox-option" data-nc-panel-pick='${escapeAttr(JSON.stringify(r))}'>`
    + `<span class="opt-name">${escapeHtml(r.name || "")}</span>`
    + `<span class="opt-mrn">${Number(r.gene_count ?? r.n_genes ?? 0)} genes</span>`
    + `</li>`
  ).join("");
  _comboClearActive(drop);
  drop.classList.remove("hidden");
}
document.addEventListener("mousedown", ev => {
  const opt = ev.target.closest("[data-new-case-lis-id], [data-nc-hpo-pick], [data-nc-panel-pick]");
  if (!opt) return;
  ev.preventDefault();
  if (opt.dataset.newCaseLisId) {
    _pickNewCaseLisOption(opt);
  } else if (opt.dataset.ncHpoPick) {
    _pickNewCaseHpoOption(opt);
  } else if (opt.dataset.ncPanelPick) {
    _pickNewCasePanelOption(opt);
  }
});
document.addEventListener("focusin", ev => {
  if (ev.target.id === "new-case-lis-id-search") {
    _renderNewCaseLisDropdown(ev.target.value, { showAll: true });
  }
});
document.addEventListener("focusout", ev => {
  // Slight delay so click on a dropdown row lands first.
  if (ev.target.id === "new-case-lis-id-search") {
    setTimeout(() => document.getElementById("new-case-lis-id-dropdown")?.classList.add("hidden"), 150);
  }
  if (ev.target.id === "new-case-hpo-search") {
    setTimeout(() => document.getElementById("new-case-hpo-search-dropdown")?.classList.add("hidden"), 150);
  }
  if (ev.target.id === "new-case-panel-search") {
    setTimeout(() => document.getElementById("new-case-panel-search-dropdown")?.classList.add("hidden"), 150);
  }
});

document.getElementById("new-case-form")?.addEventListener("submit", async (ev) => {
  ev.preventDefault();
  const form = ev.currentTarget;
  const errEl = document.getElementById("new-case-error");
  errEl.classList.add("hidden");
  errEl.textContent = "";
  const lisHidden = document.getElementById("new-case-lis-id");
  if (!lisHidden?.value) {
    const typed = (document.getElementById("new-case-lis-id-search")?.value || "").trim().toLowerCase();
    const exact = (_unregisteredList || []).find(r => String(r.lis_id || "").toLowerCase() === typed);
    if (exact) _selectNewCaseLisId(exact.lis_id);
  }
  if (!lisHidden?.value) {
    errEl.textContent = "請先從未登錄個案清單選擇一個 LIS ID";
    errEl.classList.remove("hidden");
    document.getElementById("new-case-lis-id-search")?.focus();
    _renderNewCaseLisDropdown(document.getElementById("new-case-lis-id-search")?.value || "", { showAll: true });
    return;
  }
  const fd = new FormData(form);
  // Always send the modal-edited chips; backend uses them as the
  // authoritative phenotype (overrides reviewer txt + EMR fallback).
  fd.set("hpo_json",    JSON.stringify(newCaseEdit.hpo || []));
  fd.set("panels_json", JSON.stringify(newCaseEdit.panels || []));
  showSampleLoading();
  try {
    const resp = await fetch(`${API_BASE}/samples`, {
      method: "POST",
      credentials: "same-origin",
      body: fd,
    });
    if (resp.status === 401) { showLoginModal(); throw new Error("not authenticated"); }
    if (resp.status === 409) {
      const body = await resp.json().catch(() => ({}));
      throw new Error(body.detail || "個案已登錄");
    }
    if (!resp.ok) {
      const body = await resp.json().catch(() => ({}));
      throw new Error(body.detail || `${resp.status} ${resp.statusText}`);
    }
    const out = await resp.json();
    _removeUnregisteredFromCache(lisHidden.value);
    _removeUnregisteredFromCache(out.sample_id);
    hideModal("new-case-modal");
    await loadIndex();
    await loadSample(out.sample_id);
    renderAll();
    const stEl = document.getElementById("search-status");
    if (stEl) stEl.textContent = `已登錄 ${out.sample_id}`;
    if (out.job_id) _startJobPolling(out.sample_id, out.job_id);
  } catch (e) {
    errEl.textContent = e.message;
    errEl.classList.remove("hidden");
  } finally {
    hideSampleLoading();
  }
});

// ---- Re-analyze target picker ------------------------------------

// Replaces the immediate startAnalysis() call. Pops the modal, pre-
// fills the overwrite dropdown with existing versions, and routes the
// confirmed selection back to startAnalysis().
function requestAnalysis() {
  if (!state.currentLIS) return;
  const versions = state.data?.analyses || [];
  const select   = document.getElementById("reanalyze-target");
  const active   = state.data?.active_analysis || "default";

  if (!versions.length) {
    // Brand-new sample with no analysis yet → just go.
    startAnalysis({ version: "default", mode: "overwrite" });
    return;
  }

  select.innerHTML = versions.map(v =>
    `<option value="${escapeAttr(v.name)}" ${v.name === active ? "selected" : ""}>${escapeHtml(v.name)}</option>`
  ).join("");
  document.querySelector('input[name="reanalyze-mode"][value="overwrite"]').checked = true;
  document.getElementById("reanalyze-name").value = "";
  document.getElementById("reanalyze-error")?.classList.add("hidden");
  showModal("reanalyze-modal");
}

document.getElementById("reanalyze-form")?.addEventListener("submit", (ev) => {
  ev.preventDefault();
  const errEl = document.getElementById("reanalyze-error");
  errEl.classList.add("hidden");
  const mode = document.querySelector('input[name="reanalyze-mode"]:checked')?.value;
  let version, runMode;
  if (mode === "new") {
    version = document.getElementById("reanalyze-name").value.trim();
    if (!/^[A-Za-z0-9_\-]{1,32}$/.test(version)) {
      errEl.textContent = "版本名稱必須符合 [A-Za-z0-9_-]{1,32}";
      errEl.classList.remove("hidden");
      return;
    }
    if (version === "default" && (state.data?.analyses || []).some(v => v.name === "default")) {
      errEl.textContent = "default 已存在；改用「覆蓋」或取另一個名稱";
      errEl.classList.remove("hidden");
      return;
    }
    if ((state.data?.analyses || []).some(v => v.name === version)) {
      errEl.textContent = `版本 ${version} 已存在；改用「覆蓋」`;
      errEl.classList.remove("hidden");
      return;
    }
    runMode = "new";
  } else {
    version = document.getElementById("reanalyze-target").value;
    runMode = "overwrite";
  }
  hideModal("reanalyze-modal");
  startAnalysis({ version, mode: runMode });
});

// ---- Version dropdown on the phenotype card ----------------------

function renderVersionPicker() {
  const select = document.getElementById("version-select");
  const delBtn = document.getElementById("btn-delete-version");
  if (!select) return;
  const versions = state.data?.analyses || [];
  const active   = state.data?.active_analysis || "";
  if (!versions.length) {
    select.innerHTML = `<option value="">—</option>`;
    select.disabled = true;
    if (delBtn) delBtn.hidden = true;
    return;
  }
  select.disabled = false;
  select.innerHTML = versions.map(v =>
    `<option value="${escapeAttr(v.name)}" ${v.name === active ? "selected" : ""}>${escapeHtml(v.name)}</option>`
  ).join("");
  if (delBtn) delBtn.hidden = (active === "default");
}

document.getElementById("version-select")?.addEventListener("change", async (ev) => {
  if (!state.currentLIS) return;
  const target = ev.target.value;
  if (!target || target === state.data?.active_analysis) return;

  if (state.dirty) {
    const ok = confirm("有未儲存的編輯。切換版本會丟失它們，繼續？");
    if (!ok) {
      ev.target.value = state.data?.active_analysis || "";
      return;
    }
  }
  const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
  const sid = row?.sample_id || state.currentLIS;
  await apiPut(`/samples/${encodeURIComponent(sid)}/active_analysis`, { name: target });
  await loadSample(state.currentLIS);
  renderAll();
});

document.getElementById("btn-delete-version")?.addEventListener("click", async () => {
  if (!state.currentLIS) return;
  const active = state.data?.active_analysis;
  if (!active || active === "default") return;
  if (!confirm(`刪除版本「${active}」？此操作無法復原。`)) return;
  const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
  const sid = row?.sample_id || state.currentLIS;
  const resp = await fetch(`${API_BASE}/samples/${encodeURIComponent(sid)}/analyses/${encodeURIComponent(active)}`, {
    method: "DELETE",
    credentials: "same-origin",
  });
  if (!resp.ok) {
    alert("刪除失敗：" + resp.statusText);
    return;
  }
  await loadSample(state.currentLIS);
  renderAll();
});

// ---- Multi-version picker on load --------------------------------

// When the loaded sample has more than one analysis version, pop a
// picker so the reviewer chooses which one to land on. Defaults to
// the active version. Single-version samples skip the picker.
function maybeShowVersionPicker(onPick) {
  const versions = state.data?.analyses || [];
  if (versions.length <= 1) return false;
  const active = state.data?.active_analysis || versions[0].name;
  const list = document.getElementById("version-pick-list");
  list.innerHTML = versions.map(v => {
    const meta = [
      v.updated_at ? `updated ${new Date(v.updated_at).toLocaleString()}` : "",
      `${v.n_hpo} HPO + ${v.n_panels} panels`,
      v.note ? `note: ${escapeHtml(v.note)}` : "",
    ].filter(Boolean).join(" · ");
    return `<li data-version="${escapeAttr(v.name)}" class="${v.name === active ? "active" : ""}">
              <span class="v-name">${escapeHtml(v.name)}${v.name === "default" ? " (預設)" : ""}</span>
              <span class="v-meta">${meta}</span>
            </li>`;
  }).join("");
  list.onclick = async (ev) => {
    const li = ev.target.closest("li[data-version]");
    if (!li) return;
    const target = li.dataset.version;
    hideModal("version-pick-modal");
    if (target !== active) {
      const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
      const sid = row?.sample_id || state.currentLIS;
      await apiPut(`/samples/${encodeURIComponent(sid)}/active_analysis`, { name: target });
      await loadSample(state.currentLIS);
    }
    if (onPick) onPick();
  };
  showModal("version-pick-modal");
  return true;
}

// ---------- Lightweight hover tooltip (replaces native `title`) ------
// Native `title` has a long (~1 s) delay and tiny multi-line text; for
// the `ⓘ` hints (mito FILTER / TLOD, …) we use a 0.5 s custom popup.
// Opt in with escaped text on `data-tip` ("\n" → <br>) or trusted,
// application-generated markup on `data-tip-html`.  HTML tooltips stay open
// while hovered so evidence-reference links can be clicked.
(() => {
  let tipEl = null, showTimer = null, hideTimer = null, curTarget = null;
  function ensureEl() {
    if (!tipEl) {
      tipEl = document.createElement("div");
      tipEl.className = "app-tooltip";
      tipEl.setAttribute("role", "tooltip");
      tipEl.style.display = "none";
      tipEl.addEventListener("mouseenter", cancelHide);
      tipEl.addEventListener("mouseleave", scheduleHide);
      document.body.appendChild(tipEl);
    }
    return tipEl;
  }
  function show(el) {
    const txt = el.getAttribute("data-tip");
    const html = el.getAttribute("data-tip-html");
    if (!txt && !html) return;
    const t = ensureEl();
    t.innerHTML = html || String(txt).split("\n").map(escapeHtml).join("<br>");
    t.style.display = "block";
    t.style.left = "0px"; t.style.top = "0px";
    const r = el.getBoundingClientRect();
    const tw = t.offsetWidth, th = t.offsetHeight;
    const vw = document.documentElement.clientWidth;
    const vh = document.documentElement.clientHeight;
    let left = r.left + window.scrollX;
    const maxLeft = window.scrollX + vw - tw - 8;
    if (left > maxLeft) left = Math.max(window.scrollX + 8, maxLeft);
    let top = r.bottom + window.scrollY + 6;
    if (r.bottom + 6 + th > vh) top = r.top + window.scrollY - th - 6;
    t.style.left = left + "px";
    t.style.top = top + "px";
    curTarget = el;
  }
  function cancelHide() {
    if (hideTimer) { clearTimeout(hideTimer); hideTimer = null; }
  }
  function scheduleHide() {
    cancelHide();
    hideTimer = setTimeout(hide, 140);
  }
  function hide() {
    if (showTimer) { clearTimeout(showTimer); showTimer = null; }
    cancelHide();
    curTarget = null;
    if (tipEl) tipEl.style.display = "none";
  }
  document.addEventListener("mouseover", ev => {
    const el = ev.target.closest("[data-tip], [data-tip-html]");
    if (!el || el === curTarget) return;
    cancelHide();
    if (showTimer) clearTimeout(showTimer);
    showTimer = setTimeout(() => { showTimer = null; show(el); }, 350);
  });
  document.addEventListener("mouseout", ev => {
    const el = ev.target.closest("[data-tip], [data-tip-html]");
    if (!el) return;
    if (ev.relatedTarget && el.contains(ev.relatedTarget)) return;
    scheduleHide();
  });
  document.addEventListener("focusin", ev => {
    const el = ev.target.closest("[data-tip], [data-tip-html]");
    if (el) { cancelHide(); show(el); }
  });
  document.addEventListener("focusout", ev => {
    if (ev.target.closest("[data-tip], [data-tip-html]")) scheduleHide();
  });
  document.addEventListener("mousedown", ev => {
    if (tipEl?.contains(ev.target) || ev.target.closest("[data-tip], [data-tip-html]")) return;
    hide();
  }, true);
  window.addEventListener("scroll", hide, true);
})();

// ─────────────────────────────────────────────────────────────────
// 二級分析 — FASTQ samplesheet helper (UI for /api/secondary/*)
// ─────────────────────────────────────────────────────────────────

const _SECONDARY_STATE = {
  index: null,                // {meta, wes: [...], wgs: [...]}
  selected: { wes: "", wgs: "" },
  mode: "",                  // "" | wes | wgs
  batch: [],                 // FASTQ rows to write into samplesheet.csv
};

function _secondaryCurrentList(mode = _SECONDARY_STATE.mode) {
  const idx = _SECONDARY_STATE.index;
  if (!idx) return [];
  return (mode === "wgs" ? idx.wgs : idx.wes) || [];
}

function _secondaryActiveMode() {
  if (_SECONDARY_STATE.selected.wes) return "wes";
  if (_SECONDARY_STATE.selected.wgs) return "wgs";
  return "";
}

function _secondaryRenderMeta() {
  const el = document.getElementById("secondary-index-meta");
  if (!el) return;
  const m = _SECONDARY_STATE.index?.meta;
  if (!m || !m.updated_at) { el.textContent = ""; return; }
  const when = new Date(m.updated_at).toLocaleString();
  const stale = m.stale ? "（已過期，自動背景刷新中）" : "";
  const wgsLanes = Number(m.wgs_lane_count || 0);
  const wgsMeta = wgsLanes ? `WGS ${m.wgs_count} samples / ${wgsLanes} lanes` : `WGS ${m.wgs_count}`;
  el.textContent = `上次更新 ${when} · WES ${m.wes_count} · ${wgsMeta} ${stale}`;
}

function _secondaryFastqLabel(row) {
  const laneCount = Number(row.lane_count || 0);
  const lane = laneCount
    ? ` · ${laneCount} lanes (${Number(row.fastq_file_count || laneCount * 2)} FASTQ)`
    : (row.lane ? ` · ${row.lane}` : "");
  const re = row.reanalysis ? " · reanalysis" : "";
  return `${row.sample_id || ""}${lane}${re} · ${row.run || ""} · ${_dragenFmtSize(row.size)} · ${_dragenFmtMtime(row.mtime)}`;
}

function _secondaryMatchFastqs(mode, query, { showAll = false } = {}) {
  const q = showAll ? "" : (query || "").trim().toLowerCase();
  if (!q && !showAll) return [];
  return _secondaryCurrentList(mode).filter(row => {
    if (!q) return true;
    const laneValues = (row.lanes || []).flatMap(lane => [lane.lane, lane.fastq_1, lane.fastq_2]);
    return [row.sample_id, row.source_sample_id, row.run, row.input_dir, row.fastq_1, row.fastq_2, row.lane, ...laneValues]
      .some(value => String(value || "").toLowerCase().includes(q));
  });
}

function _secondaryFastqTitle(row) {
  if (!Array.isArray(row.lanes) || !row.lanes.length) {
    return `${row.fastq_1 || ""}\n${row.fastq_2 || ""}`;
  }
  return row.lanes.map(lane => `${lane.lane || ""}: ${lane.fastq_1 || ""}\n${lane.lane || ""}: ${lane.fastq_2 || ""}`).join("\n");
}

function _secondaryRenderDropdown(mode, { showAll = false } = {}) {
  const input = document.getElementById(`secondary-fastq-${mode}`);
  const list = document.getElementById(`secondary-fastq-${mode}-dropdown`);
  if (!input || !list) return;
  const rows = _secondaryMatchFastqs(mode, input.value, { showAll });
  list.innerHTML = "";
  if (!rows.length) {
    if (input.value.trim() || showAll) {
      list.innerHTML = `<li class="combobox-option dragen-vcf-option muted">（沒有符合條件的 FASTQ）</li>`;
      list.classList.remove("hidden");
    } else {
      list.classList.add("hidden");
    }
    return;
  }
  rows.forEach(row => {
    const li = document.createElement("li");
    li.className = "combobox-option dragen-vcf-option";
    li.dataset.path = row.fastq_1 || "";
    li.innerHTML = `<span>${escapeHtml(row.sample_id || "")}</span>` +
      `<span class="opt-vcf-meta">${escapeHtml(_secondaryFastqLabel(row))}</span>`;
    li.title = _secondaryFastqTitle(row);
    li.addEventListener("mousedown", ev => {
      ev.preventDefault();
      _secondaryPickFastq(mode, row);
    });
    list.appendChild(li);
  });
  list.classList.remove("hidden");
}

function _secondaryHideDropdown(mode) {
  document.getElementById(`secondary-fastq-${mode}-dropdown`)?.classList.add("hidden");
}

function _secondaryBatchSeqType() {
  return _SECONDARY_STATE.batch[0]?.seq_type || "";
}

function _secondaryCurrentRowFromForm() {
  const mode = _secondaryActiveMode();
  if (!mode) return null;
  const path = _SECONDARY_STATE.selected[mode] || "";
  const row = _secondaryCurrentList(mode).find(v => v.fastq_1 === path);
  if (!row) return null;
  const sid = (document.getElementById("secondary-sample-id")?.value || row.sample_id || "").trim();
  if (!sid) return null;
  return { ...row, sample_id: sid, seq_type: mode.toUpperCase() };
}

function _secondaryRenderBatch() {
  const box = document.getElementById("secondary-batch-list");
  if (!box) return;
  const rows = _SECONDARY_STATE.batch;
  if (!rows.length) {
    box.hidden = true;
    box.innerHTML = "";
    return;
  }
  const uniqueSamples = new Set(rows.map(r => r.sample_id)).size;
  const body = rows.map((row, idx) => {
    const laneCount = Number(row.lane_count || 0);
    const path = laneCount
      ? `${row.input_dir || ""} (${laneCount} lanes / ${Number(row.fastq_file_count || laneCount * 2)} FASTQ)`
      : (row.lane ? `${row.fastq_1 || ""} (${row.lane})` : (row.fastq_1 || ""));
    return `
      <div class="dragen-batch-row">
        <code>${escapeHtml(row.sample_id || "")}</code>
        <code>${escapeHtml(row.seq_type || "")}</code>
        <code>${escapeHtml(row.run || "")}</code>
        <span class="dragen-batch-path" title="${escapeAttr(path)}">${escapeHtml(path)}</span>
        <button type="button" class="dragen-batch-remove" data-idx="${idx}" title="移除">×</button>
      </div>
    `;
  }).join("");
  box.innerHTML = `
    <div class="dragen-batch-head">
      <span>批次清單：${uniqueSamples} 個 sample / ${escapeHtml(rows[0].seq_type || "")}</span>
      <button type="button" class="btn btn-ghost btn-link" id="secondary-batch-clear">清空</button>
    </div>
    ${body}
  `;
  box.hidden = false;
  box.querySelector("#secondary-batch-clear")?.addEventListener("click", () => {
    _SECONDARY_STATE.batch = [];
    _secondaryRenderBatch();
  });
  box.querySelectorAll(".dragen-batch-remove").forEach(btn => {
    btn.addEventListener("click", () => {
      const idx = Number(btn.dataset.idx);
      if (Number.isInteger(idx)) {
        _SECONDARY_STATE.batch.splice(idx, 1);
        _secondaryRenderBatch();
      }
    });
  });
}

function _secondaryAddRow(row) {
  const seqType = row.seq_type || "";
  const batchSeq = _secondaryBatchSeqType();
  if (batchSeq && batchSeq !== seqType) {
    alert("同一批 samplesheet 只能包含 WES 或 WGS 其中一種。");
    return false;
  }
  const dup = _SECONDARY_STATE.batch.find(x => x.fastq_1 === row.fastq_1);
  if (dup) return false;
  _SECONDARY_STATE.batch.push(row);
  return true;
}

function _secondaryAddCurrentToBatch() {
  const row = _secondaryCurrentRowFromForm();
  if (!row) { alert("請先選一個 FASTQ 並確認 Sample ID"); return; }
  if (!_secondaryAddRow(row)) {
    alert("這列 FASTQ 已在批次清單中，或批次類型不一致。");
  }
  _secondaryRenderBatch();
}

function _secondaryAddFolderToBatch() {
  const current = _secondaryCurrentRowFromForm();
  if (!current) { alert("請先選一個 FASTQ"); return; }
  const mode = current.seq_type.toLowerCase();
  const rows = _secondaryCurrentList(mode).filter(row => row.input_dir === current.input_dir);
  let added = 0;
  rows.forEach(row => {
    if (_secondaryAddRow({ ...row, seq_type: current.seq_type })) added += 1;
  });
  if (!added) alert("同資料夾沒有可新增的其他 FASTQ。");
  _secondaryRenderBatch();
}

function _secondaryPickFastq(mode, row) {
  const batchSeq = _secondaryBatchSeqType();
  const seqType = mode.toUpperCase();
  if (batchSeq && batchSeq !== seqType) {
    alert("目前批次已選擇另一種 seq type；請先清空批次再切換。");
    return;
  }
  const otherMode = mode === "wes" ? "wgs" : "wes";
  const input = document.getElementById(`secondary-fastq-${mode}`);
  const other = document.getElementById(`secondary-fastq-${otherMode}`);
  _SECONDARY_STATE.selected[mode] = row.fastq_1 || "";
  _SECONDARY_STATE.selected[otherMode] = "";
  _SECONDARY_STATE.mode = mode;
  if (input) {
    input.value = row.sample_id || "";
    input.title = _secondaryFastqLabel(row);
  }
  if (other) {
    other.value = "";
    other.title = "";
  }
  const sid = document.getElementById("secondary-sample-id");
  if (sid) sid.value = row.sample_id || "";
  _secondaryHideDropdown(mode);
  _secondaryHideDropdown(otherMode);
}

async function loadSecondaryFastqList({ force = false } = {}) {
  const wes = document.getElementById("secondary-fastq-wes");
  const wgs = document.getElementById("secondary-fastq-wgs");
  if (!wes || !wgs) return;
  if (!force && _SECONDARY_STATE.index) {
    _secondaryRenderMeta();
    return;
  }
  wes.placeholder = "FASTQ 清單載入中…";
  wgs.placeholder = "FASTQ 清單載入中…";
  try {
    const idx = await apiFetch(force ? "/secondary/index/refresh" : "/secondary/fastqs", {
      method: force ? "POST" : "GET",
    });
    _SECONDARY_STATE.index = idx;
    _secondaryRenderMeta();
  } catch (e) {
    wes.placeholder = `載入失敗：${String(e)}`;
    wgs.placeholder = `載入失敗：${String(e)}`;
    return;
  }
  wes.placeholder = "輸入 sample / run / path 搜尋";
  wgs.placeholder = "輸入 sample / run / path 搜尋";
}

async function _secondaryCreateSamplesheet() {
  let samples = _SECONDARY_STATE.batch.slice();
  if (!samples.length) {
    const row = _secondaryCurrentRowFromForm();
    if (!row) { alert("請先挑選 FASTQ，或先加入批次"); return; }
    samples = [row];
  }
  const seqType = samples[0]?.seq_type || "";
  const batchName = (document.getElementById("secondary-batch-name")?.value || "").trim();
  const btn = document.getElementById("secondary-create-btn");
  if (btn) { btn.disabled = true; btn.textContent = "建立中…"; }
  try {
    const result = await apiFetch("/secondary/samplesheet", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ seq_type: seqType, batch_name: batchName, samples }),
    });
    _SECONDARY_STATE.batch = samples.slice();
    _secondaryRenderBatch();
    _secondaryRenderResult(result);
  } catch (e) {
    alert("建立 samplesheet 失敗：" + (e?.message || e));
  } finally {
    if (btn) { btn.disabled = false; btn.textContent = "建立 sample sheet"; }
  }
}

function _secondaryRenderResult(result) {
  const panel = document.getElementById("secondary-result-panel");
  const meta = document.getElementById("secondary-result-meta");
  const command = document.getElementById("secondary-command");
  const attach = document.getElementById("secondary-help-attach");
  if (!panel || !result) return;
  panel.hidden = false;
  if (meta) {
    const warn = (result.warnings || []).length ? `；${result.warnings.join("；")}` : "";
    meta.textContent = `${result.samplesheet_path || ""} → ${result.dgx_output_dir || ""}${warn}`;
  }
  if (command) command.textContent = result.command || "";
  if (attach) attach.textContent = `tmux attach -t ${result.tmux_session || ""}`;
}

async function _secondaryCopyCommand() {
  const text = document.getElementById("secondary-command")?.textContent || "";
  if (!text) return;
  const btn = document.getElementById("secondary-copy-command");
  const markCopied = () => {
    if (!btn) return;
    const old = btn.textContent;
    btn.textContent = "已複製";
    setTimeout(() => { btn.textContent = old; }, 1200);
  };
  try {
    if (navigator.clipboard?.writeText && window.isSecureContext) {
      await navigator.clipboard.writeText(text);
      markCopied();
      return;
    }
    const ta = document.createElement("textarea");
    ta.value = text;
    ta.setAttribute("readonly", "");
    ta.style.position = "fixed";
    ta.style.left = "-9999px";
    ta.style.top = "0";
    document.body.appendChild(ta);
    ta.focus();
    ta.select();
    const ok = document.execCommand("copy");
    document.body.removeChild(ta);
    if (!ok) throw new Error("execCommand copy returned false");
    markCopied();
  } catch (_e) {
    alert("無法使用瀏覽器剪貼簿，請手動選取指令複製。");
  }
}

function _secondaryWireCombobox(mode) {
  const input = document.getElementById(`secondary-fastq-${mode}`);
  const list = document.getElementById(`secondary-fastq-${mode}-dropdown`);
  if (!input || !list) return;
  let activeIdx = -1;
  input.addEventListener("focus", () => {
    activeIdx = -1;
    _secondaryRenderDropdown(mode, { showAll: true });
  });
  input.addEventListener("input", () => {
    activeIdx = -1;
    _SECONDARY_STATE.selected[mode] = "";
    input.title = "";
    _SECONDARY_STATE.mode = _secondaryActiveMode();
    _secondaryRenderDropdown(mode);
  });
  input.addEventListener("blur", () => {
    setTimeout(() => _secondaryHideDropdown(mode), 120);
  });
  input.addEventListener("keydown", ev => {
    const opts = Array.from(list.querySelectorAll(".dragen-vcf-option[data-path]"));
    if (ev.key === "ArrowDown") {
      ev.preventDefault();
      activeIdx = Math.min(opts.length - 1, activeIdx + 1);
      opts.forEach((el, i) => el.classList.toggle("active", i === activeIdx));
    } else if (ev.key === "ArrowUp") {
      ev.preventDefault();
      activeIdx = Math.max(0, activeIdx - 1);
      opts.forEach((el, i) => el.classList.toggle("active", i === activeIdx));
    } else if (ev.key === "Enter" && activeIdx >= 0 && opts[activeIdx]) {
      ev.preventDefault();
      const row = _secondaryCurrentList(mode).find(v => v.fastq_1 === opts[activeIdx].dataset.path);
      if (row) _secondaryPickFastq(mode, row);
    } else if (ev.key === "Escape") {
      _secondaryHideDropdown(mode);
    }
  });
}

function setupSecondaryButton() {
  let btn = document.getElementById("btn-secondary-launch");
  if (!btn) {
    const dragenBtn = document.getElementById("btn-dragen-launch");
    if (!dragenBtn?.parentElement) return;
    btn = document.createElement("button");
    btn.id = "btn-secondary-launch";
    btn.className = "btn btn-ghost";
    btn.type = "button";
    btn.hidden = dragenBtn.hidden;
    btn.title = "搜尋 FASTQ、建立二級分析 samplesheet，並產生 DGX2 tmux 執行指令";
    btn.textContent = "二級分析";
    dragenBtn.parentElement.insertBefore(btn, dragenBtn);
  }
  if (!btn) return;
  btn.addEventListener("click", async () => {
    showModal("secondary-modal");
    await loadSecondaryFastqList();
  });
  _secondaryWireCombobox("wes");
  _secondaryWireCombobox("wgs");
  document.getElementById("secondary-refresh-btn")?.addEventListener("click", async ev => {
    const b = ev.currentTarget;
    if (b) { b.disabled = true; b.textContent = "↻ 更新中…"; }
    try { await loadSecondaryFastqList({ force: true }); }
    finally { if (b) { b.disabled = false; b.textContent = "↻ 更新索引"; } }
  });
  document.getElementById("secondary-show-clean-command")?.addEventListener("click", _secondaryShowCleanupCommand);
  document.getElementById("secondary-add-batch-btn")?.addEventListener("click", _secondaryAddCurrentToBatch);
  document.getElementById("secondary-add-folder-btn")?.addEventListener("click", _secondaryAddFolderToBatch);
  document.getElementById("secondary-create-btn")?.addEventListener("click", _secondaryCreateSamplesheet);
  document.getElementById("secondary-copy-command")?.addEventListener("click", _secondaryCopyCommand);
  document.getElementById("secondary-copy-clean-command")?.addEventListener("click", _secondaryCopyCleanupCommand);
}

async function _secondaryShowCleanupCommand() {
  const btn = document.getElementById("secondary-show-clean-command");
  const panel = document.getElementById("secondary-clean-result-panel");
  const commandEl = document.getElementById("secondary-clean-command");
  const meta = document.getElementById("secondary-clean-result-meta");
  if (btn) { btn.disabled = true; btn.textContent = "產生中…"; }
  try {
    const result = await apiFetch("/secondary/nf-work/cleanup-command");
    if (!result?.command) throw new Error("後端未回傳 DGX2 清理指令");
    if (commandEl) commandEl.textContent = result.command;
    if (meta) meta.textContent = `請在 DGX2 執行 · ${result?.path || "/raid/DGM/work"}`;
    if (panel) panel.hidden = false;
  } catch (e) {
    if (commandEl) commandEl.textContent = `產生清理指令失敗：${e.message || e}`;
    if (meta) meta.textContent = "";
    if (panel) panel.hidden = false;
  } finally {
    if (btn) { btn.disabled = false; btn.textContent = "顯示 DGX2 清理指令"; }
  }
}

async function _secondaryCopyCleanupCommand() {
  const command = document.getElementById("secondary-clean-command")?.textContent || "";
  if (!command || command.startsWith("產生清理指令失敗")) return;
  try {
    await navigator.clipboard.writeText(command);
    const btn = document.getElementById("secondary-copy-clean-command");
    if (btn) {
      btn.textContent = "已複製";
      setTimeout(() => { btn.textContent = "複製"; }, 1200);
    }
  } catch (_e) {
    prompt("請複製以下 DGX2 清理指令：", command);
  }
}

// 三級分析 — DRAGEN VCF kicker (UI for /api/dragen/*)
//
// Button on the topbar opens a modal listing every hard-filtered
// VCF found under the server's DRAGEN_VCF_ROOTS, reviewer picks
// one and clicks 開始分析. Backend spawns a worker process
// (samplesheet → nextflow → copy outputs → post-processing). We poll
// /api/dragen/jobs/{job_id} every 5 s and reflect the current
// step in (a) the modal log pane and (b) a small grey status
// label next to the "成大醫院基因醫學部 NGS 分析平台" title.
// ─────────────────────────────────────────────────────────────────

const _DRAGEN_STATE = {
  // mode is derived from the combobox that has a picked VCF path.
  // After picking, the input keeps only the source sample ID so reviewers
  // can edit VAL-36 → VAL-37 without deleting run/size/date metadata.
  // selected keeps the actual path used when starting the worker.
  mode: "",             // ""    | dragen | inhouse
  index: null,          // {meta, dragen: [...], inhouse: [...]}
  selected: { inhouse: "", dragen: "" },
  seqType: { inhouse: "" },
  batch: [],            // [{mode, vcf_path, sample_id, source_sample_id, ...}]
  job: null,            // current job state, polled
  lastProgressPct: 0,   // keep progress visually monotonic across polling races
  pollTimer: null,
  recoverTimer: null,
};

const DRAGEN_INHOUSE_WGS_SIZE_THRESHOLD = 100 * 1024 * 1024;

function _dragenFmtSize(b) {
  if (!b) return "—";
  const u = ["B","KB","MB","GB","TB"];
  let i = 0, x = b;
  while (x >= 1024 && i < u.length - 1) { x /= 1024; i++; }
  return `${x.toFixed(x < 10 && i > 0 ? 1 : 0)} ${u[i]}`;
}

function _dragenFmtMtime(t) {
  if (!t) return "";
  try {
    return new Intl.DateTimeFormat("zh-TW", {
      timeZone: "Asia/Taipei",
      year: "numeric", month: "2-digit", day: "2-digit",
      hour: "2-digit", minute: "2-digit", hourCycle: "h23",
    }).format(new Date(t * 1000));
  }
  catch { return ""; }
}

function _dragenInferSeqType(row) {
  const fromIndex = String(row?.seq_type || "").trim().toUpperCase();
  if (fromIndex === "WES" || fromIndex === "WGS") return fromIndex;
  return Number(row?.size || 0) >= DRAGEN_INHOUSE_WGS_SIZE_THRESHOLD ? "WGS" : "WES";
}

function _dragenCurrentList(mode = _DRAGEN_STATE.mode) {
  const idx = _DRAGEN_STATE.index;
  if (!idx) return [];
  return (mode === "inhouse" ? idx.inhouse : idx.dragen) || [];
}

// Which of the two comboboxes currently carries a picked VCF.
function _dragenActiveMode() {
  if (_DRAGEN_STATE.selected.inhouse) return "inhouse";
  if (_DRAGEN_STATE.selected.dragen)  return "dragen";
  return "";
}

function _dragenRenderMeta() {
  const el = document.getElementById("dragen-index-meta");
  if (!el) return;
  const m = _DRAGEN_STATE.index?.meta;
  if (!m || !m.updated_at) { el.textContent = ""; return; }
  const when = new Date(m.updated_at).toLocaleString();
  const counts = `dragen: ${m.dragen_count} · in-house: ${m.inhouse_count}`;
  const stale = m.stale ? "（已過期，自動背景刷新中）" : "";
  el.textContent = `上次更新 ${when}  ·  ${counts}  ${stale}`;
}

function _dragenRenderSiblings() {
  const box = document.getElementById("dragen-siblings");
  const body = document.getElementById("dragen-siblings-body");
  if (!box || !body) return;
  box.hidden = true;
  body.innerHTML = "";
}

function _dragenRenderSeqType() {
  const rowEl = document.getElementById("dragen-seq-type-row");
  const meta = document.getElementById("dragen-seq-type-meta");
  if (!rowEl) return;
  const mode = _dragenActiveMode();
  const path = _DRAGEN_STATE.selected.inhouse;
  if (mode !== "inhouse" || !path) {
    rowEl.hidden = true;
    if (meta) meta.textContent = "";
    return;
  }
  const row = _dragenCurrentList("inhouse").find(v => v.path === path);
  const seqType = _DRAGEN_STATE.seqType.inhouse || _dragenInferSeqType(row || {});
  _DRAGEN_STATE.seqType.inhouse = seqType;
  rowEl.hidden = false;
  rowEl.querySelectorAll(".dragen-seq-option").forEach(btn => {
    btn.classList.toggle("active", btn.dataset.seqType === seqType);
    btn.setAttribute("aria-pressed", btn.dataset.seqType === seqType ? "true" : "false");
  });
  if (meta) meta.textContent = "";
}

function _dragenVcfLabel(v) {
  return `${v.sample_id} · ${v.run || ""} · ${_dragenFmtSize(v.size)} · ${_dragenFmtMtime(v.mtime)}`;
}

function _dragenMatchVcfs(mode, query, { showAll = false } = {}) {
  const q = showAll ? "" : (query || "").trim().toLowerCase();
  if (!q && !showAll) return [];
  return _dragenCurrentList(mode).filter(v => {
    if (!q) return true;
    return [v.sample_id, v.run, v.path, v.cnv_vcf, v.sv_vcf, v.mito_vcf]
      .some(value => String(value || "").toLowerCase().includes(q));
  });
}

function _dragenRenderDropdown(mode, { showAll = false } = {}) {
  const input = document.getElementById(`dragen-vcf-${mode}`);
  const list  = document.getElementById(`dragen-vcf-${mode}-dropdown`);
  if (!input || !list) return;
  const rows = _dragenMatchVcfs(mode, input.value, { showAll });
  list.innerHTML = "";
  if (!rows.length) {
    if (input.value.trim() || showAll) {
      list.innerHTML = `<li class="combobox-option dragen-vcf-option muted">（沒有符合條件的 VCF）</li>`;
      list.classList.remove("hidden");
    } else {
      list.classList.add("hidden");
    }
    return;
  }
  rows.forEach(row => {
    const li = document.createElement("li");
    li.className = "combobox-option dragen-vcf-option";
    li.dataset.path = row.path || "";
    li.innerHTML = `<span>${escapeHtml(row.sample_id || "")}</span>` +
      `<span class="opt-vcf-meta">${escapeHtml(row.run || "")} · ${escapeHtml(_dragenFmtSize(row.size))} · ${escapeHtml(_dragenFmtMtime(row.mtime))}</span>`;
    li.title = row.path || "";
    li.addEventListener("mousedown", ev => {
      ev.preventDefault();
      _dragenPickVcf(mode, row);
    });
    list.appendChild(li);
  });
  list.classList.remove("hidden");
}

function _dragenHideDropdown(mode) {
  document.getElementById(`dragen-vcf-${mode}-dropdown`)?.classList.add("hidden");
}

function _dragenBatchMode() {
  return _DRAGEN_STATE.batch[0]?.mode || "";
}

function _dragenBatchLocked() {
  const state = _DRAGEN_STATE.job?.state || "";
  return state === "queued" || state === "running" || !!_DRAGEN_STATE.job?.running;
}

function _dragenBuildSample(mode, row, sampleId) {
  const sample = {
    mode,
    vcf_path: row.path || "",
    sample_id: sampleId,
    source_sample_id: row.sample_id || "",
    seq_type: mode === "inhouse" ? (_DRAGEN_STATE.seqType.inhouse || _dragenInferSeqType(row)) : "WGS",
  };
  if (mode === "inhouse") {
    sample.cnv_vcf = row.cnv_vcf || "";
    sample.sv_vcf = row.sv_vcf || "";
    sample.mito_vcf = row.mito_vcf || "";
  }
  return sample;
}

function _dragenCurrentSampleFromForm() {
  const mode = _dragenActiveMode();
  if (!mode) return null;
  const path = _DRAGEN_STATE.selected[mode] || "";
  if (!path) return null;
  const row = _dragenCurrentList(mode).find(v => v.path === path);
  if (!row) return null;
  const sidIn = document.getElementById("dragen-sample-id");
  const sampleId = (sidIn?.value || row.sample_id || "").trim();
  if (!sampleId) return null;
  return _dragenBuildSample(mode, row, sampleId);
}

function _dragenRenderBatch() {
  const box = document.getElementById("dragen-batch-list");
  if (!box) return;
  const rows = _DRAGEN_STATE.batch;
  if (!rows.length) {
    box.hidden = true;
    box.innerHTML = "";
    return;
  }
  const modeLabel = rows[0].mode === "inhouse" ? "in-house" : "DRAGEN";
  const locked = _dragenBatchLocked();
  const body = rows.map((row, idx) => `
    <div class="dragen-batch-row">
      <code>${escapeHtml(row.sample_id || "")}</code>
      <code>${escapeHtml(row.source_sample_id || "")}</code>
      <code>${escapeHtml(row.seq_type || "")}</code>
      <span class="dragen-batch-path" title="${escapeAttr(row.vcf_path || "")}">${escapeHtml(row.vcf_path || "")}</span>
      <button type="button" class="dragen-batch-remove" data-idx="${idx}" title="移除" ${locked ? "disabled" : ""}>×</button>
    </div>
  `).join("");
  box.innerHTML = `
    <div class="dragen-batch-head">
      <span>批次清單：${rows.length} 個 ${modeLabel} sample</span>
      <button type="button" class="btn btn-ghost btn-link" id="dragen-batch-clear" ${locked ? "disabled" : ""}>清空</button>
    </div>
    ${body}
  `;
  box.hidden = false;
  box.querySelector("#dragen-batch-clear")?.addEventListener("click", () => {
    if (_dragenBatchLocked()) return;
    _DRAGEN_STATE.batch = [];
    _dragenRenderBatch();
  });
  box.querySelectorAll(".dragen-batch-remove").forEach(btn => {
    btn.addEventListener("click", () => {
      const idx = Number(btn.dataset.idx);
      if (Number.isInteger(idx)) {
        if (_dragenBatchLocked()) return;
        _DRAGEN_STATE.batch.splice(idx, 1);
        _dragenRenderBatch();
      }
    });
  });
}

function _dragenAddCurrentToBatch() {
  if (_dragenBatchLocked()) {
    alert("三級分析執行中，請等這批完成後再修改批次清單。");
    return;
  }
  const sample = _dragenCurrentSampleFromForm();
  if (!sample) { alert("請先選一個 VCF 並確認 Sample ID"); return; }
  const batchMode = _dragenBatchMode();
  if (batchMode && batchMode !== sample.mode) {
    alert("同一批 sample sheet 只能包含同一種來源：請先清空批次或分開執行。");
    return;
  }
  const dup = _DRAGEN_STATE.batch.find(row =>
    row.sample_id === sample.sample_id || row.source_sample_id === sample.source_sample_id || row.vcf_path === sample.vcf_path
  );
  if (dup) { alert("這個 sample 已在批次清單中"); return; }
  _DRAGEN_STATE.batch.push(sample);
  _dragenRenderBatch();
}

function _dragenPickVcf(mode, row) {
  if (_dragenBatchLocked()) {
    alert("三級分析執行中，請等這批完成後再修改批次清單。");
    return;
  }
  const batchMode = _dragenBatchMode();
  if (batchMode && batchMode !== mode) {
    alert("目前批次已選擇另一種來源；請先清空批次再切換來源。");
    return;
  }
  const otherMode = mode === "inhouse" ? "dragen" : "inhouse";
  const input = document.getElementById(`dragen-vcf-${mode}`);
  const other = document.getElementById(`dragen-vcf-${otherMode}`);
  _DRAGEN_STATE.selected[mode] = row.path || "";
  _DRAGEN_STATE.selected[otherMode] = "";
  _DRAGEN_STATE.mode = mode;
  if (mode === "inhouse") _DRAGEN_STATE.seqType.inhouse = _dragenInferSeqType(row);
  else _DRAGEN_STATE.seqType.inhouse = "";
  if (input) {
    input.value = row.sample_id || "";
    input.title = _dragenVcfLabel(row);
  }
  if (other) {
    other.value = "";
    other.title = "";
  }
  _dragenHideDropdown(mode);
  _dragenHideDropdown(otherMode);
  const sidIn = document.getElementById("dragen-sample-id");
  if (sidIn && row.sample_id) sidIn.value = _dragenSuggestSid(row.sample_id, mode);
  _dragenRenderSeqType();
  _dragenRenderSiblings();
}

async function loadDragenVcfList({ force = false } = {}) {
  const inhouseInput = document.getElementById("dragen-vcf-inhouse");
  const dragenInput  = document.getElementById("dragen-vcf-dragen");
  if (!inhouseInput || !dragenInput) return;
  if (!force && _DRAGEN_STATE.index) {
    _dragenRenderMeta(); return;
  }
  inhouseInput.placeholder = "VCF 清單載入中…";
  dragenInput.placeholder  = "VCF 清單載入中…";
  try {
    const idx = await apiFetch(force ? "/dragen/index/refresh" : "/dragen/vcfs", {
      method: force ? "POST" : "GET",
    });
    _DRAGEN_STATE.index = idx;
    _dragenRenderMeta();
  } catch (e) {
    inhouseInput.placeholder = `載入失敗：${String(e)}`;
    dragenInput.placeholder  = `載入失敗：${String(e)}`;
    return;
  }
  inhouseInput.placeholder = "輸入 sample / run / path 搜尋";
  dragenInput.placeholder  = "輸入 sample / run / path 搜尋";
}

function _dragenSetJob(state) {
  const isNewJob = state?.job_id && state.job_id !== _DRAGEN_STATE.job?.job_id;
  _DRAGEN_STATE.job = state;
  if (!state || isNewJob) _DRAGEN_STATE.lastProgressPct = 0;
  const rawPct = _dragenProgressPercent(state);
  const pct = state?.state === "done"
    ? 100
    : Math.max(rawPct, _DRAGEN_STATE.lastProgressPct || 0);
  _DRAGEN_STATE.lastProgressPct = pct;
  // Topbar status
  const top = document.getElementById("topbar-job-status");
  if (top) {
    if (!state || state.state === "done" || state.state === "failed" || state.state === "cancelled") {
      top.hidden = true;
      top.textContent = "";
    } else {
      top.hidden = false;
      const mode = state.mode === "inhouse" ? "in-house" : "dragen";
      const label = state.sample_count && state.sample_count > 1
        ? `${state.sample_count} samples`
        : state.sample_id;
      top.textContent = `· 三級分析 [${mode}]: ${label} (${pct}%)`;
    }
  }
  // Modal panel
  const panel = document.getElementById("dragen-job-panel");
  const stepEl  = document.getElementById("dragen-job-step");
  const stateEl = document.getElementById("dragen-job-state");
  const logEl   = document.getElementById("dragen-job-log");
  const logToggle = document.getElementById("dragen-job-log-toggle");
  const progressBar = document.getElementById("dragen-job-progress-bar");
  const progressText = document.getElementById("dragen-job-progress-text");
  const cancelBtn = document.getElementById("dragen-job-cancel-btn");
  if (panel && state) {
    panel.hidden = false;
    if (stepEl)  stepEl.textContent  = state.step || "";
    if (stateEl) stateEl.textContent = state.error
      ? `failed — ${state.error}`
      : state.state;
    if (progressBar) progressBar.style.width = `${pct}%`;
    if (progressText) progressText.textContent = `${pct}%`;
    if (progressBar) progressBar.classList.toggle("failed", state.state === "failed");
    if (cancelBtn) {
      const running = state.state === "queued" || state.state === "running" || !!state.running;
      cancelBtn.hidden = !running;
      cancelBtn.disabled = !running;
    }
    if (isNewJob && logEl) logEl.hidden = true;
    if (isNewJob && logToggle) {
      logToggle.textContent = "▶ Log";
      logToggle.setAttribute("aria-expanded", "false");
    }
    if (logEl)   logEl.textContent   = state.log_tail || "";
    if (logEl)   logEl.scrollTop     = logEl.scrollHeight;
    _dragenRenderBatch();
  }
}

function _dragenProgressPercent(state) {
  if (!state) return 0;
  if (state.state === "done") return 100;
  if (state.state === "cancelled") return _clampPct(_dragenProgressPercent({ ...state, state: "running" }));
  if (state.state === "failed") return _clampPct(_dragenProgressPercent({ ...state, state: "running" }));
  const step = String(state.step || "");
  if (step.startsWith("post-processing") || step.startsWith("sample-step")) {
    return _dragenStopgapProgressPercent(state);
  }
  if ((step === "nextflow" || step.startsWith("nextflow:")) && Number.isFinite(Number(state.nextflow_progress_pct))) {
    return _clampPct(state.nextflow_progress_pct);
  }
  const byStep = {
    queued: 1,
    "detect-pipeline-output": 2,
    samplesheet: 2,
    stage: 2,
    nextflow: 3,
    "nextflow:prepare-vcf": 4,
    "nextflow:prepare-vcf-dragen": 4,
    "nextflow:add-callers-tag": 4,
    "nextflow:filter-for-annotation": 6,
    "nextflow:vep-annotate": 20,
    "nextflow:pangolin-score": 52,
    "nextflow:parse-csq": 60,
    "nextflow:acmg-classify": 76,
    "nextflow:vep-annotate-done": 52,
    "copy-pipeline-tsv": 82,
    "prepare-postprocessing": 82,
    "post-processing": 82,
    done: 100,
  };
  return byStep[state.step] ?? 0;
}

function _clampPct(value) {
  return Math.max(0, Math.min(100, Math.round(Number(value) || 0)));
}

function _dragenStopgapProgressPercent(state) {
  const total = Math.max(1, Number(state.post_processing_sample_count || state.stopgap_sample_count || state.sample_count || 1));
  const idx = Math.max(0, Math.min(total - 1, Number(state.post_processing_sample_index ?? state.stopgap_sample_index ?? 0)));
  const sub = {
    "post-processing": 0,
    "post-processing:genebe": 0.20,
    "post-processing:extra-vep": 0.52,
    "post-processing:giab-strata": 0.66,
    "post-processing:inhouse-af": 0.70,
    "post-processing:mane-refseq": 0.72,
    "post-processing:annotsv": 0.78,
    "sample-step:snv-overlay": 0.82,
    "sample-step:review-tsv": 0.86,
    "sample-step:gene-index": 0.92,
  };
  const within = sub[state.step] ?? 0;
  return _clampPct(82 + ((idx + within) / total) * 17);
}

function _toggleDragenLog() {
  const log = document.getElementById("dragen-job-log");
  const btn = document.getElementById("dragen-job-log-toggle");
  if (!log || !btn) return;
  log.hidden = !log.hidden;
  btn.textContent = log.hidden ? "▶ Log" : "▼ Log";
  btn.setAttribute("aria-expanded", log.hidden ? "false" : "true");
  if (!log.hidden) log.scrollTop = log.scrollHeight;
}

function _dragenStartPolling(jobId) {
  if (_DRAGEN_STATE.pollTimer) clearInterval(_DRAGEN_STATE.pollTimer);
  const tick = async () => {
    try {
      const s = await apiFetch(`/dragen/jobs/${encodeURIComponent(jobId)}`);
      if (!s) return;
      _dragenSetJob(s);
      if (s.state === "done" || s.state === "failed" || s.state === "cancelled") {
        clearInterval(_DRAGEN_STATE.pollTimer);
        _DRAGEN_STATE.pollTimer = null;
        // Refresh sample list so the new SID appears in the combobox.
        try { await loadIndex?.(); } catch (_e) {}
      }
    } catch (_e) {}
  };
  tick();
  _DRAGEN_STATE.pollTimer = setInterval(tick, 5000);
}

async function _dragenRecoverActiveJob() {
  // On page load — if a job is still running in the backend, surface
  // it in the topbar (so the user sees their previous run continuing).
  try {
    const jobs = await apiFetch("/dragen/jobs") || [];
    const active = jobs.find(j => j.running || j.state === "queued" || j.state === "running");
    if (
      active?.job_id
      && (active.job_id !== _DRAGEN_STATE.job?.job_id || !_DRAGEN_STATE.pollTimer)
    ) {
      _dragenStartPolling(active.job_id);
    }
  } catch (_e) {}
}

function _dragenStartRecoveryPolling() {
  if (_DRAGEN_STATE.recoverTimer) clearInterval(_DRAGEN_STATE.recoverTimer);
  _dragenRecoverActiveJob();
  _DRAGEN_STATE.recoverTimer = setInterval(_dragenRecoverActiveJob, 15000);
}

async function _dragenCancelCurrentJob() {
  const jobId = _DRAGEN_STATE.job?.job_id;
  if (!jobId) return;
  if (!confirm("確定要終止這次三級分析？正在執行的 Nextflow / post-processing 會收到終止訊號。")) return;
  const btn = document.getElementById("dragen-job-cancel-btn");
  if (btn) { btn.disabled = true; btn.textContent = "終止中…"; }
  try {
    const state = await apiFetch(`/dragen/jobs/${encodeURIComponent(jobId)}/cancel`, {
      method: "POST",
    });
    if (state) _dragenSetJob(state);
  } catch (e) {
    alert("終止失敗：" + (e?.message || e));
  } finally {
    if (btn) { btn.textContent = "終止"; }
  }
}

async function _dragenStart() {
  const extra = document.getElementById("dragen-extra-vep");
  const pgx = document.getElementById("dragen-pgx");
  let samples = _DRAGEN_STATE.batch.slice();
  if (!samples.length) {
    const sample = _dragenCurrentSampleFromForm();
    if (!sample) { alert("請先從 In-house 或 DRAGEN 清單挑選 VCF，或先加入批次"); return; }
    samples = [sample];
  }
  const mode = samples[0]?.mode || "";
  if (!mode) { alert("請先選擇 sample"); return; }
  const body = {
    mode,
    with_extra_vep: !!extra?.checked,
    with_pgx: pgx ? !!pgx.checked : true,
    samples,
  };
  const btn = document.getElementById("dragen-start-btn");
  if (btn) { btn.disabled = true; btn.textContent = "啟動中…"; }
  try {
    const r = await fetch(`${API_BASE}/dragen/jobs`, {
      method: "POST",
      credentials: "same-origin",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
    if (!r.ok) throw new Error(await r.text());
    const { job_id } = await r.json();
    _DRAGEN_STATE.batch = samples.slice();
    _dragenSetJob({
      job_id,
      state: "queued",
      step: "queued",
      running: true,
      mode,
      sample_id: samples[0]?.sample_id || samples[0]?.source_sample_id || "",
      sample_count: samples.length,
      log_tail: "",
    });
    _dragenRenderBatch();
    _dragenStartPolling(job_id);
  } catch (e) {
    alert("啟動失敗：" + (e?.message || e));
  } finally {
    if (btn) { btn.disabled = false; btn.textContent = "開始分析"; }
  }
}

function _dragenSuggestSid(vcfSid, mode) {
  if (!vcfSid) return "";
  const suffix = mode === "dragen" ? "-dragen" : mode === "inhouse" ? "-nckuh" : "";
  if (!suffix) return vcfSid;
  return vcfSid.toLowerCase().endsWith(suffix) ? vcfSid : `${vcfSid}${suffix}`;
}

function _dragenWireCombobox(mode) {
  const input = document.getElementById(`dragen-vcf-${mode}`);
  const list = document.getElementById(`dragen-vcf-${mode}-dropdown`);
  if (!input || !list) return;
  let activeIdx = -1;
  input.addEventListener("focus", () => {
    activeIdx = -1;
    _dragenRenderDropdown(mode, { showAll: true });
  });
  input.addEventListener("input", () => {
    activeIdx = -1;
    _DRAGEN_STATE.selected[mode] = "";
    input.title = "";
    _DRAGEN_STATE.mode = _dragenActiveMode();
    if (mode === "inhouse") _DRAGEN_STATE.seqType.inhouse = "";
    _dragenRenderSeqType();
    _dragenRenderSiblings();
    _dragenRenderDropdown(mode);
  });
  input.addEventListener("blur", () => {
    setTimeout(() => _dragenHideDropdown(mode), 120);
  });
  input.addEventListener("keydown", ev => {
    const opts = Array.from(list.querySelectorAll(".dragen-vcf-option[data-path]"));
    if (ev.key === "ArrowDown") {
      ev.preventDefault();
      activeIdx = Math.min(opts.length - 1, activeIdx + 1);
      opts.forEach((el, i) => el.classList.toggle("active", i === activeIdx));
    } else if (ev.key === "ArrowUp") {
      ev.preventDefault();
      activeIdx = Math.max(0, activeIdx - 1);
      opts.forEach((el, i) => el.classList.toggle("active", i === activeIdx));
    } else if (ev.key === "Enter" && activeIdx >= 0 && opts[activeIdx]) {
      ev.preventDefault();
      const row = _dragenCurrentList(mode).find(v => v.path === opts[activeIdx].dataset.path);
      if (row) _dragenPickVcf(mode, row);
    } else if (ev.key === "Escape") {
      _dragenHideDropdown(mode);
    }
  });
}

function setupDragenButton() {
  const btn = document.getElementById("btn-dragen-launch");
  if (!btn) return;
  btn.addEventListener("click", async () => {
    showModal("dragen-modal");
    await loadDragenVcfList();
  });
  _dragenWireCombobox("inhouse");
  _dragenWireCombobox("dragen");
  document.querySelectorAll(".dragen-seq-option").forEach(btn => {
    btn.addEventListener("click", () => {
      const seqType = String(btn.dataset.seqType || "").toUpperCase();
      if (seqType !== "WES" && seqType !== "WGS") return;
      _DRAGEN_STATE.seqType.inhouse = seqType;
      _dragenRenderSeqType();
    });
  });
  document.getElementById("dragen-refresh-btn")?.addEventListener("click", async (ev) => {
    const b = ev.currentTarget;
    if (b) { b.disabled = true; b.textContent = "↻ 更新中…"; }
    try { await loadDragenVcfList({ force: true }); }
    finally { if (b) { b.disabled = false; b.textContent = "↻ 更新索引"; } }
  });
  document.getElementById("dragen-start-btn")?.addEventListener("click", _dragenStart);
  document.getElementById("dragen-add-batch-btn")?.addEventListener("click", _dragenAddCurrentToBatch);
  document.getElementById("dragen-job-log-toggle")?.addEventListener("click", _toggleDragenLog);
  document.getElementById("dragen-job-cancel-btn")?.addEventListener("click", _dragenCancelCurrentJob);
  document.getElementById("topbar-job-status")?.addEventListener("click", () => {
    showModal("dragen-modal");
  });
  _dragenStartRecoveryPolling();
}

let _pipelineListRows = [];

function _drawPipelineListRows() {
  const host = document.getElementById("pipeline-list-table");
  const status = document.getElementById("pipeline-list-status");
  if (!host) return;
  const query = String(document.getElementById("pipeline-list-search")?.value || "").trim().toLowerCase();
  const rows = query
    ? _pipelineListRows.filter(row => [
        row.sample_id,
        row.source_sample_id,
        row.pipeline_sample_id,
      ].some(v => String(v || "").toLowerCase().includes(query)))
    : _pipelineListRows;
  if (!_pipelineListRows.length) {
    host.innerHTML = `<div class="muted" style="padding:10px">（尚無三級分析資料夾）</div>`;
    return;
  }
  if (!rows.length) {
    host.innerHTML = `<div class="muted" style="padding:10px">（沒有符合搜尋的 sample）</div>`;
    if (status) status.textContent = "";
    return;
  }
  const head = `
    <tr>
      <th>Sample</th>
      <th>狀態</th>
      <th>最近更新</th>
      <th>Log</th>
      <th></th>
    </tr>`;
  const body = rows.map(row => {
    const ready = row.has_acmg ? "完成" : (row.has_output ? "未完成" : "無輸出");
    const state = row.job_state ? ` · ${row.job_state}` : "";
    const source = row.source_sample_id && row.source_sample_id !== row.sample_id
      ? `<div class="muted">source: ${escapeHtml(row.source_sample_id)}</div>`
      : "";
    return `
      <tr>
        <td>${escapeHtml(row.sample_id || "")}${source}</td>
        <td>${escapeHtml(ready + state)}</td>
        <td>${escapeHtml(_dragenFmtMtime(row.mtime) || "—")}</td>
        <td><button type="button" class="btn btn-ghost pipeline-log-view"
          data-sample-id="${escapeAttr(row.sample_id || "")}">查看 Log</button></td>
        <td><button type="button" class="btn btn-danger pipeline-output-delete"
          data-sample-id="${escapeAttr(row.sample_id || "")}"
          data-pipeline-sample-id="${escapeAttr(row.pipeline_sample_id || row.source_sample_id || row.sample_id || "")}">刪除</button></td>
      </tr>`;
  }).join("");
  host.innerHTML = `<table>${head}${body}</table>`;
  if (status) status.textContent = query ? `符合 ${rows.length} / ${_pipelineListRows.length} 筆` : "";
}

async function _renderPipelineList() {
  const host = document.getElementById("pipeline-list-table");
  if (!host) return;
  host.innerHTML = `<div class="muted" style="padding:10px">載入中…</div>`;
  try {
    _pipelineListRows = await apiFetch("/dragen/outputs") || [];
    _drawPipelineListRows();
  } catch (e) {
    host.innerHTML = `<div class="muted" style="padding:10px">載入失敗：${escapeHtml(String(e))}</div>`;
  }
}

function _fmtBytes(bytes) {
  const n = Number(bytes) || 0;
  if (n < 1024) return `${n} B`;
  const units = ["KB", "MB", "GB", "TB"];
  let v = n / 1024;
  let idx = 0;
  while (v >= 1024 && idx < units.length - 1) {
    v /= 1024;
    idx += 1;
  }
  return `${v.toFixed(v >= 10 ? 1 : 2)} ${units[idx]}`;
}

async function _cleanupNextflowWork() {
  const btn = document.getElementById("pipeline-clean-nf-work");
  const status = document.getElementById("pipeline-list-status");
  if (!confirm("確定清理 Nextflow 暫存？\n\n將刪除 /home/n102968/NGS_UI/nf_work 底下的所有內容。\n正在執行三級分析時後端會拒絕清理。")) return;
  if (btn) { btn.disabled = true; btn.textContent = "清理中…"; }
  if (status) status.textContent = "清理 Nextflow 暫存中…";
  try {
    const result = await apiFetch("/dragen/nf-work/cleanup", { method: "POST" });
    const count = Number(result?.deleted_count || 0);
    const size = _fmtBytes(result?.freed_bytes || 0);
    const path = result?.path || "/home/n102968/NGS_UI/nf_work";
    if (status) status.textContent = `已清理 ${count} 個項目，約釋放 ${size}：${path}`;
  } catch (e) {
    if (status) status.textContent = `清理失敗：${e.message || e}`;
  } finally {
    if (btn) { btn.disabled = false; btn.textContent = "清理 Nextflow 暫存"; }
  }
}

function setupPipelineList() {
  const btn = document.getElementById("btn-pipeline-list");
  const host = document.getElementById("pipeline-list-table");
  const status = document.getElementById("pipeline-list-status");
  const logPanel = document.getElementById("pipeline-list-log-panel");
  if (!btn || !host) return;
  btn.addEventListener("click", async () => {
    showModal("pipeline-list-modal");
    if (status) status.textContent = "";
    if (logPanel) logPanel.hidden = true;
    const search = document.getElementById("pipeline-list-search");
    if (search) search.value = "";
    await _renderPipelineList();
  });
  document.getElementById("pipeline-list-search")?.addEventListener("input", _drawPipelineListRows);
  document.getElementById("pipeline-clean-nf-work")?.addEventListener("click", _cleanupNextflowWork);
  host.addEventListener("click", async ev => {
    const logBtn = ev.target.closest?.(".pipeline-log-view");
    const delBtn = ev.target.closest?.(".pipeline-output-delete");
    const sid = (logBtn || delBtn)?.dataset.sampleId || "";
    if (!sid) return;
    if (logBtn) {
      const title = document.getElementById("pipeline-list-log-title");
      const log = document.getElementById("pipeline-list-log");
      if (title) title.textContent = `${sid} Log`;
      if (log) log.textContent = "載入中…";
      if (logPanel) logPanel.hidden = false;
      try {
        const data = await apiFetch(`/dragen/outputs/${encodeURIComponent(sid)}/log`);
        if (log) log.textContent = data?.log || "（沒有可用的 NGS-UI 三級分析 Log）";
      } catch (e) {
        if (log) log.textContent = `讀取失敗：${e.message || e}`;
      }
      return;
    }
    const pipelineSid = delBtn.dataset.pipelineSampleId || sid;
    if (!confirm(`確定刪除三級分析原始檔案、NGS-UI 個案資料與 job log？\n\n/home/datalake_Intermediate/pipeline/tertiary_output/${sid}/\n\n若為舊個案，也會清除：\n/home/pipeline/tertiary_output/${pipelineSid}/\nNGS_UI/tertiary_output/${sid}/\n\n此操作無法復原。`)) return;
    delBtn.disabled = true;
    if (status) status.textContent = `刪除 ${sid} 中…`;
    try {
      const resp = await fetch(`${API_BASE}/dragen/outputs/${encodeURIComponent(sid)}`, {
        method: "DELETE",
        credentials: "same-origin",
      });
      if (resp.status === 401) { showLoginModal(); throw new Error("尚未登入"); }
      const body = await resp.json().catch(() => ({}));
      if (!resp.ok) throw new Error(body.detail || `${resp.status} ${resp.statusText}`);
      if (status) status.textContent = `已刪除 ${sid}。`;
      if (logPanel) logPanel.hidden = true;
      await _renderPipelineList();
    } catch (e) {
      if (status) status.textContent = `刪除失敗：${e.message || e}`;
      delBtn.disabled = false;
    }
  });
  document.getElementById("pipeline-list-log-close")?.addEventListener("click", () => {
    if (logPanel) logPanel.hidden = true;
  });
}

// Wire it up at boot — sits alongside setupCombobox / setupGeneSearch.
document.addEventListener("DOMContentLoaded", () => {
  try { setupSecondaryButton(); } catch (_e) {}
  try { setupDragenButton(); } catch (_e) {}
  try { setupPipelineList(); } catch (_e) {}
});
