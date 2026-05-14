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
};

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
  state.index = list.map(r => ({
    LIS_ID:    r.lis_id || r.sample_id || r.LIS_ID || "",
    Name:      r.name || r.Name || "",
    MRN:       r.mrn || r.MRN || "",
    Test:      r.test_type || r.Test || "",
    Category:  r.category || r.Category || "",
    Tag:       (r.tags || []).join(",") || r.Tag || "",
    sample_id: r.sample_id || r.lis_id || r.LIS_ID || "",
  }));
  const opts = await apiFetch("/options").catch(() => null);
  state.options = opts && typeof opts === "object"
    ? opts
    : { category_options: [], tag_suggestions: [] };
}

function matchSamples(q) {
  if (!state.index) return [];
  const s = (q || "").trim().toLowerCase();
  if (!s) return state.index.slice(0, 50);
  return state.index.filter(r => {
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
  // Exact match on LIS_ID first, then MRN, then single-hit Name
  const hitLis = state.index.find(r => (r.LIS_ID || "").trim() === s);
  if (hitLis) return hitLis;
  const hitMrn = state.index.find(r => (r.MRN || "").trim() === s);
  if (hitMrn) return hitMrn;
  const hitNames = state.index.filter(r => (r.Name || "").trim() === s);
  if (hitNames.length === 1) return hitNames[0];
  if (hitNames.length > 1) throw new Error(`找到 ${hitNames.length} 筆同名樣本，請改用 LIS_ID 或病歷號`);
  // Fallback: single partial match across fields
  const partial = matchSamples(s);
  if (partial.length === 1) return partial[0];
  return null;
}

async function loadSample(LIS_ID) {
  // The combobox carries the row's `sample_id` (which equals LIS_ID for
  // legacy samples). Look it up so we use the directory name the backend
  // wants on the URL path.
  const row = (state.index || []).find(r => r.LIS_ID === LIS_ID);
  const sid = row?.sample_id || LIS_ID;
  const data = await apiFetch(`/samples/${encodeURIComponent(sid)}`);
  if (!data) throw new Error(`找不到 sample ${sid}`);
  const reports = await apiFetch(`/samples/${encodeURIComponent(sid)}/report`) || {
    status: {}, edits: {}, panels: {}, clinical_description: "",
    genetic_counseling: "", comment: "",
    category: null, yield: 0, updated_at: null,
  };
  if (reports.clinical_description == null) reports.clinical_description = "";
  if (reports.genetic_counseling   == null) reports.genetic_counseling   = data.genetic_counseling || "";
  if (reports.comment == null)               reports.comment = "";
  if (!Array.isArray(reports.tags))          reports.tags = [];
  if (!Array.isArray(reports.manual_variants)) reports.manual_variants = [];
  state.data       = data;
  state.reports    = reports;
  state.currentLIS = LIS_ID;
  state.dirty      = false;
  _saveError       = "";
  _lastSavedAt     = null;
  clearTimeout(_autoSaveTimer);
  // Reset manual-toggle tracking between samples so defaultOpen applies fresh.
  toggledBlocks.clear();

  // Staged loading: the core payload above carries empty CNV/SV + Mito
  // side-channels (aux_pending). Pull them in the background so the
  // SNV/Indel view + report sections appear immediately; each card
  // re-renders itself when its data lands. A monotonic token drops a
  // stale response that arrives after the user switched samples.
  if (data.aux_pending) {
    const token = (state._auxLoadToken = (state._auxLoadToken || 0) + 1);
    state.cnvSvPending = true;
    state.mitoPending  = true;
    apiFetch(`/samples/${encodeURIComponent(sid)}/cnv-sv`)
      .then(aux => {
        if (token !== state._auxLoadToken || !state.data) return;
        if (aux) Object.assign(state.data, aux);
        state.cnvSvPending = false;
        try { renderCnvSvTabBar(); } catch (_e) {}
      })
      .catch(() => {
        if (token !== state._auxLoadToken) return;
        state.cnvSvPending = false;
        try { renderCnvSvTabBar(); } catch (_e) {}
      });
    apiFetch(`/samples/${encodeURIComponent(sid)}/mito`)
      .then(aux => {
        if (token !== state._auxLoadToken || !state.data) return;
        if (aux) Object.assign(state.data, aux);
        state.mitoPending = false;
        try { renderMitoTabBar(); } catch (_e) {}
      })
      .catch(() => {
        if (token !== state._auxLoadToken) return;
        state.mitoPending = false;
        try { renderMitoTabBar(); } catch (_e) {}
      });
  } else {
    state.cnvSvPending = false;
    state.mitoPending  = false;
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

function formatClinvar(sig, conf, stars) {
  // Treat pipeline placeholders ('.', 'NA', '') as "no ClinVar data"
  // — otherwise the cell renders as `.(0★)` instead of `—`.
  const sigStr = (sig == null ? "" : String(sig)).trim();
  if (!sigStr || sigStr === "." || sigStr.toUpperCase() === "NA" || sigStr.toUpperCase() === "N/A") {
    return "—";
  }
  const starTxt = (stars != null && stars !== "") ? `(${stars}★)` : "";
  if (sigStr.startsWith("Conflicting") && conf) {
    const parts = String(conf).split(/[,|]/).map(p => p.trim().replace(/^_/, ""));
    const out = parts.map(p => {
      const m = p.match(/^(.+?)\((\d+)\)$/);
      if (!m) return p;
      return (CLINVAR_ABBREV[m[1]] || m[1]) + "(" + m[2] + ")";
    });
    return out.join("|") + starTxt;
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
    case "pathogenic/likely pathogenic":      case "p/lp":  return "sig-p";
    case "likely pathogenic":                 case "lp":    return "sig-lp";
    case "uncertain significance":            case "vus":   return "sig-vus";
    case "likely benign":                     case "lb":    return "sig-lb";
    case "benign":                            case "b":     return "sig-b";
    case "benign/likely benign":              case "b/lb":  return "sig-b";
    default:                                                return null;
  }
}

// Numeric classifiers for in-silico tools. Cutoffs are the ClinGen-recommended
// calibrations:
//   AlphaMissense / MetaRNN — Pejaver V et al, AJHG 2022 (PMID 36413997),
//   the same calibration the R script uses for its ps grading.
//   SpliceAI — Walker LC et al, AJHG 2023 (PMID 37352859), ClinGen SVI
//   Splicing Subgroup recommendations.
// Each entry holds the four boundary scores; classifyByThresholds() walks
// them top-down and returns the first matching sig-* class.
const TOOL_CUTOFFS = {
  alphamissense: { p: 0.990, lp: 0.792, vus: 0.170, lb: 0.071 },
  metarnn:       { p: 0.939, lp: 0.748, vus: 0.267, lb: 0.108 },
  // SpliceAI is a splice-impact probability, not a pathogenicity score —
  // a low value just means "no predicted splice change", not "benign", so
  // anything below the VUS threshold stays uncoloured (no lb / no implicit B).
  spliceai:      { p: 0.800, lp: 0.500, vus: 0.200 },
};
function classifyByThresholds(score, cutoffs) {
  if (score == null || score === "") return null;
  const x = Number(score);
  if (!Number.isFinite(x)) return null;
  if (cutoffs.p   != null && x >= cutoffs.p)   return "sig-p";
  if (cutoffs.lp  != null && x >= cutoffs.lp)  return "sig-lp";
  if (cutoffs.vus != null && x >= cutoffs.vus) return "sig-vus";
  if (cutoffs.lb  != null && x >= cutoffs.lb)  return "sig-lb";
  // Below all configured thresholds. Tools that have an `lb` cutoff
  // (real pathogenicity scores) get the dark-green B chip; tools that
  // omit `lb` (e.g. SpliceAI) leave low-end scores uncoloured.
  return cutoffs.lb == null ? null : "sig-b";
}

// LoGoFunc emits strings like "GOF (0.123)*", "LOF (0.456)", or "Neutral (...)".
// A trailing star means probability > class-specific cutoff (deeper red);
// no star but class is GOF / LOF gives a lighter red. Neutral / NA is uncoloured.
function classifyLoGoFunc(text) {
  if (text == null || text === "" || text === "—") return null;
  const s = String(text).trim();
  const m = s.match(/^(GOF|LOF|Neutral)/i);
  if (!m) return null;
  if (m[1].toUpperCase() === "NEUTRAL") return null;
  return s.endsWith("*") ? "sig-p" : "sig-lp";
}

// MaxEntScan diff cutoffs from the Pejaver-style PS3 calibration the
// R script uses internally: |diff| ≥ 7.65 strong / 5.96 moderate /
// 4.24 supporting; below that no colour.
function classifyMaxEntScan(score) {
  if (score == null || score === "") return null;
  const x = Math.abs(Number(score));
  if (!Number.isFinite(x)) return null;
  if (x >= 7.65) return "sig-p";
  if (x >= 5.96) return "sig-lp";
  if (x >= 4.24) return "sig-vus";
  return null;
}

// PDIVAS cutoffs: ≥ 0.5 high-confidence pathogenic intronic, ≥ 0.082
// (paper default) supporting; below that no colour.
// Source: Kurosawa R et al, BMC Genomics 2023.
function classifyPDIVAS(score) {
  if (score == null || score === "") return null;
  const x = Number(score);
  if (!Number.isFinite(x)) return null;
  if (x >= 0.5)   return "sig-p";
  if (x >= 0.082) return "sig-lp";
  return null;
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
function fmtAdVaf(ad, vaf) {
  const adPart  = (ad == null || ad === "") ? "—" : String(ad);
  const vafPart = (vaf == null || vaf === "" || !Number.isFinite(Number(vaf)))
    ? "—"
    : fmtNum(vaf);
  return `${adPart} (${vafPart})`;
}

function variantUrls(v) {
  const tag = `${v.CHROM}-${v.POS}-${v.REF}-${v.ALT}`;
  // Route Varsome / Franklin / GeneBe to the matching genome build so
  // hg19 samples don't land on an hg38 coordinate page. Build info
  // comes from the R webdata writer (state.data.genome_build); falls
  // back to hg38 for older samples that predate the field.
  const build = state.data?.genome_build === "hg19" ? "hg19" : "hg38";
  return {
    varsome:  `https://varsome.com/variant/${build}/${tag}`,
    franklin: `https://franklin.genoox.com/clinical-db/variant/snp/${tag}-${build}`,
    genebe:   `https://genebe.net/variant/${build}/${tag}`,
    omim:     v.OMIM_link || (v.OMIM_id ? `https://www.omim.org/entry/${v.OMIM_id}` : null),
  };
}

// ---------- Render: sample header / phenotype ----------------------

function renderSampleMeta() {
  const m = state.data.meta || {};
  document.getElementById("m-lis").textContent       = m.LIS_ID || "—";
  document.getElementById("m-name").textContent      = m.Name || "—";
  document.getElementById("m-mrn").textContent       = m.MRN || "—";
  // Generated date: only the YYYY-MM-DD prefix; the time/tz suffix
  // adds noise without telling the reviewer anything they care about.
  const gen = state.data.generated_at || "";
  document.getElementById("m-generated").textContent = gen ? gen.slice(0, 10) : "—";

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

  const sel = document.getElementById("m-category");
  const opts = (state.options && state.options.category_options) || [];
  const current = state.reports.category ?? m.Category ?? "";
  const all = Array.from(new Set(["", ...opts, current].filter(x => x !== undefined && x !== null)));
  sel.innerHTML = all.map(o => {
    const label = o === "" ? "—" : o;
    return `<option value="${escapeAttr(o)}" ${o === current ? "selected" : ""}>${escapeHtml(label)}</option>`;
  }).join("");

  document.getElementById("sample-card").classList.remove("hidden");
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
  if (ev.target.id === "m-test")  _saveSampleMeta({ test_type:    ev.target.value });
  if (ev.target.id === "m-build") _saveSampleMeta({ genome_build: ev.target.value });
  if (ev.target.id === "m-sex")   _saveSampleMeta({ sex:          ev.target.value });
});

// Generic renderer for collapsible free-text cards (Clinical presentation,
// Comment). Both default to collapsed; user-toggled state is remembered
// across re-renders via toggledBlocks.
function renderCollapsibleCard(cardId, headerId, bodyId, taId, value) {
  const card   = document.getElementById(cardId);
  const header = document.getElementById(headerId);
  const body   = document.getElementById(bodyId);
  const ta     = document.getElementById(taId);

  ta.value = value || "";
  const open = toggledBlocks.has(cardId)
    ? card.dataset.wasOpen === "1"
    : false;
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
                        "comment-text", state.reports.comment);
  renderTagPicker();
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
  // If there's a queued/running job for this sample, pick up polling.
  _resumeJobPollingIfAny();
}

async function _resumeJobPollingIfAny() {
  if (!state.currentLIS) return;
  clearInterval(_jobPollTimer);
  _setJobStatus("");
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
    li.innerHTML = `
      <span class="chip-label">${escapeHtml(row.name)}</span>
      <select class="chip-weight" data-panel-idx="${idx}" title="Weight">
        ${[1,2,3,4,5].map(n => `<option value="${n}" ${n===row.weight?"selected":""}>w=${n}</option>`).join("")}
      </select>
      <button class="chip-remove" data-panel-idx="${idx}" type="button" title="移除">×</button>`;
    ul.appendChild(li);
  });
}

// Cached panel list from /api/panels — fetched once per session.
let _panelOptions = null;
async function loadPanelOptions() {
  if (_panelOptions) return _panelOptions;
  _panelOptions = await apiFetch("/panels") || [];
  return _panelOptions;
}

// HPO + panel typeaheads use delegated listeners on `document` so they
// keep working even if the phenotype-card is re-rendered. Per-element
// addEventListener was racing the legacy global handlers; delegated
// dispatch sidesteps it entirely.
let _hpoSearchAbort = null;
let _hpoSearchTimer = null;

function _hpoOpen()  { document.getElementById("hpo-search-dropdown")?.classList.remove("hidden"); }
function _hpoClose() { document.getElementById("hpo-search-dropdown")?.classList.add("hidden"); }
function _panelOpen()  { document.getElementById("panel-search-dropdown")?.classList.remove("hidden"); }
function _panelClose() { document.getElementById("panel-search-dropdown")?.classList.add("hidden"); }

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
    _hpoOpen();
  } catch (e) {
    if (e.name !== "AbortError") console.error("HPO search failed", e);
  }
}

async function _runPanelSearch(q) {
  const dropdown = document.getElementById("panel-search-dropdown");
  if (!dropdown) return;
  const opts = await loadPanelOptions();
  const ql = (q || "").trim().toLowerCase();
  const picked = new Set(phenoEdit.panels.map(p => p.name));
  const matches = opts
    .filter(p => !picked.has(p.name) && (!ql || p.name.toLowerCase().includes(ql)))
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

// Mousedown so the option fires before the input's blur kills the dropdown.
document.addEventListener("mousedown", ev => {
  const opt = ev.target.closest(".combobox-option");
  if (!opt) return;
  const dropdown = opt.parentElement;
  if (dropdown?.id === "hpo-search-dropdown") {
    ev.preventDefault();
    addHpo(opt.dataset.id, opt.dataset.name);
    const inp = document.getElementById("hpo-search");
    if (inp) inp.value = "";
    _hpoClose();
  } else if (dropdown?.id === "panel-search-dropdown") {
    ev.preventDefault();
    addPanel(opt.dataset.name);
    const inp = document.getElementById("panel-search");
    if (inp) inp.value = "";
    _panelClose();
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
// box. Typing a symbol + Enter opens a modal listing every variant
// of that gene as cards — the same renderers (renderVariantCard /
// renderCnvSvCard) the tier tables use, so the cards are fully
// interactive (status dropdown, disease list, comment, …). The
// modal's own input lets the reviewer pivot to another gene without
// closing it.

function _geneSearchSnv(geneUpper) {
  const matches = Object.entries(state.data?.variants || {})
    .filter(([, v]) => (v.gene_symbol || "").toUpperCase() === geneUpper);
  matches.sort((a, b) => {
    const sa = Number(a[1].total_score), sb = Number(b[1].total_score);
    return (Number.isFinite(sb) ? sb : -Infinity) - (Number.isFinite(sa) ? sa : -Infinity);
  });
  return matches;
}

function _geneSearchCnvSv(geneUpper) {
  const all = [
    ...Object.entries(state.data?.cnv_variants || {}),
    ...Object.entries(state.data?.sv_variants  || {}),
  ];
  const matches = all.filter(([, v]) =>
    (v.gene_list || []).some(g => String(g).toUpperCase() === geneUpper)
  );
  matches.sort((a, b) => {
    const ra = Number(a[1].ranking_score), rb = Number(b[1].ranking_score);
    return (Number.isFinite(rb) ? rb : -Infinity) - (Number.isFinite(ra) ? ra : -Infinity);
  });
  return matches;
}

function renderGeneSearchResults(kind, geneUpper) {
  const titleEl = document.getElementById("gene-search-title");
  const host    = document.getElementById("gene-search-results");
  if (!titleEl || !host) return;
  host.innerHTML = "";
  if (!geneUpper) {
    titleEl.textContent = "基因變異搜尋";
    host.innerHTML = `<div class="muted" style="padding:12px">輸入基因名稱以搜尋。</div>`;
    return;
  }
  if (kind === "all") {
    const snvMatches = _geneSearchSnv(geneUpper);
    const cnvMatches = _geneSearchCnvSv(geneUpper);
    titleEl.textContent = `${geneUpper} 的所有變異（SNV/Indel: ${snvMatches.length}，CNV/SV: ${cnvMatches.length}）`;
    if (!snvMatches.length && !cnvMatches.length) {
      host.innerHTML = `<div class="muted" style="padding:12px">找不到 ${escapeHtml(geneUpper)} 的變異。</div>`;
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
    const matches = _geneSearchSnv(geneUpper);
    titleEl.textContent = `${geneUpper} 的 ${label} 變異（${matches.length}）`;
    if (!matches.length) {
      host.innerHTML = `<div class="muted" style="padding:12px">找不到 ${escapeHtml(geneUpper)} 的變異。</div>`;
      return;
    }
    matches.forEach(([id, v], i) => {
      host.appendChild(renderVariantCard(v, id, "candidate", { index: i + 1, diseaseCheckbox: true }));
    });
  } else {
    const matches = _geneSearchCnvSv(geneUpper);
    titleEl.textContent = `${geneUpper} 的 ${label} 變異（${matches.length}）`;
    if (!matches.length) {
      host.innerHTML = `<div class="muted" style="padding:12px">找不到涵蓋 ${escapeHtml(geneUpper)} 的 CNV/SV。</div>`;
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
  if (inp) {
    inp.style.display = "";
    inp.value = gene || "";
    inp.dataset.kind = kind;
  }
  renderGeneSearchResults(kind, (gene || "").trim().toUpperCase());
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
    renderGeneSearchResults(inp.dataset.kind || "snv", inp.value.trim().toUpperCase());
  });
}

// 上傳個案清單: pick an xlsx → POST /api/patient_list → toast the
// {added, updated, total} result. The roster it builds is what the
// 載入新個案 modal reads to auto-fill MRN / 姓名 / Test type.
function setupPatientListUpload() {
  const btn  = document.getElementById("btn-upload-list");
  const file = document.getElementById("upload-list-file");
  if (!btn || !file) return;
  btn.addEventListener("click", () => { file.value = ""; file.click(); });
  file.addEventListener("change", async () => {
    const f = file.files && file.files[0];
    if (!f) return;
    const origLabel = btn.textContent;
    btn.disabled = true;
    btn.textContent = "上傳中…";
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
      alert(`個案清單上傳完成\n\n解析 ${body.parsed} 筆 · 新增 ${body.added} · 更新 ${body.updated}\nroster 目前共 ${body.total} 筆`);
      // Refresh the unregistered-sample cache so an open 載入新個案
      // modal picks up the new MRN/name mappings next time it opens.
      _unregisteredById = {};
    } catch (e) {
      alert("個案清單上傳失敗：" + (e.message || e));
    } finally {
      btn.disabled = false;
      btn.textContent = origLabel;
    }
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

// Download diagnostic DOCX. Streams from the backend so the user
// just sees a normal browser save dialog. Filename comes from the
// Content-Disposition header (RFC 5987 UTF-8 encoded).
async function exportDiagnosticDocx() {
  if (!state.currentLIS) return;
  const row = (state.index || []).find(r => r.LIS_ID === state.currentLIS);
  const sid = row?.sample_id || state.currentLIS;
  try {
    const resp = await fetch(`${API_BASE}/samples/${encodeURIComponent(sid)}/report.docx`, {
      credentials: "same-origin",
    });
    if (resp.status === 401) { showLoginModal(); return; }
    if (!resp.ok) throw new Error(`匯出失敗 (${resp.status})`);
    const blob = await resp.blob();
    downloadBlob(blob, `${sid}_diagnosis.docx`);
  } catch (e) {
    alert("匯出失敗：" + e.message);
  }
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

// 「開始分析」: instant in-house pheno_score + queued Exomiser/LIRICAL.
// The phenotype POST returns immediately so the cards' IN_PANEL badges
// flip right away; the Exomiser/LIRICAL job runs in the background and
// the polling loop refreshes the sample once it lands.
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
    // Refresh so cards see the freshly written IN_PANEL flag right
    // away, before the slower Exomiser job finishes.
    await loadSample(state.currentLIS);
    renderAll();
  } catch (e) {
    hint.textContent = "失敗：" + e.message;
    _setJobStatus("失敗", false);
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
          // Reload sample so cards pick up the new score columns.
          await loadSample(state.currentLIS);
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
  // kind: "candidate" → causative/other/candidate/skip/reject mapped to
  //                     1 / 2 / C / 0 / X. "C" routes the variant into
  //                     the Candidate variants report section.
  //       "panel"     → ACMG SF / Proactive / Carrier → V/0/X.
  // "0" = reviewed but kept on the page (surfaces nowhere in Report).
  if (kind === "panel") return ["", "V", "0", "X"];
  return ["", "1", "2", "C", "0", "X"];
}

function getStatus(id) {
  return (state.reports.status && state.reports.status[id]) || "";
}

function setStatus(id, val) {
  state.reports.status = state.reports.status || {};
  if (val) state.reports.status[id] = val;
  else     delete state.reports.status[id];
  state.dirty = true;
  renderAll();
}

// Panel-specific status (per panel category, so V in proactive doesn't surface in carrier)
function getPanelStatus(id, panel) {
  const m = state.reports.panels && state.reports.panels[id];
  return (m && m[panel]) || "";
}

function setPanelStatus(id, panel, val) {
  state.reports.panels = state.reports.panels || {};
  state.reports.panels[id] = state.reports.panels[id] || {};
  if (val) state.reports.panels[id][panel] = val;
  else     delete state.reports.panels[id][panel];
  if (Object.keys(state.reports.panels[id]).length === 0) {
    delete state.reports.panels[id];
  }
  state.dirty = true;
  renderAll();
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
        <select class="status-select" data-id="${escapeAttr(id)}" ${panelAttr}>
          ${options.map(s => `<option value="${s}" ${s===curStatus?"selected":""}>${s||"—"}</option>`).join("")}
        </select>
      </div>`;
    return card;
  }

  const urls = variantUrls(v);
  const card = document.createElement("div");
  card.className = "variant-card";
  card.dataset.inPanel = v.in_panel ? "true" : "false";

  const links = [
    `<a href="${urls.varsome}"  target="_blank" rel="noopener">Varsome</a>`,
    `<a href="${urls.franklin}" target="_blank" rel="noopener">Franklin</a>`,
    `<a href="${urls.genebe}"   target="_blank" rel="noopener">GeneBe</a>`,
    urls.omim ? `<a href="${urls.omim}"     target="_blank" rel="noopener">OMIM</a>` : "",
  ].join("");

  const editAcmgClass = getEdit(id, "ACMG_classification") ?? v.ACMG_classification ?? "";
  const editAcmgCrit  = getEdit(id, "ACMG_criteria")       ?? v.ACMG_criteria       ?? "";
  const editAcmgScore = getEdit(id, "ACMG_score")          ?? (v.ACMG_score ?? "");
  const editComment   = getEdit(id, "comment")             ?? "";

  const clinvarDate = formatClinvarDate(state.data?.clinvar_date);
  const clinvarLabel = clinvarDate ? `ClinVar (${escapeHtml(clinvarDate)})` : "ClinVar";

  // Extras shown only when the user clicks the "More" button. Each row is
  // pushed only when the underlying field is present in the webdata, so a
  // sample with none of these fields has the More button hidden too.
  const extras = [];
  if (v.in_silico_prediction != null && v.in_silico_prediction !== "") {
    extras.push({ key: "In silico prediction",
                  html: fmtInSilico(v.in_silico_prediction) });
  }
  if (v.LoGoFunc != null && v.LoGoFunc !== "" && v.LoGoFunc !== "NA") {
    extras.push({ key: "LoGoFunc", text: String(v.LoGoFunc),
                  cls: classifyLoGoFunc(v.LoGoFunc) });
  }
  if (v.MaxEntScan_diff != null && v.MaxEntScan_diff !== "") {
    extras.push({ key: "MaxEntScan", text: fmtNum(v.MaxEntScan_diff),
                  cls: classifyMaxEntScan(v.MaxEntScan_diff) });
  }
  if (v.PDIVAS_score != null && v.PDIVAS_score !== "") {
    extras.push({ key: "PDIVAS", text: fmtNum(v.PDIVAS_score),
                  cls: classifyPDIVAS(v.PDIVAS_score) });
  }
  // Tertiary-spec extras: P-KNN / REVEL / BayesDel / ESM2 / Evo2 / CADD.
  // P-KNN is the primary missense PP3/BP4 metric per the 2026 spec; the
  // others are secondary references shown only when populated.
  if (v.PKNN_LLR != null && v.PKNN_LLR !== "") {
    extras.push({ key: "P-KNN LLR", text: fmtNum(v.PKNN_LLR) });
  }
  if (v.REVEL != null && v.REVEL !== "") {
    extras.push({ key: "REVEL", text: fmtNum(v.REVEL) });
  }
  if (v.BayesDel != null && v.BayesDel !== "") {
    extras.push({ key: "BayesDel", text: fmtNum(v.BayesDel) });
  }
  if (v.CADD_phred != null && v.CADD_phred !== "") {
    extras.push({ key: "CADD", text: fmtNum(v.CADD_phred) });
  }
  if (v.ESM2_score != null && v.ESM2_score !== "") {
    extras.push({ key: "ESM-2", text: fmtNum(v.ESM2_score) });
  }
  if (v.Evo2_score != null && v.Evo2_score !== "") {
    extras.push({ key: "Evo2", text: fmtNum(v.Evo2_score) });
  }
  if (v.loftee_hc || v.loftee_filter || v.loftee_flags) {
    const parts = [v.loftee_hc, v.loftee_filter, v.loftee_flags]
      .filter(Boolean).join(" / ");
    extras.push({ key: "LOFTEE", text: parts });
  }
  const extrasHtml = extras.map(x => {
    const valHtml = x.html != null ? x.html : escapeHtml(x.text);
    return `<span class="k">${escapeHtml(x.key)}</span>`
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
      <select class="status-select" data-id="${escapeAttr(id)}" ${panelAttr}>
        ${options.map(s => `<option value="${s}" ${s===curStatus?"selected":""}>${s||"—"}</option>`).join("")}
      </select>
      <span class="hgvs">${v.clinvar_upgrade ? `<span class="clinvar-upgrade-arrow" title="ClinVar 升級">${escapeHtml(v.clinvar_upgrade)}</span> ` : ""}${escapeHtml(v.HGVS || id)}<button class="btn-copy" data-copy="${escapeAttr(v.HGVS || id)}" title="複製 HGVS">${COPY_ICON_SVG}</button> <span class="variant-tag">([${escapeHtml(state.data?.genome_build || "hg38")}] ${escapeHtml(id)}<button class="btn-copy" data-copy="${escapeAttr(id)}" title="複製 chr-pos-ref-alt">${COPY_ICON_SVG}</button>)</span></span>
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
        <span class="k">Score</span><span class="v">${escapeHtml(scoreLine)}</span>
        <span class="k">Exomiser</span><span class="v">${escapeHtml(fmtScoreRank(v.total_score_exomiser_variant, v.rank_exomiser_variant))}</span>
        <span class="k">LIRICAL</span><span class="v">${escapeHtml(fmtScoreRank(v.lirical_variant_score, v.rank_lirical_variant))}</span>
      </div>
      <div>
        <span class="k">Zygosity</span><span class="v">${fmtTxt(v.zygosity)}</span>
        <span class="k">Read depth (VAF)</span><span class="v">${escapeHtml(fmtAdVaf(v.AD, v.alt_af))}</span>
        <span class="k">Consequence</span><span class="v">${fmtTxt(v.Consequence)}</span>
        <div class="more-extras hidden">
          <span class="k">Exon / Intron</span><span class="v">${fmtExonIntron(v)}</span>
          <span class="k">Phase</span><span class="v">${fmtPhase(v)}</span>
        </div>
      </div>
      <div>
        <span class="k">${clinvarLabel}</span><span class="v ${classifySignificance(v.CLNSIG) || ""}">${escapeHtml(formatClinvar(v.CLNSIG, v.CLNSIGCONF, v.clinvar_stars))}${v.clinvar_upgrade && v.CLNSIG_old ? ` <span class="clinvar-old" title="原 ClinVar 分類">(was: ${escapeHtml(formatClinvar(v.CLNSIG_old, v.CLNSIGCONF_old, v.clinvar_stars_old))})</span>` : ""}</span>
        <span class="k">ACMG</span>
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
        <span class="k">AlphaMissense</span><span class="v ${classifyByThresholds(v.AlphaMissense_score, TOOL_CUTOFFS.alphamissense) || ''}">${fmtNum(v.AlphaMissense_score)}</span>
        <span class="k">MetaRNN</span><span class="v ${classifyByThresholds(v.MetaRNN_score, TOOL_CUTOFFS.metarnn) || ''}">${fmtNum(v.MetaRNN_score)}</span>
        <span class="k">SpliceAI</span><span class="v ${classifyByThresholds(v.SpliceAI_score, TOOL_CUTOFFS.spliceai) || ''}">${fmtNum(v.SpliceAI_score)}</span>
        ${extras.length ? `<div class="more-extras hidden">${extrasHtml}</div>` : ""}
      </div>
      <div>
        <span class="k">AF</span><span class="v">${fmtNum(v.AF, 5)}</span>
        <span class="k">AF_eas</span><span class="v">${fmtNum(v.AF_eas, 5)}</span>
        ${Number.isFinite(Number(v.TaiwanBioBank))
          ? `<span class="k">TWB</span><span class="v">${fmtNum(v.TaiwanBioBank, 5)}</span>`
          : ""}
      </div>
    </div>
    <button class="btn-more" type="button">▾ More</button>
    <div class="more-extras hidden">${renderManeAll(v)}</div>
    ${renderDiseaseList(v, id, !!opts.diseaseCheckbox)}
  `;

  return card;
}

// ---- helpers used by renderVariantCard (Phase 4) -----------------

// Top-of-card chip row: TRANSCRIPT_TYPE / CALLERS / panel/ROH/blacklist
// hits / LOFTEE HC. Empty-string entries are filtered out so the row
// hides itself when nothing is worth showing.
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
function fmtExonIntron(v) {
  const e = String(v.exon || "").trim();
  if (e) return `exon ${e}`;
  const i = String(v.intron || "").trim();
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

// Inline <details> block listing every MANE transcript carried by
// MANE_ALL. The pipeline emits this as an array of
// {transcript, transcript_type, hgvs_c, hgvs_p, consequence}; older
// rows leave the list empty and we hide the block in that case.
function renderManeAll(v) {
  const rows = Array.isArray(v.MANE_ALL) ? v.MANE_ALL : [];
  if (!rows.length) return "";
  const cells = rows.map(r => `
    <tr>
      <td><span class="badge badge-${(r.transcript_type || "").toLowerCase().replace(/_/g, "-")}">${escapeHtml(r.transcript_type || "")}</span></td>
      <td class="mane-tx">${escapeHtml(r.transcript || "")}</td>
      <td>${escapeHtml(r.hgvs_c || "")}</td>
      <td>${escapeHtml(r.hgvs_p || "")}</td>
      <td>${escapeHtml(r.consequence || "")}</td>
    </tr>`).join("");
  return `
    <details class="mane-all">
      <summary>MANE transcripts (${rows.length})</summary>
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

function renderDiseaseList(v, id, withCheckbox) {
  const rows = [];
  const picked = (getEdit(id, "report_diseases") || {});
  for (let i = 1; i <= 5; i++) {
    const d = v[`Disease${i}`];
    if (!d || d === "NA") continue;
    const summary = (String(d).split("\n")[0] || "").slice(0, 120);
    const checked = picked[i] ? "checked" : "";
    const checkbox = withCheckbox
      ? `<input type="checkbox" class="disease-pick" data-id="${escapeAttr(id)}" data-idx="${i}" ${checked} title="報告要發這個疾病" />`
      : "";
    rows.push(`
      <details class="disease-row">
        <summary>${checkbox}<span class="disease-summary-text">${escapeHtml(summary)}</span></summary>
        <div class="disease-detail">${escapeHtml(String(d))}<button type="button" class="disease-collapse">▴ 收合</button></div>
      </details>`);
  }
  if (!rows.length) return "";
  return `<div class="disease-list">${rows.join("")}</div>`;
}

// ---------- Render: sections ---------------------------------------

const REPORT_SECTION_DEFS = [
  { el: "sec-causative", title: "Causative variants", match: id => getStatus(id) === "1", dropdown: "candidate", defaultOpen: true, diseaseCheckbox: true, manualStatus: "1" },
  { el: "sec-other",     title: "Other variants",     match: id => getStatus(id) === "2", dropdown: "candidate", defaultOpen: true, diseaseCheckbox: true, manualStatus: "2" },
  { el: "sec-candidate", title: "Candidate variants", match: id => getStatus(id) === "C", dropdown: "candidate", defaultOpen: true, diseaseCheckbox: true, manualStatus: "C" },
  // ACMG SF / Proactive / Carrier / PharmCat all live inside the
  // Secondary findings collapsible group in the HTML; they render
  // the same way as before, just nested in a different container.
  { el: "sec-acmg-sf",   title: "ACMG SF",            category: "acmg_sf",   dropdown: "panel", diseaseCheckbox: true },
  { el: "sec-proactive", title: "Proactive",          category: "proactive", dropdown: "panel", diseaseCheckbox: true },
  { el: "sec-carrier",   title: "Carrier screening",  category: "carrier",   dropdown: "panel", diseaseCheckbox: true },
];

// Tier sections per 三級輸出計畫.md §2.3. Backend categorises each variant
// into 1A / 1B / 2 / 3 / 4 / 5 based on ClinVar / LOFTEE / ACMG points.
const CANDIDATE_SECTION_DEFS = [
  { el: "cat-tier-1a", title: "1A — ClinVar P/LP ≥ 1★",        category: "1A", dropdown: "candidate", tier: "1A", defaultOpen: true },
  { el: "cat-tier-1b", title: "1B — Frameshift / Nonsense (LOFTEE HC)", category: "1B", dropdown: "candidate", tier: "1B", defaultOpen: true },
  { el: "cat-tier-1c", title: "1C — ACMG points ≥ 4",          category: "1C", dropdown: "candidate", tier: "1C", defaultOpen: true },
  { el: "cat-tier-2",  title: "2 — ClinVar P/LP 0★ or Conflicting (含 P)", category: "2", dropdown: "candidate", tier: "2" },
  { el: "cat-tier-3",  title: "3 — ACMG points < 4",           category: "3",  dropdown: "candidate", tier: "3"  },
  { el: "cat-acmg-sf-c",   title: "ACMG SF",          category: "acmg_sf",   dropdown: "panel" },
  { el: "cat-proactive-c", title: "Proactive",        category: "proactive", dropdown: "panel" },
  { el: "cat-carrier-c",   title: "Carrier screening", category: "carrier",  dropdown: "panel" },
];

function idsForReportSection(def) {
  const known         = Object.keys(state.data.variants || {});
  const reported      = Object.keys(state.reports.status || {});
  const panelReported = Object.keys(state.reports.panels || {});
  const all = Array.from(new Set([...known, ...reported, ...panelReported]));

  if (def.match) {
    // Causative / Other / Candidate report sections: sort by
    // total_score desc, then cluster same-gene variants together.
    // The gene with the highest-scored variant leads; its lower-
    // scored siblings get pulled up directly behind it instead of
    // scattering down the list. Manual entries (no gene_symbol)
    // stay put as singleton clusters.
    const sorted = all.filter(def.match).sort((a, b) => {
      const sa = Number(state.data.variants[a]?.total_score);
      const sb = Number(state.data.variants[b]?.total_score);
      const va = Number.isFinite(sa) ? sa : -Infinity;
      const vb = Number.isFinite(sb) ? sb : -Infinity;
      return vb - va;
    });
    const groups = new Map();
    for (const id of sorted) {
      const key = state.data.variants[id]?.gene_symbol || `__${id}`;
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(id);
    }
    return Array.from(groups.values()).flat();
  }
  if (def.category) {
    const inCat = new Set(state.data.categories?.[def.category] || []);
    return all.filter(id => inCat.has(id) && getPanelStatus(id, def.category) === "V");
  }
  return [];
}

function idsForCandidateSection(def) {
  const ids = state.data.categories?.[def.category] || [];
  // Common-variant filter: gnomAD genome global AF > 0.01 means the
  // variant is too frequent to be causative. Tier 1A keeps these (those
  // are already-curated P/LP★ calls, kept for visibility), every other
  // section drops them. Missing AF is treated as rare.
  if (def.tier === "1A") return ids;
  const COMMON_AF = 0.01;
  return ids.filter(id => {
    const af = Number(state.data.variants?.[id]?.AF);
    return !Number.isFinite(af) || af <= COMMON_AF;
  });
}

function renderBlock(def, ids, openKey) {
  const host = document.getElementById(def.el);
  host.innerHTML = "";
  host.dataset.openKey = openKey;

  const isPanel = def.dropdown === "panel" && !!def.category;
  // Skip X-marked variants — panel sections use panel-specific X, others use global
  const visibleIds = ids.filter(id => isPanel
    ? getPanelStatus(id, def.category) !== "X"
    : getStatus(id) !== "X");
  const wasOpen = toggledBlocks.has(def.el)
    ? host.dataset.wasOpen === "1"
    : !!def.defaultOpen;
  host.dataset.wasOpen = wasOpen ? "1" : "0";

  const header = document.createElement("div");
  header.className = "block-header" + (wasOpen ? " open" : "");
  // Counts: "In panel X / Total Y" so reviewers can see at a glance how
  // much of the section overlaps the requested panel.
  const inPanelCount = visibleIds.filter(
    id => state.data.variants?.[id]?.in_panel
  ).length;
  const countLabel = `In panel ${inPanelCount} / Total ${visibleIds.length}`;
  header.innerHTML = `
    <span><span class="arrow"></span><span class="title">${escapeHtml(def.title)}</span></span>
    <span class="count">${escapeHtml(countLabel)}</span>`;
  host.appendChild(header);

  const body = document.createElement("div");
  body.className = "block-body" + (wasOpen ? " open" : "");
  const manuals = def.manualStatus
    ? (state.reports.manual_variants || []).filter(m => m.status === def.manualStatus)
    : [];
  if (!visibleIds.length && !manuals.length && !def.manualStatus) {
    body.innerHTML = `<div class="muted">（無符合點位）</div>`;
  } else {
    visibleIds.forEach((id, i) => {
      const v = state.data.variants[id] || null;
      body.appendChild(renderVariantCard(v, id, def.dropdown, {
        category: def.category,
        index: i + 1,
        diseaseCheckbox: !!def.diseaseCheckbox,
      }));
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
  for (const def of CANDIDATE_SECTION_DEFS) {
    renderBlock(def, idsForCandidateSection(def), def.el);
  }
  renderPharmcatBlock("cat-pharmcat-c");
  document.getElementById("category-sections").classList.remove("hidden");
  updateInPanelCount();
  renderTierTabBar();
  renderCnvSvTabBar();
  renderMitoTabBar();
}

// Build the SNV/Indel tier tab bar from the same defs / counts as the
// panels themselves. Each tab carries the tier title + 'In panel X /
// Total Y' so the collapsed view still surfaces those numbers. The
// active panel is whatever was active before, falling back to the
// first tier with any visible variants, falling back to 1A.
const TIER_ORDER = ["1A", "1B", "1C", "2", "3"];
let activeTierTab = null;

function renderTierTabBar() {
  const bar = document.getElementById("tier-tab-bar");
  if (!bar) return;

  const counts = {};
  for (const tier of TIER_ORDER) {
    const def = CANDIDATE_SECTION_DEFS.find(d => d.tier === tier);
    const ids = def ? idsForCandidateSection(def) : [];
    const visible = ids.filter(id => getStatus(id) !== "X");
    const inPanel = visible.filter(
      id => state.data.variants?.[id]?.in_panel
    ).length;
    counts[tier] = { total: visible.length, inPanel };
  }

  // Pick the active tier: keep what was active if it still exists,
  // else first tier with variants, else 1A so the bar is never blank.
  if (!TIER_ORDER.includes(activeTierTab)) activeTierTab = null;
  if (!activeTierTab) {
    activeTierTab = TIER_ORDER.find(t => counts[t].total > 0) || "1A";
  }

  const titles = {
    "1A": "1A — ClinVar P/LP ≥ 1★",
    "1B": "1B — Frameshift / Nonsense",
    "1C": "1C — ACMG points ≥ 4",
    "2":  "2 — ClinVar P/LP 0★ or CONF",
    "3":  "3 — ACMG points < 4",
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
  return (state.data?.cnv_variants?.[id])
      || (state.data?.sv_variants?.[id])
      || null;
}

function _cnvSvIdsForTier(tier) {
  const cats = tier.startsWith("CNV-")
    ? state.data?.cnv_categories
    : state.data?.sv_categories;
  return (cats && cats[tier]) || [];
}

function renderCnvSvTabBar() {
  const bar = document.getElementById("cnv-sv-tab-bar");
  if (!bar) return;
  if (!CNV_SV_TIER_ORDER.includes(activeCnvSvTab)) activeCnvSvTab = null;
  if (!activeCnvSvTab) activeCnvSvTab = "CNV-1A";

  const loading = !!state.cnvSvPending;
  bar.innerHTML = CNV_SV_TIER_ORDER.map(t => {
    const active = t === activeCnvSvTab ? " active" : "";
    const ids = _cnvSvIdsForTier(t);
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
    ids.forEach((id, i) => {
      const v = _cnvSvVariantById(id);
      if (!v) return;
      body.appendChild(renderCnvSvCard(v, id, { tier, index: i + 1 }));
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
// Only disease-relevant variants are shown — MITO-1 = pathogenic (per
// MITOMAP/MitoTIP), MITO-2 = anything else with a MITOMAP disease
// association. Polymorphisms / haplogroup variants are dropped server-side.
const MITO_TIER_ORDER = ["MITO-1", "MITO-2"];
const MITO_TITLES = {
  "MITO-1": "1 Pathogenic",
  "MITO-2": "2 Disease-associated",
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
        tier === "MITO-1" ? "（無致病性 mtDNA 變異）" : "（無 disease 相關 mtDNA 變異）"
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

// Click dispatch for all three tier-tab groups (SNV, CNV/SV, Mito).
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
  } else {
    barId = "tier-tab-bar"; current = activeTierTab;
    setActive = t => { activeTierTab = t; }; applyActive = applyTierTabActive;
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
  // for SNVs (gnomAD normalises). MITOMAP search page as a fallback.
  const links = [];
  if (v.REF && v.ALT) {
    links.push({
      label: "gnomAD-MT",
      href: `https://gnomad.broadinstitute.org/variant/M-${v.POS}-${escapeAttr(v.REF)}-${escapeAttr(v.ALT)}?dataset=gnomad_r3`,
    });
  }
  links.push({ label: "MITOMAP", href: "https://www.mitomap.org/foswiki/bin/view/Main/SearchAllele" });
  return links;
}

function _mitoReviewerStatusSel(id) {
  const status = (state.reports?.status?.[id]) || "";
  const opts = ["", "1", "2", "C", "0", "X"];
  return `<select class="status-select" data-id="${escapeAttr(id)}" title="判定">${
    opts.map(s => `<option value="${s}" ${s===status?"selected":""}>${s||"—"}</option>`).join("")
  }</select>`;
}

function _mitoHeteroplasmy(v) {
  const h = v.heteroplasmy;
  if (h == null || !Number.isFinite(Number(h))) return "—";
  return `${(Number(h) * 100).toFixed(1)}%`;
}

function _mitoConsequenceLabel(v) {
  // The MITOMAP-only annotator already emits readable labels
  // (missense / synonymous / stop_gained / non-coding (tRNA) / …).
  return (v.consequence || "—");
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

function _renderMitoDetailBox(v) {
  const filt = v.filter && v.filter !== "." ? v.filter : "PASS";
  const tlod = (v.TLOD != null) ? Number(v.TLOD).toFixed(2) : "—";
  const ad   = v.AD || "—";
  const dp   = (v.depth != null) ? v.depth : "—";
  const consL = _mitoConsequenceLabel(v);
  // MITOMAP reasoning block: disease / status / plasmy / GenBank
  // freq / refs / MitoTIP. Folded — most variants are polymorphisms
  // with nothing here.
  const hasMm = v.mitomap_disease || v.mitomap_status || v.mitomap_plasmy
             || v.mitomap_gb_freq || v.mitotip_score;
  const mmItems = [];
  if (v.mitomap_disease) mmItems.push(`<li><strong>Disease:</strong> ${escapeHtml(v.mitomap_disease)}</li>`);
  if (v.mitomap_status)  mmItems.push(`<li><strong>Status:</strong> ${escapeHtml(v.mitomap_status)}</li>`);
  if (v.mitomap_plasmy)  mmItems.push(`<li><strong>Plasmy (homo/hetero reports):</strong> ${escapeHtml(v.mitomap_plasmy)}</li>`);
  if (v.mitomap_gb_freq) mmItems.push(`<li><strong>GenBank freq FL(CR):</strong> ${escapeHtml(v.mitomap_gb_freq)}${v.mitomap_gb_seqs ? ` <span class="muted">(${escapeHtml(v.mitomap_gb_seqs)} seqs)</span>` : ""}</li>`);
  if (v.mitotip_score)   mmItems.push(`<li><strong>MitoTIP:</strong> ${escapeHtml(v.mitotip_score)}</li>`);
  if (v.mitomap_refs)    mmItems.push(`<li><strong>References:</strong> ${escapeHtml(v.mitomap_refs)}</li>`);
  const mmBlock = hasMm
    ? `<details class="cnv-sv-reasoning" open>
         <summary>MITOMAP${v.mitomap_allele ? ` <span class="muted" style="font-weight:400">(${escapeHtml(v.mitomap_allele)})</span>` : ""}</summary>
         <ul class="cnv-sv-reasoning-list">${mmItems.join("")}</ul>
       </details>`
    : `<div class="muted" style="font-size:12px;margin-top:6px">MITOMAP 無此變異紀錄（多為 polymorphism / haplogroup 變異）</div>`;
  return `<div class="cnv-sv-detail-box">
    <div class="cnv-sv-detail-row">
      <span><strong>變化:</strong> ${escapeHtml(v.REF || "?")}→${escapeHtml(v.ALT || "?")}</span>
      <span><strong>類型:</strong> ${escapeHtml(MITO_LOCUS_LABELS[v.locus_type] || v.locus_type || "—")}</span>
      <span><strong>Heteroplasmy:</strong> ${_mitoHeteroplasmy(v)} <span class="muted">(AD ${escapeHtml(ad)} · DP ${dp})</span></span>
      <span data-tip="${escapeAttr(_mitoFilterTitle(filt))}"><strong>Filter:</strong> ${escapeHtml(filt)} <span class="muted" style="cursor:help">ⓘ</span></span>
    </div>
    <div class="cnv-sv-detail-row">
      <span><strong>Consequence:</strong> ${escapeHtml(consL)}</span>
      ${v.aa_change ? `<span><strong>Protein change:</strong> ${escapeHtml(v.aa_change)}</span>` : ""}
      <span data-tip="${escapeAttr(_MITO_TLOD_TITLE)}"><strong>TLOD:</strong> ${tlod} <span class="muted" style="cursor:help">ⓘ</span></span>
    </div>
    ${mmBlock}
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
  const comment = (getEdit(id, "comment") || "");
  card.innerHTML = `
    <div class="variant-head">
      ${idxTxt}
      ${_mitoReviewerStatusSel(id)}
      ${locusPill}
      <span class="cnv-sv-pos">${escapeHtml(hgvs)}<button class="btn-copy" data-copy="${escapeAttr(hgvs)}" title="複製">${COPY_ICON_SVG}</button>
        ${v.gene_symbol ? ` <span class="muted" style="font-size:12px">${escapeHtml(v.gene_symbol)}</span>` : ""}
      </span>
      <span class="mito-het-badge" title="heteroplasmy fraction">${_mitoHeteroplasmy(v)}</span>
      <span style="flex:1"></span>
      <span class="ext-links">${links}</span>
    </div>
    ${_renderMitoDetailBox(v)}
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
  const options = ["", "1", "2", "C", "0", "X"];
  const statusSel = `<select class="status-select" data-id="${escapeAttr(id)}">${
    options.map(s => `<option value="${s}" ${s===status?"selected":""}>${s||"—"}</option>`).join("")
  }</select>`;
  const idxTxt = opts.index ? `<span class="card-idx">#${opts.index}</span>` : "";
  const links = _cnvSvExternalLinks(v).map(l =>
    `<a href="${escapeAttr(l.href)}" target="_blank" rel="noopener">${escapeHtml(l.label)}</a>`
  ).join("");

  return `<div class="variant-head">
    ${idxTxt}
    ${statusSel}
    <span class="cnv-sv-source-tag">${sourceLabel}</span>
    ${typeChip}
    <span class="cnv-sv-pos">${escapeHtml(region)}<button class="btn-copy" data-copy="${escapeAttr(regionRaw)}" title="複製座標">${COPY_ICON_SVG}</button>
      ${lengthPart ? ` <span class="muted" style="font-size:11px">· ${escapeHtml(lengthPart)}</span>` : ""}
      ${cytoBoth ? ` <span class="muted" style="font-size:11px">· ${escapeHtml(cytoBoth)}</span>` : ""}
    </span>
    <span class="ext-links">${links}</span>
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
    </tr>`;
  };

  const tableHead = `<thead><tr>
    <th></th><th>Gene</th><th>Tx</th><th>Location</th><th>CDS%</th>
    <th>Inheritance</th><th>OMIM</th><th>Phenotype</th><th>Pheno</th>
  </tr></thead>`;
  const visibleRows = genes.map(rowHtml).join("");

  // Overflow body is rendered lazily on first <details> open. For
  // SVs that span 1500+ genes, eagerly building the chip DOM was
  // adding ~100 ms per card even though the panel was hidden.
  const genesOverflow = v.genes_overflow || [];
  const overflowCount = genesOverflow.length + genesCompact.length;
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

function _renderCnvSvSameGeneRow(v) {
  // Mirror the SNV "搜尋同基因" button. AnnotSV gives us a gene_list;
  // pick the first as the search target (typically the primary gene
  // of the SV). Multi-gene SVs can still pivot through the gene-search
  // box in the modal once it's open.
  const genes = Array.isArray(v.gene_list) ? v.gene_list.filter(Boolean) : [];
  if (!genes.length) return "";
  const g = String(genes[0]);
  const hint = genes.length > 1
    ? `搜尋 ${g}（此 SV 共涵蓋 ${genes.length} 個基因；modal 內可改搜其他基因）`
    : `搜尋 ${g} 的所有 SNV/Indel + CNV/SV 變異`;
  return `<div class="variant-badges cnv-sv-same-gene-row">
    <span class="variant-badges-chips"></span>
    <button class="same-gene-btn" data-gene="${escapeAttr(g)}" type="button" title="${escapeAttr(hint)}">搜尋同基因</button>
  </div>`;
}

function renderCnvSvCard(v, id, opts = {}) {
  const card = document.createElement("div");
  card.className = "variant-card cnv-sv-card";
  card.dataset.id = id;
  card.innerHTML = `
    ${_renderCnvSvHeader(v, id, opts)}
    ${_renderCnvSvSameGeneRow(v)}
    ${_renderCnvSvDetailBox(v, id)}
    ${_renderCnvSvGeneTable(v, id)}
    ${_renderCnvSvOverlap(v)}
    ${_renderCnvSvBenign(v)}
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
  const card = t.closest?.(".cnv-sv-card");
  if (!card) return;
  const id = card.dataset.id;
  if (!id) return;
  if (t.matches(".gene-pick")) {
    const picked = { ...(getEdit(id, "report_genes") || {}) };
    const gene = t.dataset.gene;
    if (t.checked) picked[gene] = true; else delete picked[gene];
    setEdit(id, "report_genes", picked);
  } else if (t.matches(".cnv-sv-acmg-select")) {
    setEdit(id, "ACMG_class_sv", t.value);
    // Refresh the sig-* colour class so the field repaints in place
    // without a full card re-render.
    t.classList.remove("sig-p","sig-lp","sig-vus","sig-lb","sig-b");
    const next = SV_ACMG_SIG_CLASS[Number(t.value)];
    if (next) t.classList.add(next);
  }
});

document.addEventListener("input", ev => {
  const t = ev.target;
  if (!t.matches?.(".cnv-sv-comment-text")) return;
  const id = t.dataset.id;
  if (!id) return;
  setEdit(id, "comment", t.value);
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
  const overflowFull = v.genes_overflow || [];
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
    </tr>`;
  };
  const tableHead = `<thead><tr>
    <th></th><th>Gene</th><th>Tx</th><th>Location</th><th>CDS%</th>
    <th>Inheritance</th><th>OMIM</th><th>Phenotype</th><th>Pheno</th>
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

// Tally how many of the currently-loaded SNV variants flagged the in_panel
// bit, then render that next to the "In panel only" toggle so the user
// can tell at a glance whether the filter is doing anything useful.
function updateInPanelCount() {
  const variants = state.data?.variants || {};
  const total = Object.keys(variants).length;
  const inPanel = Object.values(variants).filter(v => v.in_panel).length;
  const el = document.getElementById("in-panel-count");
  if (el) el.textContent = total ? `(${inPanel} / ${total})` : "";
}

function renderPharmcatBlock(hostId) {
  const host = document.getElementById(hostId);
  host.innerHTML = "";
  host.style.display = "";

  const pc = state.data?.pharmcat;
  const genes = pc?.genes ? Object.values(pc.genes) : [];
  genes.sort((a, b) => String(a.gene || "").localeCompare(String(b.gene || "")));

  // Mark the block as data-bearing so the CSS only paints the gray header
  // background when PharmCAT actually returned something.
  host.classList.toggle("has-data", genes.length > 0);

  const wasOpen = toggledBlocks.has(hostId)
    ? host.dataset.wasOpen === "1"
    : false;
  host.dataset.wasOpen = wasOpen ? "1" : "0";

  const header = document.createElement("div");
  header.className = "block-header" + (wasOpen ? " open" : "");
  header.innerHTML = `
    <span><span class="arrow"></span><span class="title">Pharmacogenomics</span></span>
    <span class="count">${genes.length}</span>`;
  host.appendChild(header);

  const body = document.createElement("div");
  body.className = "block-body" + (wasOpen ? " open" : "");
  if (!genes.length) {
    body.innerHTML = `<div class="muted">尚無 PharmCAT 結果。</div>`;
  } else {
    body.innerHTML = renderPharmcatBody(pc, genes);
  }
  host.appendChild(body);

  header.addEventListener("click", () => {
    const open = body.classList.toggle("open");
    header.classList.toggle("open", open);
    host.dataset.wasOpen = open ? "1" : "0";
    toggledBlocks.add(hostId);
  });
}

// PharmCAT genes are split into "Clinically actionable" (the call carries a
// drug-relevant signal) and "Routine" (the call is benign / reference / a
// negative HLA screen / Uncertain Susceptibility / etc.).
//
// Phenotype-based vetoes run BEFORE the drug-rec presence check so that
// e.g. RYR1 (Uncertain Susceptibility) and HLA-B (only *15:02/*57:01/*58:01
// negative) stay in Routine even though PharmCAT lists them as implicated
// for desflurane / allopurinol / etc.
function isPharmcatActionable(gene, recsByGene) {
  const phenotype = String(gene.phenotype || "").trim().toLowerCase();
  const a1 = String(gene.allele1_function || "").trim().toLowerCase();
  const a2 = String(gene.allele2_function || "").trim().toLowerCase();
  const sym = String(gene.gene || "");

  // Vetoes — never actionable regardless of drug rec presence.
  if (phenotype.includes("uncertain susceptibility")) return false;
  if (sym.startsWith("HLA") && !/\bpositive\b/.test(phenotype)) return false;
  if (phenotype.includes("metabolizer") &&
      (phenotype.includes("normal metabolizer") ||
       phenotype.includes("indeterminate") ||
       phenotype.includes("unknown"))) return false;
  // G6PD writes the phenotype as bare "Normal" rather than "Normal Metabolizer";
  // every CPIC drug rec for G6PD Normal is "No reason to avoid based on G6PD
  // status" so it stays in Routine alongside the other Normal Metabolizers.
  if (phenotype === "normal") return false;
  if (phenotype.includes("function") && phenotype.includes("normal function")) return false;
  if (sym === "CFTR" && /non-respons/.test(phenotype)) return false;
  if (sym === "MT-RNR1" && !/increased risk/.test(phenotype)) return false;

  // Positive signals.
  if ((recsByGene[gene.gene] || []).length > 0) return true;
  if (phenotype.includes("metabolizer")) return true;        // already filtered Normal/Indet/Unk above
  if (phenotype.includes("function"))    return true;        // already filtered Normal function
  if (sym.startsWith("HLA"))             return true;        // contains "positive"
  if (sym === "MT-RNR1")                 return true;        // contains "increased risk"
  if (sym === "CFTR")                    return /respons/.test(phenotype);
  if (sym === "VKORC1") {
    return /\baa\b/.test(phenotype)
        || /coumarin sensitivity/.test(a1)
        || /coumarin sensitivity/.test(a2);
  }
  if (phenotype.includes("susceptibility")) return true;     // already filtered Uncertain
  return false;
}

function renderPharmcatBody(pc, genes) {
  const recs = Array.isArray(pc.recommendations) ? pc.recommendations : [];
  const recsByGene = {};
  for (const r of recs) {
    // Older webdata writes implicated as a bare string when there's exactly
    // one gene (jsonlite auto_unbox); normalise so for…of doesn't iterate
    // characters of "CYP3A5".
    const imp = Array.isArray(r.implicated)
      ? r.implicated
      : (r.implicated ? [r.implicated] : []);
    for (const g of imp) {
      if (!recsByGene[g]) recsByGene[g] = [];
      recsByGene[g].push(r);
    }
  }

  const actionable = [];
  const routine = [];
  for (const g of genes) {
    if (isPharmcatActionable(g, recsByGene)) actionable.push(g);
    else routine.push(g);
  }

  const ts = pc.timestamp
    ? `<div class="muted pharmcat-ts">Generated ${escapeHtml(pc.timestamp)}</div>`
    : "";

  const actionableHtml = actionable.length
    ? `<div class="pgx-section">
         <h4 class="pgx-heading">Clinically actionable</h4>
         ${actionable.map(g => renderPharmcatActionableGene(g, recsByGene[g.gene] || [])).join("")}
       </div>`
    : "";

  const routineHtml = routine.length
    ? `<div class="pgx-section pgx-routine">
         <h4 class="pgx-heading">Routine (no special action)</h4>
         ${renderPharmcatRoutine(routine)}
       </div>`
    : "";

  const empty = !actionable.length && !routine.length
    ? `<div class="muted">尚無 PharmCAT 結果。</div>`
    : "";

  return `${ts}${actionableHtml}${routineHtml}${empty}`;
}

function renderPharmcatActionableGene(g, recs) {
  const heading = `${escapeHtml(g.gene || "")} — `
    + `${escapeHtml(displayPhenotype(g))}`
    + (g.label && g.label !== g.phenotype
        ? ` <span class="muted">(${escapeHtml(g.label)})</span>`
        : "");
  // PharmCAT's drugRecommendation text is already valid HTML (with <ul>,
  // <li>, <a>, &quot;, &gt;, …); pass it straight to innerHTML. Source label
  // and drug name are plain strings so they still get escapeHtml.
  const recItems = recs.map(r => {
    const items = (r.items || []).map(it => {
      const src = it.source ? ` <span class="muted">(${escapeHtml(it.source)})</span>` : "";
      return `<div class="pgx-rec-item">${it.text || ""}${src}</div>`;
    }).join("");
    return `<div class="pgx-rec">
              <div class="pgx-drug">${escapeHtml(r.drug)}</div>
              ${items}
            </div>`;
  }).join("");
  const recsBlock = recs.length
    ? `<div class="pgx-recs">${recItems}</div>`
    : `<div class="muted">No drug-level recommendation in PharmCAT report.</div>`;
  return `<div class="pgx-gene">
            <div class="pgx-gene-head">${heading}</div>
            ${recsBlock}
            ${renderPharmcatGeneDetails(g)}
          </div>`;
}

// "n/a" / empty phenotype gets replaced with the allele's clinical-function
// label (so IFNL3 reads "Favorable response allele" instead of "n/a"); G6PD's
// "Normal" is rewritten to "Normal Metabolizer" so it merges into the same
// routine bucket as the CYPs.
function displayPhenotype(g) {
  const raw = String(g.phenotype || "").trim();
  const lower = raw.toLowerCase();
  if (!raw || lower === "n/a" || lower === "—") {
    const fn = (g.allele1_function || g.allele2_function || "").trim();
    if (fn) return fn;
    if (g.label) return g.label;  // CYP4F2 etc.: show diplotype if nothing else
  }
  if (lower === "normal") return "Normal Metabolizer";
  return raw || "—";
}

// HLA rows always get their own line because the diplotype label differs per
// gene; everything else is grouped by displayPhenotype, so e.g. RYR1 and
// CACNA1S (both Uncertain Susceptibility) collapse into a single row.
function renderPharmcatRoutine(genes) {
  const groups   = new Map();   // key (lowercased label) → { label, items }
  const hlaRows  = [];

  for (const g of genes) {
    const sym = String(g.gene || "");
    const label = displayPhenotype(g);
    if (sym.startsWith("HLA")) {
      hlaRows.push({ gene: sym, dip: g.label || "—", phenotype: g.phenotype || "" });
      continue;
    }
    const key = label.toLowerCase();
    if (!groups.has(key)) groups.set(key, { label, items: [] });
    groups.get(key).items.push({ gene: sym, dip: g.label || "" });
  }

  const order = key => {
    if (key.includes("normal metabolizer")) return 0;
    if (key.includes("uncertain"))          return 1;
    return 2;
  };
  const sortedKeys = [...groups.keys()].sort((a, b) => {
    const oa = order(a), ob = order(b);
    return oa !== ob ? oa - ob : a.localeCompare(b);
  });

  const groupRows = sortedKeys.map(key => {
    const { label, items } = groups.get(key);
    const list = items.map(it => {
      const dip = it.dip
        ? ` <span class="muted">(${escapeHtml(it.dip)})</span>`
        : "";
      return `${escapeHtml(it.gene)}${dip}`;
    }).join(", ");
    return `<div class="pgx-routine-row"><strong>${escapeHtml(label)}:</strong> ${list}</div>`;
  }).join("");

  const hlaHtml = hlaRows.map(h => {
    const ph = h.phenotype
      ? ` <span class="muted">(${escapeHtml(h.phenotype)})</span>`
      : "";
    return `<div class="pgx-routine-row"><strong>${escapeHtml(h.gene)}</strong> — ${escapeHtml(h.dip)}${ph}</div>`;
  }).join("");

  return `${groupRows}${hlaHtml}`;
}

function renderPharmcatGeneDetails(g) {
  // Re-uses the per-gene allele / variant detail block from the prior
  // implementation, now nested inside <details> for actionable genes.
  const fnLine = (name, fn) => {
    if (!name && !fn) return "";
    const left  = name ? escapeHtml(name) : "—";
    const right = fn   ? ` <span class="muted">(${escapeHtml(fn)})</span>` : "";
    return `<div>${left}${right}</div>`;
  };
  const vrows = (g.variants || []).map(v => `
    <tr>
      <td>${escapeHtml(v.rsid || "—")}</td>
      <td>${escapeHtml(v.chr || "")}:${escapeHtml(String(v.pos ?? ""))}</td>
      <td>${escapeHtml(v.call || "")}</td>
      <td>${escapeHtml(v.alleles || "")}</td>
    </tr>`).join("");
  const uncalled = (g.uncalled && g.uncalled.length)
    ? `<div class="muted">Uncalled: ${g.uncalled.map(escapeHtml).join(", ")}</div>`
    : "";
  const variants = vrows
    ? `<table class="pharmcat-variants"><thead><tr><th>rsID</th><th>Position</th><th>Call</th><th>Alleles</th></tr></thead><tbody>${vrows}</tbody></table>`
    : "";
  return `<details class="pgx-detail"><summary>詳細</summary>
            <div class="pharmcat-alleles">
              ${fnLine(g.allele1_name, g.allele1_function)}
              ${fnLine(g.allele2_name, g.allele2_function)}
            </div>
            ${uncalled}
            ${variants}
          </details>`;
}

function renderAll() {
  if (!state.data) return;
  renderSampleMeta();
  renderGeneticCounseling();
  renderClinicalDescription();
  renderPhenotype();
  renderVersionPicker();
  renderComment();
  renderReportSections();
  renderCandidateSections();
  updateSaveHint();
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
  } else if (t.matches(".btn-export-html, .btn-export-screening")) {
    // Other export targets (analysis HTML / screening PDF) are not yet
    // ported from the legacy GitHub-Pages tool.
    alert("此匯出格式尚未實作。");
  } else if (t.closest(".clinical-header")) {
    toggleCollapsibleCard(t.closest(".clinical-header"));
  } else if (t.closest(".btn-copy")) {
    ev.stopPropagation();
    copyToClipboard(t.closest(".btn-copy"));
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
  } else if (t.matches(".disease-pick")) {
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

document.addEventListener("change", ev => {
  const t = ev.target;
  if (t.matches(".status-select")) {
    const panel = t.dataset.panel;
    if (panel) setPanelStatus(t.dataset.id, panel, t.value);
    else       setStatus(t.dataset.id, t.value);
  } else if (t.matches("#m-category")) {
    state.reports.category = t.value || null;
    state.dirty = true;
    updateSaveHint();
  } else if (t.matches(".acmg-class")) {
    setEdit(t.dataset.id, "ACMG_classification", t.value);
    updateSaveHint();
    // Re-apply significance color to match the edited value
    t.classList.remove(...SIG_CLASSES);
    const cls = classifySignificance(t.value);
    if (cls) t.classList.add(cls);
  } else if (t.matches(".acmg-score")) {
    setEdit(t.dataset.id, "ACMG_score", t.value);
    updateSaveHint();
  } else if (t.matches(".acmg-crit")) {
    setEdit(t.dataset.id, "ACMG_criteria", t.value);
    updateSaveHint();
  } else if (t.matches(".variant-comment")) {
    setEdit(t.dataset.id, "comment", t.value);
    updateSaveHint();
  } else if (t.matches(".disease-pick")) {
    const picked = { ...(getEdit(t.dataset.id, "report_diseases") || {}) };
    const idx = t.dataset.idx;
    if (t.checked) picked[idx] = true;
    else           delete picked[idx];
    setEdit(t.dataset.id, "report_diseases", picked);
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
  "Likely pathogenic":      "可能致病性",
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
  XL:  "性染色體遺傳",
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
  return String(raw).replace(/_/g, " ");
}
function clinvarText(v) {
  const sig = variantClinSig(v);
  if (!sig) return "此變異位點未在疾病資料庫 (ClinVar) 中報導。";
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
  const isWGS = String(meta.Test || "").toUpperCase() === "WGS";
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
    `2. ClinVar 及 ACMG&AMP 指引: 引用 ClinVar 資料庫截至 ${clinvarDate} 更新的註解，及美國醫學遺傳學暨基因體學學會 (ACMG) 與分子病理學學會 (AMP) 2015 年頒佈的指引，並且主要列入致病 (Pathogenic) 及可能致病 (Likely pathogenic) 變異；其他類別變異經醫師判斷認為與疾病相關時亦可列入。`,
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
  { title: "主動性篩檢基因清單",       url: "./proactive.txt" },
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
    dataVariants?.[id] && getPanelStatus(id, panelKey) === "V"
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
  if (!pc || !pc.genes) {
    w.para("（無 PharmCAT 結果）", { indent: 4 });
    w.gap(4);
    return;
  }
  // Build the same actionable/routine split the web UI uses.
  const recs = Array.isArray(pc.recommendations) ? pc.recommendations : [];
  const recsByGene = {};
  for (const r of recs) {
    const imp = Array.isArray(r.implicated)
      ? r.implicated
      : (r.implicated ? [r.implicated] : []);
    for (const g of imp) {
      if (!recsByGene[g]) recsByGene[g] = [];
      recsByGene[g].push(r);
    }
  }
  const allGenes = Object.values(pc.genes);
  allGenes.sort((a, b) => String(a.gene || "").localeCompare(String(b.gene || "")));
  const actionable = [], routine = [];
  for (const g of allGenes) {
    if (isPharmcatActionable(g, recsByGene)) actionable.push(g);
    else routine.push(g);
  }

  // Actionable
  w.doc.setFont(PDF_FONT_NAME, "bold");
  w.doc.setFontSize(12);
  w.doc.setTextColor(180, 50, 50);
  if (w.y + 14 > w.pageH - w.margin) { w.doc.addPage(); w.y = w.margin; }
  w.doc.text("與用藥相關", w.margin, w.y + 12 * 0.7);
  w.y += 12 * 0.7 + 4;
  w.doc.setTextColor(20, 20, 20);

  for (const g of actionable) {
    const dispPheno = displayPhenotype(g);
    const labelSuffix = (g.label && g.label !== g.phenotype) ? `  (${g.label})` : "";
    w.para(`${g.gene} — ${dispPheno}${labelSuffix}`, { size: 10 });
    const recsList = recsByGene[g.gene] || [];
    if (!recsList.length) {
      w.para("（PharmCAT 報告中尚無藥物層級建議）", { size: 9, indent: 6 });
    } else {
      for (const r of recsList) {
        w.doc.setFont(PDF_FONT_NAME, "normal");
        w.doc.setFontSize(10);
        const dh = 10 * 0.55 + 2;
        if (w.y + dh > w.pageH - w.margin) { w.doc.addPage(); w.y = w.margin; }
        w.doc.setTextColor(70, 70, 70);
        w.doc.text("• " + r.drug, w.margin + 6, w.y + dh * 0.7);
        w.doc.setTextColor(20, 20, 20);
        w.y += dh;
        for (const it of (r.items || [])) {
          // Strip HTML tags from the rec text for the PDF (jsPDF can't render them).
          const txt = String(it.text || "")
            .replace(/<[^>]+>/g, " ")
            .replace(/&quot;/g, '"').replace(/&gt;/g, ">").replace(/&lt;/g, "<")
            .replace(/&amp;/g, "&").replace(/\s+/g, " ").trim();
          if (!txt) continue;
          w.para(`(${it.source || ""}) ${txt}`, { size: 9, indent: 14 });
        }
      }
    }
    w.gap(2);
  }
  if (!actionable.length) w.para("（無）", { size: 9, indent: 6 });
  w.gap(4);

  // Routine — collapse same-phenotype groups, HLA per-row, just like the UI.
  w.doc.setFont(PDF_FONT_NAME, "bold");
  w.doc.setFontSize(12);
  w.doc.setTextColor(60, 90, 60);
  if (w.y + 14 > w.pageH - w.margin) { w.doc.addPage(); w.y = w.margin; }
  w.doc.text("標準處方", w.margin, w.y + 12 * 0.7);
  w.y += 12 * 0.7 + 4;
  w.doc.setTextColor(20, 20, 20);

  const groups = new Map();
  const hlaRows = [];
  for (const g of routine) {
    const sym = String(g.gene || "");
    const label = displayPhenotype(g);
    if (sym.startsWith("HLA")) {
      hlaRows.push({ gene: sym, dip: g.label || "—", phenotype: g.phenotype || "" });
      continue;
    }
    const key = label.toLowerCase();
    if (!groups.has(key)) groups.set(key, { label, items: [] });
    groups.get(key).items.push({ gene: sym, dip: g.label || "" });
  }
  const order = key => key.includes("normal metabolizer") ? 0 : key.includes("uncertain") ? 1 : 2;
  const sortedKeys = [...groups.keys()].sort((a, b) => {
    const oa = order(a), ob = order(b);
    return oa !== ob ? oa - ob : a.localeCompare(b);
  });
  for (const key of sortedKeys) {
    const { label, items } = groups.get(key);
    const list = items.map(it => it.dip ? `${it.gene} (${it.dip})` : it.gene).join(", ");
    w.para(`${label}: ${list}`, { size: 10, indent: 4 });
  }
  for (const h of hlaRows) {
    const ph = h.phenotype ? `  (${h.phenotype})` : "";
    w.para(`${h.gene} — ${h.dip}${ph}`, { size: 10, indent: 4 });
  }
  if (!groups.size && !hlaRows.length) w.para("（無）", { size: 9, indent: 6 });
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
  w.kv("檢驗項目", meta.Test ? `次世代定序${meta.Test === "WGS" ? "全基因組" : "全外顯子"}定序檢測` : "");
  w.kv("產生日期", todayYmd().replace(/^(\d{4})(\d{2})(\d{2})$/, "$1-$2-$3"));
  w.gap(6);

  pdfWriteSection(w, "重大疾病風險篩檢（美國遺傳醫學會 ACMG 次要發現基因）", cats.acmg_sf,   data.variants, "acmg_sf");
  pdfWriteSection(w, "主動性篩檢",                                          cats.proactive, data.variants, "proactive");
  pdfWriteSection(w, "帶因者篩檢",                                          cats.carrier,   data.variants, "carrier");
  pdfWritePharmacogenomics(w, data.pharmcat);

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
    `2. ClinVar 及 ACMG&AMP 指引: 引用 ClinVar 資料庫截至 ${clinvarDate} 更新的註解，及美國醫學遺傳學暨基因體學學會 (ACMG) 與分子病理學學會 (AMP) 2015 年頒佈的指引，並且主要列入致病 (Pathogenic) 及可能致病 (Likely pathogenic) 變異；其他類別變異經醫師判斷認為與疾病相關時亦可列入。`,
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
  try {
    await loadIndex();
    const n = state.index ? state.index.length : 0;
    document.getElementById("search-status").textContent =
      n ? `索引已載入：${n} 筆樣本` : "索引為空";
  } catch (e) {
    document.getElementById("search-status").textContent = "載入索引失敗: " + e.message;
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
  setupInPanelFilter();
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

// "In panel only" toggle: pure presentational filter that hides variant
// cards whose data-in-panel attribute is "false". Re-render is not
// needed because cards are tagged at render time; flipping the checkbox
// just toggles a class on #category-sections that drives a CSS rule.
function setupInPanelFilter() {
  const cb = document.getElementById("filter-in-panel-only");
  if (!cb) return;
  // Apply once at boot so the initial render honours the checked default.
  document.getElementById("category-sections")
    ?.classList.toggle("filter-in-panel-only", cb.checked);
  cb.addEventListener("change", () => {
    document.getElementById("category-sections")
      .classList.toggle("filter-in-panel-only", cb.checked);
  });
}

// OMIM-display toggle: unchecked → hide the .disease-list block under
// every SNV variant card via a class on #category-sections.
function setupOmimFilter() {
  const cb = document.getElementById("filter-omim");
  if (!cb) return;
  const apply = () => document.getElementById("category-sections")
    ?.classList.toggle("hide-omim", !cb.checked);
  apply();
  cb.addEventListener("change", apply);
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
  // itself.
  if (t.matches?.(".modal")) {
    t.classList.add("hidden");
  }
});
document.addEventListener("keydown", ev => {
  if (ev.key === "Escape") {
    document.querySelectorAll(".modal:not(.hidden)").forEach(m => m.classList.add("hidden"));
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

// Editable HPO/panel state for the load-new-case modal. Mirrors
// phenoEdit on the analysis page but kept separate so the analysis
// page's running session isn't disturbed while the modal is open.
const newCaseEdit = {
  hpo: [],
  panels: [],
  emrPhenotype: null,   // raw EMR phenotype payload (read-only ref)
  source: "",           // 'reviewer-txt' / 'EMR' / 'edited' — for the source line
};

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

  // Populate the LIS_ID dropdown from /api/samples/unregistered. Newest
  // first so the just-finished pipeline output sits at the top.
  const select = document.getElementById("new-case-lis-id");
  select.innerHTML = `<option value="">— 載入中… —</option>`;
  showModal("new-case-modal");
  try {
    const list = await apiFetch("/samples/unregistered") || [];
    _unregisteredById = {};
    list.forEach(r => { _unregisteredById[r.lis_id] = r; });
    if (!list.length) {
      select.innerHTML = `<option value="">（沒有未登錄的個案）</option>`;
      return;
    }
    const fmt = ts => ts ? new Date(ts * 1000).toLocaleString() : "";
    const fmtKB = b => `${(b / 1024).toFixed(0)} KB`;
    select.innerHTML = `<option value="">— 選擇 —</option>` +
      list.map(r =>
        `<option value="${escapeAttr(r.lis_id)}">`
        + `${escapeHtml(r.lis_id)}  ·  ${fmt(r.mtime)}`
        + (r.tsv_size ? `  ·  ${fmtKB(r.tsv_size)}` : "")
        + `</option>`
      ).join("");
  } catch (e) {
    select.innerHTML = `<option value="">（讀取失敗：${escapeHtml(e.message)}）</option>`;
  }
});

// Picking a sample preloads phenotype from the reviewer txt + auto-
// fills the MRN that was embedded in the phenotype filename.
document.getElementById("new-case-lis-id")?.addEventListener("change", (ev) => {
  const lis_id = ev.target.value;
  const entry = _unregisteredById[lis_id];
  if (!entry) {
    newCaseEdit.hpo = [];
    newCaseEdit.panels = [];
    newCaseEdit.source = "";
    renderNewCasePhenoEditor();
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
  if (testSel && roster && roster.test_type) testSel.value = roster.test_type;
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
});

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
      return `<li class="chip chip-panel">`
        + `<span class="chip-label">${escapeHtml(p.name || "")}</span>`
        + `<select class="chip-weight" data-nc-panel-idx="${i}" title="Weight">${opts}</select>`
        + `<button type="button" class="chip-remove" data-nc-panel-idx="${i}" title="移除">×</button>`
        + `</li>`;
    }).join("");
  }
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
  drop.classList.remove("hidden");
}
async function _ncRunPanelSearch(q) {
  const drop = document.getElementById("new-case-panel-search-dropdown");
  if (!drop) return;
  q = (q || "").trim();
  let rows = [];
  try { rows = await apiFetch("/panels") || []; }
  catch { rows = []; }
  if (q) {
    const ql = q.toLowerCase();
    rows = rows.filter(r => (r.name || "").toLowerCase().includes(ql));
  }
  rows = rows.slice(0, 15);
  if (!rows.length) { drop.classList.add("hidden"); drop.innerHTML = ""; return; }
  drop.innerHTML = rows.map(r =>
    `<li class="combobox-option" data-nc-panel-pick='${escapeAttr(JSON.stringify(r))}'>`
    + `<span class="opt-name">${escapeHtml(r.name || "")}</span>`
    + (r.n_genes ? `<span class="opt-mrn">${r.n_genes} genes</span>` : "")
    + `</li>`
  ).join("");
  drop.classList.remove("hidden");
}
document.addEventListener("mousedown", ev => {
  const opt = ev.target.closest("[data-nc-hpo-pick], [data-nc-panel-pick]");
  if (!opt) return;
  ev.preventDefault();
  if (opt.dataset.ncHpoPick) {
    const r = JSON.parse(opt.dataset.ncHpoPick);
    const id = r.hpo_id || r.phenotype || "";
    if (!id) return;
    if (!newCaseEdit.hpo.some(h => h.phenotype === id)) {
      newCaseEdit.hpo.push({phenotype: id, label: r.name || id, weight: 1});
      _markNewCaseEdited();
    }
    document.getElementById("new-case-hpo-search").value = "";
    document.getElementById("new-case-hpo-search-dropdown").classList.add("hidden");
    renderNewCasePhenoEditor();
  } else if (opt.dataset.ncPanelPick) {
    const r = JSON.parse(opt.dataset.ncPanelPick);
    const name = r.name || "";
    if (!name) return;
    if (!newCaseEdit.panels.some(p => p.name === name)) {
      newCaseEdit.panels.push({name, weight: 1});
      _markNewCaseEdited();
    }
    document.getElementById("new-case-panel-search").value = "";
    document.getElementById("new-case-panel-search-dropdown").classList.add("hidden");
    renderNewCasePhenoEditor();
  }
});
document.addEventListener("focusout", ev => {
  // Slight delay so click on a dropdown row lands first.
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
  const fd = new FormData(form);
  // Always send the modal-edited chips; backend uses them as the
  // authoritative phenotype (overrides reviewer txt + EMR fallback).
  fd.set("hpo_json",    JSON.stringify(newCaseEdit.hpo || []));
  fd.set("panels_json", JSON.stringify(newCaseEdit.panels || []));
  // run_analysis=true so backend enqueues exomiser/lirical right
  // after register, regardless of whether chips were edited.
  fd.set("run_analysis", "true");
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
    hideModal("new-case-modal");
    await loadIndex();
    await loadSample(out.sample_id);
    renderAll();
    const stEl = document.getElementById("search-status");
    if (stEl) {
      const job = out.job_id ? `；分析已排入 (${out.job_id})` : "";
      stEl.textContent = `已登錄 ${out.sample_id}${job}`;
    }
    if (out.job_id) _startJobPolling(out.sample_id, out.job_id);
  } catch (e) {
    errEl.textContent = e.message;
    errEl.classList.remove("hidden");
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
// Opt in by putting the text on `data-tip` ("\n" → <br>).
(() => {
  let tipEl = null, timer = null, curTarget = null;
  function ensureEl() {
    if (!tipEl) {
      tipEl = document.createElement("div");
      tipEl.className = "app-tooltip";
      tipEl.style.display = "none";
      document.body.appendChild(tipEl);
    }
    return tipEl;
  }
  function show(el) {
    const txt = el.getAttribute("data-tip");
    if (!txt) return;
    const t = ensureEl();
    t.innerHTML = String(txt).split("\n").map(escapeHtml).join("<br>");
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
  function hide() {
    if (timer) { clearTimeout(timer); timer = null; }
    curTarget = null;
    if (tipEl) tipEl.style.display = "none";
  }
  document.addEventListener("mouseover", ev => {
    const el = ev.target.closest("[data-tip]");
    if (!el || el === curTarget) return;
    if (timer) clearTimeout(timer);
    timer = setTimeout(() => { timer = null; show(el); }, 500);
  });
  document.addEventListener("mouseout", ev => {
    const el = ev.target.closest("[data-tip]");
    if (!el) return;
    if (ev.relatedTarget && el.contains(ev.relatedTarget)) return;
    hide();
  });
  document.addEventListener("mousedown", hide, true);
  window.addEventListener("scroll", hide, true);
})();
