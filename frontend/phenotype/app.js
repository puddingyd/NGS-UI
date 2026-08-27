// ============================================================
// 輸入臨床表徵 (HPO / Panel) — standalone tool served at /phenotype/
//
// Originally a GitHub-backed single-page app (hpo-docs/). Reworked
// to talk to the NGS-UI backend instead:
//   • panel list  ← GET  /api/phenotype-tool/panels   (public)
//   • save txt    → POST /api/phenotype-tool/save      (public)
//   • load txt    ← GET  /api/phenotype-tool/load?...  (public)
//   • clinical presentation sidecar uses /clinical-presentation/{save,load}
// Phenotype editing is public on the hospital intranet. Documents and GC
// counseling records still use the analysis app's authenticated APIs. Output
// txt lands in
// NGS_UI/patient_phenotype/ so 載入新個案 picks it up automatically.
// ============================================================

let hpoData = [];
let hpoById = {};            // lookup by HP ID
let fuseInstance = null;
let panelList = [];          // gene panel names
let panelFuse = null;
let generatedContent = "";
let currentGeneList = [];
let currentGeneMemberships = null;
let currentGeneListContext = null;
let currentGeneListView = "genes";
let loadedClinicalPresentationSidecar = false;
let loadedPhenotypeSidecar = false;
let clinicalPresentationLastSaved = "";
let clinicalPresentationLastSavedPath = "";
let clinicalAutosaveTimer = null;
let clinicalAutosaveDirty = false;
let clinicalAutosaveInflight = false;
const deadZoneCache = new Map();

// Fixed panel index from /api/phenotype-tool/fixed-panels — WES-I /
// WES-II / WGS tabs are driven entirely by this. Keys look like
// "WES-I__皮膚科__EB" and match a file written into GENE_PANELS_DIR,
// so phenotype_scorer consumes them with no extra plumbing.
let fixedPanelIndex = { series: [] };
const fixedPanelKeys = new Set();          // all known keys (for load() dispatch)
const selectedFixedPanels = new Map();      // key → weight (default 1)

// ============================================================
// Data loading
// ============================================================

async function loadHPOData() {
  const resp = await fetch("hpo_data.json");
  hpoData = await resp.json();
  hpoData.forEach((t) => { hpoById[t.id] = t; });
  fuseInstance = new Fuse(hpoData, {
    keys: [{ name: "name", weight: 2 }, { name: "syn", weight: 1 }],
    threshold: 0.3, distance: 100, minMatchCharLength: 2, includeScore: true,
  });
  document.getElementById("loading-overlay").classList.add("hidden");
  initRows();
  initPanelRows();
  initPanelTabs();
  loadPanelList();
  loadFixedPanels();
}

async function loadPanelList() {
  try {
    const resp = await fetch("/api/phenotype-tool/panels");
    if (resp.ok) {
      const data = await resp.json();
      // Endpoint returns either ["name", ...] or [{name, ...}, ...].
      panelList = (data || []).map(x => (typeof x === "string" ? x : x.name)).filter(Boolean);
      buildPanelFuse();
    }
  } catch { /* offline / not reachable — panel search just stays empty */ }
}

function buildPanelFuse() {
  panelFuse = new Fuse(panelList.map((name) => ({ name })), {
    keys: ["name"], threshold: 0.4, distance: 50,
  });
}

// ============================================================
// HPO term rows
// ============================================================

let rowCount = 0;

function createRow() {
  rowCount++;
  const num = rowCount;
  const container = document.getElementById("phenotype-rows");
  const row = document.createElement("div");
  row.className = "phenotype-row";
  row.id = `row-${num}`;
  row.dataset.hpId = "";
  row.dataset.hpName = "";
  row.innerHTML = `
    <span class="row-num">${num}</span>
    <div class="search-wrapper">
      <input type="text" class="search-input" placeholder="搜尋 HPO 名稱或輸入 HP 數字…"
             oninput="onSearchInput(${num}, this.value)"
             onfocus="onSearchFocus(${num})"
             onkeydown="onSearchKeydown(event, ${num})">
      <div class="selected-term" id="selected-${num}"></div>
      <div class="dropdown" id="dropdown-${num}"></div>
    </div>
    <input type="number" class="weight-input" value="1" min="0" step="1" placeholder="W">
    <button class="btn-remove" onclick="removeRow(${num})" title="移除">&times;</button>
  `;
  container.appendChild(row);
  return row;
}

function initRows() { for (let i = 0; i < 5; i++) createRow(); }
function addRow() { const r = createRow(); r.querySelector(".search-input")?.focus(); return r; }
function removeRow(num) {
  document.getElementById(`row-${num}`)?.remove();
  renumberRows();
  updatePreview();
}
function renumberRows() {
  document.querySelectorAll(".phenotype-row:not(.panel-row)").forEach((row, i) => {
    row.querySelector(".row-num").textContent = i + 1;
  });
}
function clearAllRows() { document.getElementById("phenotype-rows").innerHTML = ""; rowCount = 0; }

// ============================================================
// HPO search & dropdown
// ============================================================

let activeDropdownRow = null;
let dropdownHighlight = -1;

function onSearchInput(rowNum, query) {
  const dropdown = document.getElementById(`dropdown-${rowNum}`);
  query = query.trim();
  const row = document.getElementById(`row-${rowNum}`);
  if (row.dataset.hpId) {
    row.dataset.hpId = ""; row.dataset.hpName = "";
    row.querySelector(".search-input").classList.remove("selected");
    document.getElementById(`selected-${rowNum}`).textContent = "";
    updatePreview();
  }
  if (query.length < 2) { dropdown.classList.remove("visible"); return; }
  let results;
  if (/^\d+$/.test(query)) {
    results = hpoData
      .filter((t) => t.n === parseInt(query, 10) || t.id.includes(query))
      .slice(0, 20).map((t) => ({ item: t }));
  } else {
    results = fuseInstance.search(query, { limit: 20 });
  }
  renderDropdown(rowNum, results);
}

function onSearchFocus(rowNum) { activeDropdownRow = rowNum; dropdownHighlight = -1; }

function onSearchKeydown(event, rowNum) {
  const dropdown = document.getElementById(`dropdown-${rowNum}`);
  const items = dropdown.querySelectorAll(".dropdown-item[data-hp-id]");
  if (event.key === "ArrowDown") {
    event.preventDefault();
    event.stopPropagation();
    dropdownHighlight = Math.min(dropdownHighlight + 1, items.length - 1);
    updateHighlight(items);
  } else if (event.key === "ArrowUp") {
    event.preventDefault();
    event.stopPropagation();
    dropdownHighlight = Math.max(dropdownHighlight - 1, 0);
    updateHighlight(items);
  } else if (event.key === "Enter") {
    event.preventDefault();
    event.stopPropagation();
    const idx = dropdownHighlight >= 0 ? dropdownHighlight : 0;
    if (items[idx]) items[idx].click();
  } else if (event.key === "Escape") {
    event.preventDefault();
    event.stopPropagation();
    dropdown.classList.remove("visible");
  }
}

function updateHighlight(items) {
  items.forEach((item, i) => { item.style.background = i === dropdownHighlight ? "#e8f0fe" : ""; });
  items[dropdownHighlight]?.scrollIntoView({ block: "nearest" });
}

function renderDropdown(rowNum, results) {
  const dropdown = document.getElementById(`dropdown-${rowNum}`);
  dropdownHighlight = -1;
  if (results.length === 0) {
    dropdown.innerHTML = '<div class="dropdown-item" style="color:#999;">找不到符合的</div>';
    dropdown.classList.add("visible");
    return;
  }
  let html = "";
  const seen = new Set();
  results.forEach((r) => {
    const t = r.item;
    if (seen.has(t.id)) return;
    seen.add(t.id);
    const genes = t.g || 0;
    html += `<div class="dropdown-item" data-hp-id="${t.id}" onclick="selectTerm(${rowNum}, '${t.id}', '${escapeHtml(t.name)}', ${genes})">
      <span class="hp-id">${t.id}</span>
      <span class="hp-name">${t.name}</span>
      <span class="hp-genes">(${genes} genes)</span>
    </div>`;
    if (t.par && !seen.has(t.par)) {
      seen.add(t.par);
      const parent = hpoById[t.par];
      if (parent) {
        const pg = parent.g || 0;
        html += `<div class="dropdown-item dropdown-parent" data-hp-id="${t.par}" onclick="selectTerm(${rowNum}, '${t.par}', '${escapeHtml(parent.name)}', ${pg})">
          <span class="hp-parent-arrow">⤴</span>
          <span class="hp-id">${t.par}</span>
          <span class="hp-name">${parent.name}</span>
          <span class="hp-genes">(${pg} genes)</span>
        </div>`;
      }
    }
  });
  dropdown.innerHTML = html;
  dropdown.classList.add("visible");
}

function selectTerm(rowNum, hpId, hpName, genes) {
  const row = document.getElementById(`row-${rowNum}`);
  const input = row.querySelector(".search-input");
  row.dataset.hpId = hpId;
  row.dataset.hpName = hpName;
  input.value = `${hpId} ${hpName}`;
  input.classList.add("selected");
  document.getElementById(`selected-${rowNum}`).innerHTML = `
    <span>${escapeText(hpId)} ${escapeText(hpName)} (${Number(genes) || 0} genes)</span>
    <button type="button" class="gene-list-btn" onclick="openGeneListDrawer('hpo', '${escapeJs(hpId)}', '${escapeJs(`${hpId} ${hpName}`)}')">查看</button>
  `;
  document.getElementById(`dropdown-${rowNum}`).classList.remove("visible");
  updatePreview();
}

function escapeHtml(str) {
  return String(str || "").replace(/'/g, "\\'").replace(/"/g, "&quot;");
}

document.addEventListener("click", (e) => {
  if (!e.target.closest(".search-wrapper")) {
    document.querySelectorAll(".dropdown").forEach((d) => d.classList.remove("visible"));
  }
});

// ============================================================
// Gene panel rows
// ============================================================

let panelRowCount = 0;

function createPanelRow() {
  panelRowCount++;
  const num = panelRowCount;
  const container = document.getElementById("panel-rows");
  const row = document.createElement("div");
  row.className = "phenotype-row panel-row";
  row.id = `panel-row-${num}`;
  row.dataset.panelName = "";
  row.innerHTML = `
    <span class="row-num">P${num}</span>
    <div class="search-wrapper">
      <input type="text" class="search-input panel-search" placeholder="搜尋 gene panel…"
             oninput="onPanelSearchInput(${num}, this.value)"
             onfocus="onPanelSearchFocus(${num})"
             onkeydown="onPanelSearchKeydown(event, ${num})">
      <div class="selected-term" id="panel-selected-${num}"></div>
      <div class="dropdown" id="panel-dropdown-${num}"></div>
    </div>
    <input type="number" class="weight-input" value="1" min="0" step="1" placeholder="W">
    <button class="btn-remove" onclick="removePanelRow(${num})" title="移除">&times;</button>
  `;
  container.appendChild(row);
  return row;
}

function clearAllPanelRows() { document.getElementById("panel-rows").innerHTML = ""; panelRowCount = 0; }
function initPanelRows() { createPanelRow(); }
function addPanelRow() { const r = createPanelRow(); r.querySelector(".search-input")?.focus(); }
function removePanelRow(num) {
  document.getElementById(`panel-row-${num}`)?.remove();
  updatePreview();
}

function onPanelSearchInput(rowNum, query) {
  const dropdown = document.getElementById(`panel-dropdown-${rowNum}`);
  const row = document.getElementById(`panel-row-${rowNum}`);
  query = query.trim();
  if (row.dataset.panelName) {
    row.dataset.panelName = "";
    row.querySelector(".search-input").classList.remove("selected");
    document.getElementById(`panel-selected-${rowNum}`).textContent = "";
    updatePreview();
  }
  if (query.length < 2 || !panelFuse) { dropdown.classList.remove("visible"); return; }
  const results = panelFuse.search(query, { limit: 10 });
  dropdown.innerHTML = results.length === 0
    ? '<div class="dropdown-item" style="color:#999;">找不到 panel</div>'
    : results.map((r) =>
        `<div class="dropdown-item" data-hp-id="panel" onclick="selectPanel(${rowNum}, '${escapeHtml(r.item.name)}')">
          <span class="hp-name">${r.item.name}</span>
        </div>`).join("");
  dropdown.classList.add("visible");
}

function onPanelSearchFocus(rowNum) {}

function onPanelSearchKeydown(event, rowNum) {
  const dropdown = document.getElementById(`panel-dropdown-${rowNum}`);
  const items = dropdown.querySelectorAll(".dropdown-item[data-hp-id]");
  if (event.key === "ArrowDown") {
    event.preventDefault();
    event.stopPropagation();
    dropdownHighlight = Math.min(dropdownHighlight + 1, items.length - 1);
    updateHighlight(items);
  } else if (event.key === "ArrowUp") {
    event.preventDefault();
    event.stopPropagation();
    dropdownHighlight = Math.max(dropdownHighlight - 1, 0);
    updateHighlight(items);
  } else if (event.key === "Enter") {
    event.preventDefault();
    event.stopPropagation();
    const idx = dropdownHighlight >= 0 ? dropdownHighlight : 0;
    if (items[idx]) items[idx].click();
  } else if (event.key === "Escape") {
    event.preventDefault();
    event.stopPropagation();
    dropdown.classList.remove("visible");
  }
}

function selectPanel(rowNum, name) {
  const row = document.getElementById(`panel-row-${rowNum}`);
  const input = row.querySelector(".search-input");
  row.dataset.panelName = name;
  input.value = name;
  input.classList.add("selected");
  document.getElementById(`panel-selected-${rowNum}`).innerHTML = `
    <span>${escapeText(name)}</span>
    <button type="button" class="gene-list-btn" onclick="openGeneListDrawer('panel', '${escapeJs(name)}', '${escapeJs(name)}')">查看</button>
  `;
  document.getElementById(`panel-dropdown-${rowNum}`).classList.remove("visible");
  updatePreview();
}

// ============================================================
// Gene-list drawer
// ============================================================

async function openGeneListDrawer(kind, key, label) {
  const drawer = document.getElementById("gene-list-drawer");
  const backdrop = document.getElementById("gene-list-backdrop");
  const title = document.getElementById("gene-list-title");
  const typeEl = document.getElementById("gene-list-type");
  const sourceEl = document.getElementById("gene-list-source");
  const summaryEl = document.getElementById("gene-list-summary");
  const contentEl = document.getElementById("gene-list-content");
  const filterEl = document.getElementById("gene-list-filter");
  const geneQueryEl = document.getElementById("gene-membership-query");
  if (!drawer || !backdrop) return;
  currentGeneList = [];
  currentGeneMemberships = null;
  currentGeneListContext = { kind, key, label: label || key };
  currentGeneListView = "genes";
  updateDeadZoneButtons();
  drawer.hidden = false;
  backdrop.hidden = false;
  title.textContent = label || key;
  typeEl.textContent = kind === "hpo" ? "HPO term" : "Panel";
  sourceEl.textContent = "";
  summaryEl.textContent = "載入中…";
  contentEl.innerHTML = "";
  if (filterEl) filterEl.value = "";
  if (geneQueryEl) geneQueryEl.value = "";
  try {
    const params = new URLSearchParams({ kind, key });
    const resp = await fetch(`/api/phenotype-tool/gene-list?${params}`);
    const data = await resp.json().catch(() => ({}));
    if (!resp.ok) throw new Error(data.detail || "gene list 載入失敗");
    currentGeneList = Array.isArray(data.genes) ? data.genes : [];
    sourceEl.textContent = data.source ? `source: ${data.source}` : "";
    summaryEl.textContent = `${currentGeneList.length} genes`;
    renderCurrentGeneListView();
  } catch (e) {
    summaryEl.textContent = e.message || String(e);
    contentEl.innerHTML = "";
    currentGeneListContext = null;
    updateDeadZoneButtons();
  }
}

function openGeneMembershipSearch() {
  const drawer = document.getElementById("gene-list-drawer");
  const backdrop = document.getElementById("gene-list-backdrop");
  const title = document.getElementById("gene-list-title");
  const typeEl = document.getElementById("gene-list-type");
  const sourceEl = document.getElementById("gene-list-source");
  const summaryEl = document.getElementById("gene-list-summary");
  const contentEl = document.getElementById("gene-list-content");
  const filterEl = document.getElementById("gene-list-filter");
  const geneQueryEl = document.getElementById("gene-membership-query");
  if (!drawer || !backdrop) return;
  currentGeneList = [];
  currentGeneMemberships = null;
  currentGeneListContext = null;
  currentGeneListView = "genes";
  updateDeadZoneButtons();
  drawer.hidden = false;
  backdrop.hidden = false;
  if (title) title.textContent = "搜尋基因";
  if (typeEl) typeEl.textContent = "Gene lookup";
  if (sourceEl) sourceEl.textContent = "";
  if (summaryEl) summaryEl.textContent = "輸入 gene symbol 後搜尋有哪些 HPO term / panel 包含它。";
  if (contentEl) contentEl.innerHTML = "";
  if (filterEl) filterEl.value = "";
  setTimeout(() => geneQueryEl?.focus(), 0);
}

async function searchGeneMemberships() {
  const geneQueryEl = document.getElementById("gene-membership-query");
  const summaryEl = document.getElementById("gene-list-summary");
  const contentEl = document.getElementById("gene-list-content");
  const q = (geneQueryEl?.value || "").trim();
  if (!q) {
    if (summaryEl) summaryEl.textContent = "請輸入 gene symbol。";
    return;
  }
  currentGeneList = [];
  currentGeneMemberships = null;
  currentGeneListContext = null;
  currentGeneListView = "genes";
  updateDeadZoneButtons();
  if (summaryEl) summaryEl.textContent = "搜尋中…";
  if (contentEl) contentEl.innerHTML = "";
  try {
    const params = new URLSearchParams({ gene: q });
    const resp = await fetch(`/api/phenotype-tool/gene-memberships?${params}`);
    const data = await resp.json().catch(() => ({}));
    if (!resp.ok) throw new Error(data.detail || "gene 搜尋失敗");
    currentGeneMemberships = data;
    renderGeneMemberships(data);
  } catch (e) {
    if (summaryEl) summaryEl.textContent = e.message || String(e);
  }
}

function closeGeneListDrawer() {
  const drawer = document.getElementById("gene-list-drawer");
  const backdrop = document.getElementById("gene-list-backdrop");
  if (drawer) drawer.hidden = true;
  if (backdrop) backdrop.hidden = true;
  currentGeneList = [];
  currentGeneMemberships = null;
  currentGeneListContext = null;
  currentGeneListView = "genes";
  updateDeadZoneButtons();
}

function renderGeneList(genes) {
  const contentEl = document.getElementById("gene-list-content");
  if (!contentEl) return;
  contentEl.classList.remove("dead-zone-mode");
  if (!genes.length) {
    contentEl.innerHTML = '<div class="muted">沒有符合的 gene</div>';
    return;
  }
  contentEl.innerHTML = genes.map((g) => `<span class="gene-pill">${escapeText(g)}</span>`).join("");
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
    return (deadZoneBucket(pctA) - deadZoneBucket(pctB))
      || (pctB - pctA)
      || String(a?.gene || "").localeCompare(String(b?.gene || ""));
  });
}

function geneListFilterQuery() {
  return (document.getElementById("gene-list-filter")?.value || "").trim().toUpperCase();
}

function filteredCurrentGenes() {
  const q = geneListFilterQuery();
  return q ? currentGeneList.filter((g) => String(g).toUpperCase().includes(q)) : currentGeneList;
}

function currentDeadZoneCacheKey(testType) {
  const ctx = currentGeneListContext;
  if (!ctx || !currentGeneList.length) return "";
  return `${testType}:${ctx.kind}:${ctx.key}:${currentGeneList.join(",")}`;
}

function updateDeadZoneButtons() {
  const actions = document.getElementById("gene-list-dead-zone-actions");
  if (!actions) return;
  const show = Boolean(currentGeneListContext && currentGeneList.length);
  actions.hidden = !show;
  actions.querySelectorAll("[data-dead-zone-mode]").forEach((btn) => {
    btn.classList.toggle("is-active", currentGeneListView === btn.dataset.deadZoneMode);
  });
}

function renderCurrentGeneListView() {
  updateDeadZoneButtons();
  if (currentGeneListView === "WES" || currentGeneListView === "WGS") {
    renderDeadZoneList(currentGeneListView);
    return;
  }
  const filtered = filteredCurrentGenes();
  const summaryEl = document.getElementById("gene-list-summary");
  if (summaryEl) {
    const q = geneListFilterQuery();
    summaryEl.textContent = q
      ? `${filtered.length} / ${currentGeneList.length} genes`
      : `${currentGeneList.length} genes`;
  }
  renderGeneList(filtered);
}

async function loadDeadZoneEntries(testType) {
  const cacheKey = currentDeadZoneCacheKey(testType);
  if (!cacheKey) return null;
  if (deadZoneCache.has(cacheKey)) return deadZoneCache.get(cacheKey);
  const resp = await fetch("/api/phenotype-tool/dead-zone", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ test_type: testType, genes: currentGeneList }),
  });
  const data = await resp.json().catch(() => ({}));
  if (!resp.ok) throw new Error(data.detail || "Dead zone 載入失敗");
  deadZoneCache.set(cacheKey, data);
  return data;
}

async function renderDeadZoneList(testType) {
  const contentEl = document.getElementById("gene-list-content");
  const summaryEl = document.getElementById("gene-list-summary");
  if (!contentEl) return;
  const requestKey = currentDeadZoneCacheKey(testType);
  contentEl.classList.add("dead-zone-mode");
  if (summaryEl) summaryEl.textContent = `${testType} dead zone 載入中…`;
  contentEl.innerHTML = "";
  try {
    const data = await loadDeadZoneEntries(testType);
    if (currentGeneListView !== testType || currentDeadZoneCacheKey(testType) !== requestKey) return;
    const threshold = data?.threshold || "";
    let entries = sortDeadZoneEntries(Array.isArray(data?.entries) ? data.entries : []);
    const q = geneListFilterQuery();
    if (q) entries = entries.filter((e) => String(e.gene || "").toUpperCase().includes(q));
    if (summaryEl) {
      summaryEl.textContent = q
        ? `${testType} dead zone: ${entries.length} / ${(data.entries || []).length} genes，threshold ${threshold}X`
        : `${testType} dead zone: ${entries.length} genes，threshold ${threshold}X`;
    }
    if (!entries.length) {
      contentEl.innerHTML = '<div class="muted">目前這組 gene list 沒有 cohort dead-zone 註記。</div>';
      return;
    }
    contentEl.innerHTML = `<ul class="drawer-dead-zone-list">${entries.map((e) => {
      const pct = Number(e.cds_dead_pct || 0);
      const pctLabel = Number.isFinite(pct) ? `${pct.toFixed(1).replace(/\\.0$/, "")}%` : "";
      const label = e.exons_label || (Array.isArray(e.exons) ? e.exons.join(", ") : "");
      return `<li class="${deadZonePctClass(pct)}">
        <span class="dead-zone-gene">${escapeText(e.gene || "")}</span>
        <span class="dead-zone-exons">exon ${escapeText(label)}</span>
        ${pctLabel ? `<span class="dead-zone-cds">CDS ${escapeText(pctLabel)}</span>` : ""}
      </li>`;
    }).join("")}</ul>`;
  } catch (e) {
    if (summaryEl) summaryEl.textContent = e.message || String(e);
    contentEl.innerHTML = "";
  }
}

function renderGeneMemberships(data) {
  const contentEl = document.getElementById("gene-list-content");
  const summaryEl = document.getElementById("gene-list-summary");
  if (!contentEl) return;
  contentEl.classList.remove("dead-zone-mode");
  const hpo = Array.isArray(data.hpo) ? data.hpo : [];
  const panels = Array.isArray(data.panels) ? data.panels : [];
  const hpoTotal = Number(data.hpo_total) || hpo.length;
  const panelTotal = Number(data.panel_total) || panels.length;
  if (summaryEl) {
    summaryEl.textContent = `${data.query || ""} → ${data.canonical_gene || ""}: ${hpoTotal} HPO terms, ${panelTotal} panels`;
  }
  updateDeadZoneButtons();
  if (!hpo.length && !panels.length) {
    contentEl.innerHTML = '<div class="muted">找不到包含這個 gene 的 HPO term 或 panel。</div>';
    return;
  }
  const hpoHtml = hpo.map((item) => `
    <div class="membership-item">
      <div class="membership-item-main">
        <span class="membership-id">${escapeText(item.id)}</span>
        <span class="membership-name">${escapeText(item.name || "")}</span>
      </div>
    </div>
  `).join("");
  const panelHtml = panels.map((item) => `
    <div class="membership-item">
      <div class="membership-item-main">
        <span class="membership-name">${escapeText(item.name)}</span>
        <span class="membership-meta">${Number(item.gene_count) || 0} genes</span>
      </div>
      ${item.source ? `<div class="membership-meta">source: ${escapeText(item.source)}</div>` : ""}
    </div>
  `).join("");
  contentEl.innerHTML = `
    <div class="membership-results">
      <section class="membership-section">
        <h3>HPO terms${hpoTotal > hpo.length ? `（顯示前 ${hpo.length} / ${hpoTotal}）` : ""}</h3>
        <div class="membership-list">${hpoHtml || '<div class="muted">沒有 HPO term</div>'}</div>
      </section>
      <section class="membership-section">
        <h3>Panels${panelTotal > panels.length ? `（顯示前 ${panels.length} / ${panelTotal}）` : ""}</h3>
        <div class="membership-list">${panelHtml || '<div class="muted">沒有 panel</div>'}</div>
      </section>
    </div>
  `;
}

function filterGeneList() {
  if (!currentGeneList.length && currentGeneMemberships) {
    const summaryEl = document.getElementById("gene-list-summary");
    if (summaryEl) summaryEl.textContent = "這個欄位只篩選目前 HPO / panel 的 gene list；gene lookup 結果請用上方搜尋基因。";
    return;
  }
  renderCurrentGeneListView();
}

async function copyGeneList() {
  let text = currentGeneList.join("\n");
  if (!text && currentGeneMemberships) {
    const hpo = (currentGeneMemberships.hpo || []).map((item) =>
      `${item.id}\t${item.name || ""}`.trim());
    const panels = (currentGeneMemberships.panels || []).map((item) => item.name);
    text = [
      `gene\t${currentGeneMemberships.canonical_gene || currentGeneMemberships.query || ""}`,
      "",
      "HPO terms",
      ...hpo,
      "",
      "Panels",
      ...panels,
    ].join("\n");
  }
  if (!text) return;
  try {
    if (navigator.clipboard?.writeText) {
      await navigator.clipboard.writeText(text);
    } else {
      throw new Error("Clipboard API unavailable");
    }
    showStatus("已複製目前清單。", "success");
  } catch {
    const textarea = document.createElement("textarea");
    textarea.value = text;
    textarea.setAttribute("readonly", "");
    textarea.style.position = "fixed";
    textarea.style.left = "-9999px";
    document.body.appendChild(textarea);
    textarea.select();
    const ok = document.execCommand("copy");
    textarea.remove();
    showStatus(ok ? "已複製目前清單。" : "瀏覽器不允許直接複製，請在 gene list 中手動選取。", ok ? "success" : "error");
  }
}

document.getElementById("gene-list-backdrop")?.addEventListener("click", closeGeneListDrawer);
document.getElementById("gene-list-close")?.addEventListener("click", closeGeneListDrawer);
document.getElementById("btn-open-gene-search")?.addEventListener("click", openGeneMembershipSearch);
document.getElementById("gene-membership-search-btn")?.addEventListener("click", searchGeneMemberships);
document.querySelectorAll("[data-dead-zone-mode]").forEach((btn) => {
  btn.addEventListener("click", () => {
    const mode = btn.dataset.deadZoneMode;
    currentGeneListView = currentGeneListView === mode ? "genes" : mode;
    const contentEl = document.getElementById("gene-list-content");
    contentEl?.classList.toggle("dead-zone-mode", currentGeneListView !== "genes");
    renderCurrentGeneListView();
  });
});
document.getElementById("gene-membership-query")?.addEventListener("keydown", (e) => {
  if (e.key === "Enter") {
    e.preventDefault();
    searchGeneMemberships();
  }
});
document.addEventListener("keydown", (e) => {
  if (e.key === "Escape" && !document.getElementById("gc-record-modal")?.hidden) {
    closeGcRecordModal();
    return;
  }
  if (e.key === "Escape" && !document.getElementById("gene-list-drawer")?.hidden) {
    closeGeneListDrawer();
  }
});

// ============================================================
// Fixed panels (WES-I / WES-II / WGS tabs)
// ============================================================

function initPanelTabs() {
  document.querySelectorAll(".panel-tab").forEach((btn) => {
    btn.addEventListener("click", () => {
      const target = btn.dataset.tab;
      document.querySelectorAll(".panel-tab").forEach((b) =>
        b.classList.toggle("is-active", b === btn));
      document.querySelectorAll(".panel-tab-body").forEach((body) =>
        body.classList.toggle("is-active", body.dataset.tabBody === target));
    });
  });
}

async function loadFixedPanels() {
  try {
    const resp = await fetch("/api/phenotype-tool/fixed-panels");
    if (!resp.ok) throw new Error(resp.statusText);
    fixedPanelIndex = await resp.json();
  } catch {
    fixedPanelIndex = { series: [] };
  }
  fixedPanelKeys.clear();
  for (const s of (fixedPanelIndex.series || [])) {
    for (const g of (s.groups || [])) {
      for (const p of (g.panels || [])) fixedPanelKeys.add(p.key);
    }
  }
  renderFixedPanelHosts();
}

function renderFixedPanelHosts() {
  const seriesByKey = {};
  for (const s of (fixedPanelIndex.series || [])) seriesByKey[s.key] = s;
  document.querySelectorAll(".fixed-panel-host").forEach((host) => {
    const skey = host.dataset.series;
    const s = seriesByKey[skey];
    if (!s || !(s.groups || []).length) {
      host.innerHTML = '<div class="muted">尚未匯入此系列的 panel（dev 機執行 <code>scripts/import_fixed_panels.py</code>）</div>';
      return;
    }
    host.innerHTML = s.groups.map((g) => `
      <div class="fp-group">
        <div class="fp-group-title">${escapeText(g.category)}</div>
        <div class="fp-chips">
          ${(g.panels || []).map((p) => `
            <label class="fp-chip" data-key="${escapeAttr(p.key)}" title="${escapeAttr(p.key)}">
              <input type="checkbox" class="fp-chip-cb" value="${escapeAttr(p.key)}">
              <span class="fp-chip-label">${escapeText(p.name)}</span>
              <span class="fp-chip-count">(${p.gene_count || 0})</span>
              <button type="button" class="gene-list-btn fp-gene-list-btn"
                onclick="event.preventDefault(); event.stopPropagation(); openGeneListDrawer('panel', '${escapeJs(p.key)}', '${escapeJs(`${s.key} · ${g.category} · ${p.name}`)}')">查看</button>
            </label>
          `).join("")}
        </div>
      </div>
    `).join("");
  });
  // Wire chip toggles (one listener per checkbox is simplest here).
  document.querySelectorAll(".fp-chip-cb").forEach((cb) => {
    cb.addEventListener("change", () => {
      if (cb.checked) selectedFixedPanels.set(cb.value, 1);
      else selectedFixedPanels.delete(cb.value);
      cb.closest(".fp-chip").classList.toggle("is-selected", cb.checked);
      updatePreview();
    });
  });
  syncFixedPanelUiFromState();
}

function syncFixedPanelUiFromState() {
  document.querySelectorAll(".fp-chip-cb").forEach((cb) => {
    const on = selectedFixedPanels.has(cb.value);
    cb.checked = on;
    cb.closest(".fp-chip").classList.toggle("is-selected", on);
  });
}

function clearSelectedFixedPanels() {
  selectedFixedPanels.clear();
  syncFixedPanelUiFromState();
}

function resetPhenotypeSelections() {
  clearAllRows();
  clearAllPanelRows();
  clearSelectedFixedPanels();
  while (document.querySelectorAll(".phenotype-row:not(.panel-row)").length < 5) createRow();
  while (document.querySelectorAll(".panel-row").length < 1) createPanelRow();
}

function escapeText(s) {
  return String(s || "").replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}
function escapeAttr(s) {
  return String(s || "").replace(/&/g, "&amp;").replace(/"/g, "&quot;").replace(/</g, "&lt;");
}
function escapeJs(s) {
  return String(s || "")
    .replace(/&/g, "&amp;")
    .replace(/"/g, "&quot;")
    .replace(/</g, "&lt;")
    .replace(/\\/g, "\\\\")
    .replace(/'/g, "\\'")
    .replace(/\n/g, "\\n")
    .replace(/\r/g, "");
}

// ============================================================
// Clinical presentation autosave
// ============================================================

function _clinicalPresentationFields() {
  return {
    code: document.getElementById("patient-code").value.trim(),
    mrn: document.getElementById("patient-mrn").value.trim(),
    content: document.getElementById("clinical-presentation-text")?.value || "",
  };
}

async function saveClinicalPresentationSidecar() {
  const { code, mrn, content } = _clinicalPresentationFields();
  if (!code && !mrn) throw new Error("請先填 病歷號 或 檢體編號");
  if (!content.trim() && !loadedClinicalPresentationSidecar) return {};
  if (content === clinicalPresentationLastSaved && loadedClinicalPresentationSidecar) {
    return clinicalPresentationLastSavedPath ? { path: clinicalPresentationLastSavedPath, skipped: true } : {};
  }

  const clinicalResp = await fetch("/api/phenotype-tool/clinical-presentation/save", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ mrn: mrn || "", code: code || "", content }),
  });
  const clinicalBody = await clinicalResp.json().catch(() => ({}));
  if (!clinicalResp.ok) throw new Error(clinicalBody.detail || `${clinicalResp.status} ${clinicalResp.statusText}`);
  loadedClinicalPresentationSidecar = true;
  clinicalPresentationLastSaved = content;
  clinicalPresentationLastSavedPath = clinicalBody.path || clinicalPresentationLastSavedPath;
  return clinicalBody;
}

function scheduleClinicalPresentationAutosave() {
  clinicalAutosaveDirty = true;
  clearTimeout(clinicalAutosaveTimer);
  const { code, mrn } = _clinicalPresentationFields();
  if (!code && !mrn) {
    const msg = "Clinical presentation 尚未自動儲存：請先填病歷號；若沒有病歷號，可先填檢體編號暫存。";
    showStatus(msg, "", { clinical: true });
    return;
  }
  clinicalAutosaveTimer = setTimeout(flushClinicalPresentationAutosave, 1200);
}

async function flushClinicalPresentationAutosave() {
  if (clinicalAutosaveInflight || !clinicalAutosaveDirty) return;
  const { code, mrn } = _clinicalPresentationFields();
  if (!code && !mrn) return;

  clinicalAutosaveDirty = false;
  clinicalAutosaveInflight = true;
  try {
    const contentBeforeSave = document.getElementById("clinical-presentation-text")?.value || "";
    const body = await saveClinicalPresentationSidecar();
    if (body.path && !body.skipped) showStatus(`Clinical presentation 已自動儲存：\n${body.path}`, "success", { clinical: true });
    const contentAfterSave = document.getElementById("clinical-presentation-text")?.value || "";
    if (contentAfterSave !== contentBeforeSave) scheduleClinicalPresentationAutosave();
  } catch (e) {
    clinicalAutosaveDirty = true;
    showStatus("Clinical presentation 自動儲存失敗：" + (e.message || e), "error", { clinical: true });
  } finally {
    clinicalAutosaveInflight = false;
  }
}

function initClinicalPresentationAutosave() {
  document.getElementById("clinical-presentation-text")?.addEventListener("input", scheduleClinicalPresentationAutosave);
  document.getElementById("patient-code")?.addEventListener("input", () => {
    if (clinicalAutosaveDirty || (document.getElementById("clinical-presentation-text")?.value || "").trim()) {
      scheduleClinicalPresentationAutosave();
    }
  });
  document.getElementById("patient-mrn")?.addEventListener("input", () => {
    updatePatientRecordActions();
    if (clinicalAutosaveDirty || (document.getElementById("clinical-presentation-text")?.value || "").trim()) {
      scheduleClinicalPresentationAutosave();
    }
  });
}

// ============================================================
// Patient EMR / genetic-counseling actions
// ============================================================

function emrPatientUrl(mrn) {
  return `http://hisweb.hosp.ncku/Emrquery/autologin.aspx?chartno=${encodeURIComponent(String(mrn || "").trim())}`;
}

function updatePatientRecordActions() {
  const mrn = document.getElementById("patient-mrn")?.value.trim() || "";
  const emr = document.getElementById("btn-emr-link");
  const gc = document.getElementById("btn-gc-records");
  if (emr) {
    if (mrn) {
      emr.href = emrPatientUrl(mrn);
      emr.setAttribute("aria-disabled", "false");
      emr.removeAttribute("tabindex");
    } else {
      emr.removeAttribute("href");
      emr.setAttribute("aria-disabled", "true");
      emr.setAttribute("tabindex", "-1");
    }
  }
  if (gc) gc.disabled = !mrn;
}

let gcRequestSequence = 0;

function closeGcRecordModal() {
  gcRequestSequence += 1;
  const modal = document.getElementById("gc-record-modal");
  if (modal) modal.hidden = true;
  const btn = document.getElementById("btn-gc-records");
  if (btn) btn.textContent = "GC紀錄";
  updatePatientRecordActions();
}

function renderGcRecord(record) {
  const fields = [];
  if (record.reason) {
    fields.push(`<div class="gc-record-field"><div class="gc-record-label">Reason</div><div class="gc-record-text">${escapeText(record.reason)}</div></div>`);
  }
  if (record.diagnosis) {
    fields.push(`<div class="gc-record-field"><div class="gc-record-label">Diagnosis</div><div class="gc-record-text">${escapeText(record.diagnosis)}</div></div>`);
  }
  fields.push(`<div class="gc-record-field"><div class="gc-record-label">Counseling record</div><div class="gc-record-text">${escapeText(record.record || "")}</div></div>`);
  return `<article class="gc-record-card">`
    + `<div class="gc-record-meta"><span>看診日期：${escapeText(record.date_of_consult || "—")}</span></div>`
    + fields.join("")
    + `</article>`;
}

async function openGcRecordModal() {
  const mrn = document.getElementById("patient-mrn")?.value.trim() || "";
  if (!mrn) {
    showStatus("請先填病歷號。", "error");
    return;
  }
  const modal = document.getElementById("gc-record-modal");
  const patient = document.getElementById("gc-record-patient");
  const status = document.getElementById("gc-record-status");
  const list = document.getElementById("gc-record-list");
  const btn = document.getElementById("btn-gc-records");
  if (!modal || !status || !list) return;
  const sequence = ++gcRequestSequence;
  modal.hidden = false;
  if (patient) patient.textContent = `病歷號：${mrn}`;
  status.textContent = "查詢 EMR counseling 紀錄中…";
  status.className = "gc-modal-status";
  list.innerHTML = "";
  const original = btn?.textContent || "GC紀錄";
  if (btn) {
    btn.disabled = true;
    btn.textContent = "查詢中…";
  }
  try {
    const resp = await fetch(`/api/emr/${encodeURIComponent(mrn)}/consultation`, {
      credentials: "same-origin",
    });
    const body = await resp.json().catch(() => ({}));
    if (sequence !== gcRequestSequence) return;
    if (resp.status === 401) {
      status.textContent = "GC 紀錄需要先登入分析系統。";
      status.className = "gc-modal-status status-error";
      list.innerHTML = `<a class="btn btn-secondary" href="/" target="_blank" rel="noopener">開啟分析系統登入</a>`;
      return;
    }
    if (!resp.ok) throw new Error(body.detail || `${resp.status} ${resp.statusText}`);
    const consultation = body.consultation || {};
    if (consultation.error) throw new Error(consultation.error);
    const records = Array.isArray(consultation.records) ? consultation.records : [];
    if (!records.length) {
      status.textContent = "EMR 沒有 GC counseling 紀錄。";
      return;
    }
    status.textContent = `共 ${records.length} 筆紀錄；查詢時間 ${new Date().toLocaleString()}`;
    list.innerHTML = records.map(renderGcRecord).join("");
  } catch (error) {
    if (sequence !== gcRequestSequence) return;
    status.textContent = `GC 紀錄讀取失敗：${error.message || error}`;
    status.className = "gc-modal-status status-error";
  } finally {
    if (sequence === gcRequestSequence && btn) {
      btn.disabled = false;
      btn.textContent = original;
      updatePatientRecordActions();
    }
  }
}

document.getElementById("btn-gc-records")?.addEventListener("click", openGcRecordModal);
document.querySelectorAll("[data-gc-close]").forEach((button) => {
  button.addEventListener("click", closeGcRecordModal);
});
document.getElementById("gc-record-modal")?.addEventListener("click", (event) => {
  if (event.target.id === "gc-record-modal") closeGcRecordModal();
});

// ============================================================
// Load existing phenotype from the server (by LIS_ID then MRN)
// ============================================================

async function loadPatient() {
  const code = document.getElementById("patient-code").value.trim();
  const mrn  = document.getElementById("patient-mrn").value.trim();
  if (!code && !mrn) { showStatus("請先填 LIS_ID 或 MRN。", "error"); return; }
  showStatus("查詢中…", "");
  let loadedClinical = false;
  let loadedPhenotype = false;
  let statusParts = [];
  try {
    const params = new URLSearchParams();
    if (code) params.set("code", code);
    if (mrn)  params.set("mrn", mrn);
    loadedClinicalPresentationSidecar = false;
    loadedPhenotypeSidecar = false;
    clinicalPresentationLastSaved = "";
    clinicalPresentationLastSavedPath = "";
    clinicalAutosaveDirty = false;
    clearTimeout(clinicalAutosaveTimer);
    document.getElementById("clinical-presentation-text").value = "";
    resetPhenotypeSelections();
    updatePreview();
    const clinicalResp = await fetch(`/api/phenotype-tool/clinical-presentation/load?${params}`);
    if (clinicalResp.ok) {
      const clinicalBody = await clinicalResp.json();
      document.getElementById("clinical-presentation-text").value = clinicalBody.content || "";
      loadedClinical = true;
      loadedClinicalPresentationSidecar = true;
      clinicalPresentationLastSaved = clinicalBody.content || "";
      clinicalPresentationLastSavedPath = clinicalBody.path || "";
      statusParts.push(`Clinical presentation（${clinicalBody.filename}）`);
      showInlineClinicalStatus(`已載入 Clinical presentation：\n${clinicalBody.path || clinicalBody.filename}`, "success");
      if (clinicalBody.code && !code) document.getElementById("patient-code").value = clinicalBody.code;
      if (clinicalBody.mrn  && !mrn)  document.getElementById("patient-mrn").value  = clinicalBody.mrn;
    } else if (clinicalResp.status !== 404) {
      showStatus("Clinical presentation 讀取失敗。", "error");
      showInlineClinicalStatus("Clinical presentation 讀取失敗。", "error");
      return;
    } else {
      showInlineClinicalStatus("", "");
    }

    const resp = await fetch(`/api/phenotype-tool/load?${params}`);
    if (resp.ok) {
      const body = await resp.json();
      const content = body.content || "";
      const lines = content.trim().split("\n");
      clearAllRows();
      clearAllPanelRows();
      clearSelectedFixedPanels();
      let termCount = 0, panelCount = 0, fixedCount = 0;
      for (let i = 1; i < lines.length; i++) {
        const parts = lines[i].split("\t");
        if (parts.length < 1 || !parts[0]) continue;
        const col1 = parts[0].trim();
        const col2 = (parts[1] || "").trim();
        const weight = (parts[2] || "1").trim();
        if (col1.startsWith("HP:")) {
          const row = createRow();
          const num = parseInt(row.id.replace("row-", ""), 10);
          const known = hpoById[col1];
          const label = col2 || (known ? known.name : col1);
          const genes = known ? (known.g || 0) : 0;
          selectTerm(num, col1, label, genes);
          row.querySelector(".weight-input").value = weight;
          termCount++;
        } else if (fixedPanelKeys.has(col1)) {
          // Pre-imported WES-I / WES-II / WGS panel → toggle the chip on.
          const w = parseFloat(weight); selectedFixedPanels.set(col1, Number.isFinite(w) ? w : 1);
          fixedCount++;
        } else {
          const row = createPanelRow();
          const num = parseInt(row.id.replace("panel-row-", ""), 10);
          selectPanel(num, col1);
          row.querySelector(".weight-input").value = weight;
          panelCount++;
        }
      }
      syncFixedPanelUiFromState();
      // Pad back to the default empty-row count for convenience.
      while (document.querySelectorAll(".phenotype-row:not(.panel-row)").length < 5) createRow();
      while (document.querySelectorAll(".panel-row").length < 1) createPanelRow();
      if (body.code && !code) document.getElementById("patient-code").value = body.code;
      if (body.mrn  && !mrn)  document.getElementById("patient-mrn").value  = body.mrn;
      updatePatientRecordActions();
      loadedPhenotype = true;
      loadedPhenotypeSidecar = true;
      statusParts.push(`${termCount} 個 HPO term、${fixedCount} 個 fixed panel、${panelCount} 個自由 panel（${body.filename}）`);
    } else if (resp.status !== 404) {
      showStatus("phenotype.txt 讀取失敗。", "error");
      return;
    }
    if (!loadedClinical && !loadedPhenotype) {
      showStatus("找不到既有檔案，可以直接開始輸入。", "");
    } else {
      showStatus(`已載入：${statusParts.join("；")}`, "success");
    }
    updatePatientRecordActions();
    updatePreview();
  } catch (e) {
    showStatus("讀取失敗：" + (e.message || e), "error");
  }
}

// ============================================================
// Custom panel rows
// ============================================================

let customPanelRowCount = 0;

function createCustomPanelRow() {
  customPanelRowCount++;
  const num = customPanelRowCount;
  const container = document.getElementById("custom-panel-rows");
  const row = document.createElement("div");
  row.className = "custom-panel-row";
  row.id = `custom-panel-row-${num}`;
  row.innerHTML = `
    <div class="cp-fields">
      <input type="text" class="cp-name" placeholder="自訂 panel 名稱（例：MyPanel）">
      <input type="text" class="cp-source" placeholder="來源（例：PanelApp / PMID / Lab curated）">
      <textarea class="cp-genes" rows="3" placeholder="基因清單（逗號或換行分隔，例：BRCA1, TP53, C7orf50）"></textarea>
    </div>
    <input type="number" class="weight-input" value="1" min="0" step="1" placeholder="W" title="weight">
    <button class="btn-remove" onclick="removeCustomPanelRow(${num})" title="移除">&times;</button>
  `;
  container.appendChild(row);
  return row;
}
function initCustomPanelRows() { /* none by default — added on demand */ }
function addCustomPanelRow() { const r = createCustomPanelRow(); r.querySelector(".cp-name")?.focus(); }
function removeCustomPanelRow(num) {
  document.getElementById(`custom-panel-row-${num}`)?.remove();
  updatePreview();
}

function _collectCustomPanels() {
  // → [{rowEl, name, genes, weight}] for rows that have a name + genes.
  const out = [];
  document.querySelectorAll(".custom-panel-row").forEach((row) => {
    const name = (row.querySelector(".cp-name")?.value || "").trim();
    const source = (row.querySelector(".cp-source")?.value || "").trim();
    const genesRaw = (row.querySelector(".cp-genes")?.value || "").trim();
    const weight = row.querySelector(".weight-input")?.value || "1";
    if (name && genesRaw) out.push({ rowEl: row, name, source, genes: genesRaw, weight });
    else if (name || source || genesRaw) out.push({ rowEl: row, name, source, genes: genesRaw, weight, incomplete: true });
  });
  return out;
}

// ============================================================
// Generate (creates custom panels + writes phenotype.txt to server)
// ============================================================

function _collectHpoAndPanelLines() {
  const lines = [];
  document.querySelectorAll(".phenotype-row:not(.panel-row)").forEach((row) => {
    const hpId = row.dataset.hpId, hpName = row.dataset.hpName;
    const weight = row.querySelector(".weight-input").value || "1";
    if (hpId && hpName) lines.push(`${hpId}\t${hpName}\t${weight}`);
  });
  // Selected WES-I / WES-II / WGS chips — same wire format as a
  // free-text panel row, the key matches the file in GENE_PANELS_DIR.
  for (const [key, weight] of selectedFixedPanels.entries()) {
    lines.push(`${key}\t\t${weight || 1}`);
  }
  document.querySelectorAll(".panel-row").forEach((row) => {
    const panelName = row.dataset.panelName;
    const weight = row.querySelector(".weight-input").value || "1";
    if (panelName) lines.push(`${panelName}\t\t${weight}`);
  });
  return lines;
}

function _fixedPanelDisplayName(key) {
  for (const s of (fixedPanelIndex.series || [])) {
    for (const g of (s.groups || [])) {
      for (const p of (g.panels || [])) {
        if (p.key === key) return `${s.key} · ${g.category} · ${p.name}`;
      }
    }
  }
  return key;
}

// Preview is reviewer-facing — no tabs, no weights, just the
// human-readable signal: "Seizure HP:0001250" for HPO rows (name
// first, then id), plain "<panel_name>" for panel rows (incl. custom
// panels by their server-sanitised name).
function _collectPreviewLines() {
  const lines = [];
  document.querySelectorAll(".phenotype-row:not(.panel-row)").forEach((row) => {
    const hpId = row.dataset.hpId, hpName = row.dataset.hpName;
    if (hpId && hpName) lines.push(`${hpName} ${hpId}`);
  });
  for (const key of selectedFixedPanels.keys()) {
    lines.push(_fixedPanelDisplayName(key));
  }
  document.querySelectorAll(".panel-row").forEach((row) => {
    const panelName = row.dataset.panelName;
    if (panelName) lines.push(panelName);
  });
  document.querySelectorAll(".custom-panel-row").forEach((row) => {
    const name = (row.querySelector(".cp-name")?.value || "").trim();
    const genes = (row.querySelector(".cp-genes")?.value || "").trim();
    if (name && genes) lines.push(name);
  });
  return lines;
}

function updatePreview() {
  const preview = document.getElementById("output-preview");
  const content = document.getElementById("output-content");
  if (!preview || !content) return;
  const lines = _collectPreviewLines();
  content.textContent = lines.join("\n");
  preview.style.display = lines.length ? "block" : "none";
}

async function generateFile() {
  const saveBtns = document.querySelectorAll(".js-btn-save-phenotype");
  const mrn  = document.getElementById("patient-mrn").value.trim();
  const code = document.getElementById("patient-code").value.trim();
  if (!mrn && !code) { showStatus("請至少填 病歷號 或 檢體編號 其中一個。", "error", { clinical: true }); return; }
  const clinicalPresentation = document.getElementById("clinical-presentation-text")?.value || "";

  const customPanels = _collectCustomPanels();
  const incomplete = customPanels.find(c => c.incomplete);
  if (incomplete) { showStatus("有自訂 panel 列只填了名稱或基因其中一項，請補齊或移除該列。", "error"); return; }

  const baseLines = _collectHpoAndPanelLines();
  if (baseLines.length === 0 && customPanels.length === 0 && selectedFixedPanels.size === 0 && !clinicalPresentation.trim() && !loadedClinicalPresentationSidecar && !loadedPhenotypeSidecar) {
    showStatus("尚未選擇任何 HPO term、panel、自訂 panel，或輸入 Clinical presentation。", "error"); return;
  }

  saveBtns.forEach(btn => { btn.disabled = true; });
  showStatus("處理中…", "", { clinical: true });
  try {
    // 1) Create each custom panel on the server; remember the
    //    server-sanitised name so the phenotype.txt references it
    //    correctly.
    const customLines = [];
    for (const cp of customPanels) {
      const resp = await fetch("/api/phenotype-tool/custom-panel", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ name: cp.name, source: cp.source, genes: cp.genes }),
      });
      const body = await resp.json().catch(() => ({}));
      if (!resp.ok) {
        throw new Error(`自訂 panel「${cp.name}」建立失敗：${body.detail || resp.statusText}`);
      }
      // Reflect the sanitised name back into the row so the user sees it.
      if (cp.rowEl) {
        const nameInp = cp.rowEl.querySelector(".cp-name");
        if (nameInp && nameInp.value.trim() !== body.name) nameInp.value = body.name;
      }
      customLines.push(`${body.name}\t\t${cp.weight}`);
    }

    // 2) Save Clinical presentation sidecar first; the main reviewer UI
    //    also writes back to the same file when it is edited there.
    let clinicalBody = {};
    if (clinicalPresentation.trim() || loadedClinicalPresentationSidecar) {
      clinicalBody = await saveClinicalPresentationSidecar();
      clinicalAutosaveDirty = false;
      clearTimeout(clinicalAutosaveTimer);
    }

    // 3) Build the phenotype.txt body and save it when phenotype content exists.
    let body = {};
    if (baseLines.length || customLines.length || loadedPhenotypeSidecar) {
      const all = ["phenotype\thpo_name\tweight", ...baseLines, ...customLines];
      generatedContent = all.join("\n") + "\n";
      const resp = await fetch("/api/phenotype-tool/save", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ mrn: mrn || "", code: code || "", content: generatedContent }),
      });
      body = await resp.json().catch(() => ({}));
      if (!resp.ok) throw new Error(body.detail || `${resp.status} ${resp.statusText}`);
      loadedPhenotypeSidecar = true;
    } else {
      generatedContent = "";
    }

    // The file on disk keeps the full tab-separated format with
    // weights (parsed by phenotype_io.parse); the on-screen preview
    // strips weights + tabs so reviewers see a clean human-readable
    // list: "HP:0001250 Seizure" / "<panel_name>".
    updatePreview();
    // Refresh the panel autocomplete so the new custom panels show up.
    if (customPanels.length) loadPanelList();
    const cpNote = customPanels.length ? `（含 ${customPanels.length} 個自訂 panel）` : "";
    const savedTargets = [];
    if (body.path) savedTargets.push(body.path);
    if (clinicalBody.path) savedTargets.push(clinicalBody.path);
    const savedMessage = ["已存到伺服器：", ...savedTargets].join("\n") + cpNote;
    showStatus(savedMessage, "success", { clinical: true });
  } catch (e) {
    showStatus(e.message || String(e), "error", { clinical: true });
  } finally {
    saveBtns.forEach(btn => { btn.disabled = false; });
  }
}

function showInlineClinicalStatus(msg, type) {
  const el = document.getElementById("clinical-status-inline");
  if (!el) return;
  el.textContent = msg || "";
  el.className = "inline-status" + (type ? ` status-${type}` : "");
}

function showStatus(msg, type, opts = {}) {
  const el = document.getElementById("status-bar");
  if (el) {
    el.textContent = msg;
    el.className = type ? `status-${type}` : "";
  }
  if (opts.clinical) showInlineClinicalStatus(msg, type);
}

// ============================================================
// Boot
// ============================================================
initClinicalPresentationAutosave();
updatePatientRecordActions();
window.PatientDocuments?.init({
  button: "#btn-patient-documents",
  getContext: () => ({
    mrn: document.getElementById("patient-mrn")?.value.trim() || "",
    sourceSampleId: document.getElementById("patient-code")?.value.trim() || "",
  }),
});
document.querySelector("main")?.addEventListener("input", (event) => {
  if (event.target.closest("#custom-panel-rows")) updatePreview();
});
loadHPOData();
