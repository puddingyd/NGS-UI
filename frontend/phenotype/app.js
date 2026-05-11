// ============================================================
// 輸入臨床表徵 (HPO / Panel) — standalone tool served at /phenotype/
//
// Originally a GitHub-backed single-page app (hpo-docs/). Reworked
// to talk to the NGS-UI backend instead:
//   • panel list  ← GET  /api/phenotype-tool/panels   (public)
//   • save txt    → POST /api/phenotype-tool/save      (public)
//   • load txt    ← GET  /api/phenotype-tool/load?...  (public)
// No login required — the tool runs on the hospital intranet; only
// the analysis app gates behind auth. Output txt lands in
// NGS_UI/patient_phenotype/ so 載入新個案 picks it up automatically.
// ============================================================

let hpoData = [];
let hpoById = {};            // lookup by HP ID
let fuseInstance = null;
let panelList = [];          // gene panel names
let panelFuse = null;
let generatedContent = "";

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
  loadPanelList();
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
function removeRow(num) { document.getElementById(`row-${num}`)?.remove(); renumberRows(); }
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
    dropdownHighlight = Math.min(dropdownHighlight + 1, items.length - 1);
    updateHighlight(items);
  } else if (event.key === "ArrowUp") {
    event.preventDefault();
    dropdownHighlight = Math.max(dropdownHighlight - 1, 0);
    updateHighlight(items);
  } else if (event.key === "Enter") {
    event.preventDefault();
    if (dropdownHighlight >= 0 && items[dropdownHighlight]) items[dropdownHighlight].click();
  } else if (event.key === "Escape") {
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
  document.getElementById(`selected-${rowNum}`).textContent = `${hpId} ${hpName} (${genes} genes)`;
  document.getElementById(`dropdown-${rowNum}`).classList.remove("visible");
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
function initPanelRows() { for (let i = 0; i < 3; i++) createPanelRow(); }
function addPanelRow() { const r = createPanelRow(); r.querySelector(".search-input")?.focus(); }
function removePanelRow(num) { document.getElementById(`panel-row-${num}`)?.remove(); }

function onPanelSearchInput(rowNum, query) {
  const dropdown = document.getElementById(`panel-dropdown-${rowNum}`);
  const row = document.getElementById(`panel-row-${rowNum}`);
  query = query.trim();
  if (row.dataset.panelName) {
    row.dataset.panelName = "";
    row.querySelector(".search-input").classList.remove("selected");
    document.getElementById(`panel-selected-${rowNum}`).textContent = "";
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
    dropdownHighlight = Math.min(dropdownHighlight + 1, items.length - 1);
    updateHighlight(items);
  } else if (event.key === "ArrowUp") {
    event.preventDefault();
    dropdownHighlight = Math.max(dropdownHighlight - 1, 0);
    updateHighlight(items);
  } else if (event.key === "Enter") {
    event.preventDefault();
    if (dropdownHighlight >= 0 && items[dropdownHighlight]) items[dropdownHighlight].click();
  } else if (event.key === "Escape") {
    dropdown.classList.remove("visible");
  }
}

function selectPanel(rowNum, name) {
  const row = document.getElementById(`panel-row-${rowNum}`);
  const input = row.querySelector(".search-input");
  row.dataset.panelName = name;
  input.value = name;
  input.classList.add("selected");
  document.getElementById(`panel-selected-${rowNum}`).textContent = name;
  document.getElementById(`panel-dropdown-${rowNum}`).classList.remove("visible");
}

// ============================================================
// Load existing phenotype from the server (by LIS_ID then MRN)
// ============================================================

async function loadPatient() {
  const code = document.getElementById("patient-code").value.trim();
  const mrn  = document.getElementById("patient-mrn").value.trim();
  if (!code && !mrn) { showStatus("請先填 LIS_ID 或 MRN。", "error"); return; }
  showStatus("查詢中…", "");
  try {
    const params = new URLSearchParams();
    if (code) params.set("code", code);
    if (mrn)  params.set("mrn", mrn);
    const resp = await fetch(`/api/phenotype-tool/load?${params}`);
    if (resp.status === 404) { showStatus("找不到既有檔案，可以直接開始輸入。", ""); return; }
    if (!resp.ok) { showStatus("讀取失敗。", "error"); return; }
    const body = await resp.json();
    const content = body.content || "";
    const lines = content.trim().split("\n");
    clearAllRows();
    clearAllPanelRows();
    let termCount = 0, panelCount = 0;
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
      } else {
        const row = createPanelRow();
        const num = parseInt(row.id.replace("panel-row-", ""), 10);
        selectPanel(num, col1);
        row.querySelector(".weight-input").value = weight;
        panelCount++;
      }
    }
    // Pad back to the default empty-row count for convenience.
    while (document.querySelectorAll(".phenotype-row:not(.panel-row)").length < 5) createRow();
    while (document.querySelectorAll(".panel-row").length < 3) createPanelRow();
    if (body.code && !code) document.getElementById("patient-code").value = body.code;
    if (body.mrn  && !mrn)  document.getElementById("patient-mrn").value  = body.mrn;
    showStatus(`已載入：${termCount} 個 HPO term、${panelCount} 個 panel（來源 ${body.filename}）`, "success");
  } catch (e) {
    showStatus("讀取失敗：" + (e.message || e), "error");
  }
}

// ============================================================
// Generate + save
// ============================================================

function _collectLines() {
  const lines = ["phenotype\thpo_name\tweight"];
  document.querySelectorAll(".phenotype-row:not(.panel-row)").forEach((row) => {
    const hpId = row.dataset.hpId, hpName = row.dataset.hpName;
    const weight = row.querySelector(".weight-input").value || "1";
    if (hpId && hpName) lines.push(`${hpId}\t${hpName}\t${weight}`);
  });
  document.querySelectorAll(".panel-row").forEach((row) => {
    const panelName = row.dataset.panelName;
    const weight = row.querySelector(".weight-input").value || "1";
    if (panelName) lines.push(`${panelName}\t\t${weight}`);
  });
  return lines;
}

function generateFile() {
  const mrn = document.getElementById("patient-mrn").value.trim();
  if (!mrn) { showStatus("請填病歷號 MRN（必填）。", "error"); return; }
  const lines = _collectLines();
  if (lines.length <= 1) { showStatus("尚未選擇任何 HPO term 或 panel。", "error"); return; }
  generatedContent = lines.join("\n") + "\n";
  document.getElementById("output-preview").style.display = "block";
  document.getElementById("output-content").textContent = generatedContent;
  document.getElementById("btn-save").style.display = "inline-block";
  document.getElementById("btn-download").style.display = "inline-block";
  showStatus("已產生內容，可存到伺服器或下載。", "success");
}

function _filename() {
  const code = document.getElementById("patient-code").value.trim();
  const mrn  = document.getElementById("patient-mrn").value.trim();
  return code ? `${code}_${mrn}_phenotype.txt` : `${mrn}_phenotype.txt`;
}

async function saveToServer() {
  if (!generatedContent) { showStatus("請先按「產生 phenotype.txt」。", "error"); return; }
  const code = document.getElementById("patient-code").value.trim();
  const mrn  = document.getElementById("patient-mrn").value.trim();
  if (!mrn) { showStatus("請填病歷號 MRN（必填）。", "error"); return; }
  showStatus("存檔中…", "");
  try {
    const resp = await fetch("/api/phenotype-tool/save", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ mrn, code: code || "", content: generatedContent }),
    });
    const body = await resp.json().catch(() => ({}));
    if (!resp.ok) throw new Error(body.detail || `${resp.status} ${resp.statusText}`);
    showStatus(`已存到伺服器：${body.path}`, "success");
  } catch (e) {
    showStatus("存檔失敗：" + (e.message || e), "error");
  }
}

function downloadFile() {
  if (!generatedContent) { showStatus("請先按「產生 phenotype.txt」。", "error"); return; }
  const blob = new Blob([generatedContent], { type: "text/plain" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = _filename();
  a.click();
  URL.revokeObjectURL(a.href);
}

function showStatus(msg, type) {
  const el = document.getElementById("status-bar");
  if (!el) return;
  el.textContent = msg;
  el.className = type ? `status-${type}` : "";
}

// ============================================================
// Boot
// ============================================================
loadHPOData();
