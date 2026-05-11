// ============================================================
// HPO Term Manager - Web App
// ============================================================

const CONFIG = {
  REPO_OWNER: "puddingyd",
  REPO_NAME: "hpo-translator",
  REPO_BRANCH: "claude/patient-hpo-term-manager-OB4fk",
  PHENOTYPE_DIR: "patient_phenotype",
  PHENOTYPE_VIP_DIR: "patient_phenotype_VIP",
  // GitHub OAuth
  OAUTH_CLIENT_ID: "Ov23liiYOOhGhhzLYvIy",
  OAUTH_PROXY_URL: "https://hpo-oauth-proxy.puddingyd.workers.dev",
};

let hpoData = [];
let hpoById = {};  // lookup by HP ID
let fuseInstance = null;
let panelList = [];  // gene panel filenames
let panelFuse = null;
let githubToken = localStorage.getItem("github_token") || "";
let generatedContent = "";
let generatedVipContent = "";

// ============================================================
// Data Loading
// ============================================================

async function loadHPOData() {
  const resp = await fetch("hpo_data.json");
  hpoData = await resp.json();

  // Build ID lookup map
  hpoData.forEach((t) => { hpoById[t.id] = t; });

  // Build Fuse.js index
  fuseInstance = new Fuse(hpoData, {
    keys: [
      { name: "name", weight: 2 },
      { name: "syn", weight: 1 },
    ],
    threshold: 0.3,
    distance: 100,
    minMatchCharLength: 2,
    includeScore: true,
  });

  document.getElementById("loading-overlay").classList.add("hidden");
  initRows();
  initPanelRows();
  checkGitHubAuth();
  loadPanelList();
}

// ============================================================
// Phenotype Rows
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
      <input type="text" class="search-input" placeholder="Search HPO term or enter HP number..."
             oninput="onSearchInput(${num}, this.value)"
             onfocus="onSearchFocus(${num})"
             onkeydown="onSearchKeydown(event, ${num})">
      <div class="selected-term" id="selected-${num}"></div>
      <div class="dropdown" id="dropdown-${num}"></div>
    </div>
    <input type="number" class="weight-input" value="1" min="0" step="1" placeholder="W">
    <button class="btn-remove" onclick="removeRow(${num})" title="Remove">&times;</button>
  `;

  container.appendChild(row);
  return row;
}

function initRows() {
  for (let i = 0; i < 5; i++) {
    createRow();
  }
}

function addRow() {
  const row = createRow();
  const input = row.querySelector(".search-input");
  if (input) input.focus();
  return row;
}

function removeRow(num) {
  const row = document.getElementById(`row-${num}`);
  if (row) row.remove();
  renumberRows();
}

function renumberRows() {
  const rows = document.querySelectorAll(".phenotype-row");
  rows.forEach((row, i) => {
    row.querySelector(".row-num").textContent = i + 1;
  });
}

function clearAllRows() {
  document.getElementById("phenotype-rows").innerHTML = "";
  rowCount = 0;
}

// ============================================================
// Search & Dropdown
// ============================================================

let activeDropdownRow = null;
let dropdownHighlight = -1;

function onSearchInput(rowNum, query) {
  const dropdown = document.getElementById(`dropdown-${rowNum}`);
  query = query.trim();

  // Clear selection if user edits
  const row = document.getElementById(`row-${rowNum}`);
  if (row.dataset.hpId) {
    row.dataset.hpId = "";
    row.dataset.hpName = "";
    row.querySelector(".search-input").classList.remove("selected");
    document.getElementById(`selected-${rowNum}`).textContent = "";
  }

  if (query.length < 2) {
    dropdown.classList.remove("visible");
    return;
  }

  let results;

  // Check if input is a number (direct HP ID lookup)
  if (/^\d+$/.test(query)) {
    const num = parseInt(query, 10);
    results = hpoData
      .filter((t) => t.n === num || t.id.includes(query))
      .slice(0, 20)
      .map((t) => ({ item: t }));
  } else {
    results = fuseInstance.search(query, { limit: 20 });
  }

  renderDropdown(rowNum, results);
}

function onSearchFocus(rowNum) {
  activeDropdownRow = rowNum;
  dropdownHighlight = -1;
}

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
    if (dropdownHighlight >= 0 && items[dropdownHighlight]) {
      items[dropdownHighlight].click();
    }
  } else if (event.key === "Escape") {
    dropdown.classList.remove("visible");
  }
}

function updateHighlight(items) {
  items.forEach((item, i) => {
    item.style.background = i === dropdownHighlight ? "#e8f0fe" : "";
  });
  if (items[dropdownHighlight]) {
    items[dropdownHighlight].scrollIntoView({ block: "nearest" });
  }
}

function renderDropdown(rowNum, results) {
  const dropdown = document.getElementById(`dropdown-${rowNum}`);
  dropdownHighlight = -1;

  if (results.length === 0) {
    dropdown.innerHTML = '<div class="dropdown-item" style="color:#999;">No results found</div>';
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

    // Add parent as selectable option
    if (t.par && !seen.has(t.par)) {
      seen.add(t.par);
      const parent = hpoById[t.par];
      if (parent) {
        const pg = parent.g || 0;
        html += `<div class="dropdown-item dropdown-parent" data-hp-id="${t.par}" onclick="selectTerm(${rowNum}, '${t.par}', '${escapeHtml(parent.name)}', ${pg})">
          <span class="hp-parent-arrow">\u2934</span>
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
  const selectedLabel = document.getElementById(`selected-${rowNum}`);
  const dropdown = document.getElementById(`dropdown-${rowNum}`);

  row.dataset.hpId = hpId;
  row.dataset.hpName = hpName;

  input.value = `${hpId} ${hpName}`;
  input.classList.add("selected");
  selectedLabel.textContent = `${hpId} ${hpName} (${genes} genes)`;
  dropdown.classList.remove("visible");
}

function escapeHtml(str) {
  return str.replace(/'/g, "\\'").replace(/"/g, "&quot;");
}

// Close dropdown when clicking outside
document.addEventListener("click", (e) => {
  if (!e.target.closest(".search-wrapper")) {
    document.querySelectorAll(".dropdown").forEach((d) => d.classList.remove("visible"));
  }
});

// ============================================================
// Gene Panel Rows
// ============================================================

let panelRowCount = 0;

async function loadPanelList() {
  // Try loading from a local JSON first
  try {
    const resp = await fetch("panel_list.json");
    if (resp.ok) {
      panelList = await resp.json();
      buildPanelFuse();
      return;
    }
  } catch {}

  // Fallback: fetch from GitHub API
  if (!githubToken) return;
  try {
    const url = `https://api.github.com/repos/${CONFIG.REPO_OWNER}/${CONFIG.REPO_NAME}/contents/data/gene_panels?ref=${CONFIG.REPO_BRANCH}`;
    const resp = await fetch(url, { headers: { Authorization: `Bearer ${githubToken}` } });
    if (resp.ok) {
      const files = await resp.json();
      panelList = files
        .filter((f) => f.name.endsWith(".txt") && f.name !== ".gitkeep")
        .map((f) => f.name.replace(".txt", ""));
      buildPanelFuse();
    }
  } catch {}
}

function buildPanelFuse() {
  panelFuse = new Fuse(panelList.map((name) => ({ name })), {
    keys: ["name"],
    threshold: 0.4,
    distance: 50,
  });
}

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
      <input type="text" class="search-input panel-search" placeholder="Search gene panel..."
             oninput="onPanelSearchInput(${num}, this.value)"
             onfocus="onPanelSearchFocus(${num})"
             onkeydown="onPanelSearchKeydown(event, ${num})">
      <div class="selected-term" id="panel-selected-${num}"></div>
      <div class="dropdown" id="panel-dropdown-${num}"></div>
    </div>
    <input type="number" class="weight-input" value="1" min="0" step="1" placeholder="W">
    <button class="btn-remove" onclick="removePanelRow(${num})" title="Remove">&times;</button>
  `;

  container.appendChild(row);
  return row;
}

function clearAllPanelRows() {
  document.getElementById("panel-rows").innerHTML = "";
  panelRowCount = 0;
}

function initPanelRows() {
  for (let i = 0; i < 3; i++) {
    createPanelRow();
  }
}

function addPanelRow() {
  const row = createPanelRow();
  const input = row.querySelector(".search-input");
  if (input) input.focus();
}

function removePanelRow(num) {
  const row = document.getElementById(`panel-row-${num}`);
  if (row) row.remove();
}

function onPanelSearchInput(rowNum, query) {
  const dropdown = document.getElementById(`panel-dropdown-${rowNum}`);
  const row = document.getElementById(`panel-row-${rowNum}`);
  query = query.trim();

  if (row.dataset.panelName) {
    row.dataset.panelName = "";
    row.querySelector(".search-input").classList.remove("selected");
    document.getElementById(`panel-selected-${rowNum}`).textContent = "";
  }

  if (query.length < 2 || !panelFuse) {
    dropdown.classList.remove("visible");
    return;
  }

  const results = panelFuse.search(query, { limit: 10 });

  if (results.length === 0) {
    dropdown.innerHTML = '<div class="dropdown-item" style="color:#999;">No panels found</div>';
  } else {
    dropdown.innerHTML = results
      .map((r) => `<div class="dropdown-item" data-hp-id="panel" onclick="selectPanel(${rowNum}, '${r.item.name}')">
        <span class="hp-name">${r.item.name}</span>
      </div>`)
      .join("");
  }
  dropdown.classList.add("visible");
}

function onPanelSearchFocus(rowNum) {}

function onPanelSearchKeydown(event, rowNum) {
  const dropdown = document.getElementById(`panel-dropdown-${rowNum}`);
  const items = dropdown.querySelectorAll(".dropdown-item[data-hp-id]");
  // Reuse same keyboard logic
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
  const dropdown = document.getElementById(`panel-dropdown-${rowNum}`);

  row.dataset.panelName = name;
  input.value = name;
  input.classList.add("selected");
  document.getElementById(`panel-selected-${rowNum}`).textContent = name;
  dropdown.classList.remove("visible");
}

// ============================================================
// Load Patient from GitHub
// ============================================================

async function loadPatient() {
  const code = document.getElementById("patient-code").value.trim();
  const mrn = document.getElementById("patient-mrn").value.trim();

  if (!code && !mrn) {
    showStatus("Please enter a patient code or MRN to search.", "error");
    return;
  }

  if (!githubToken) {
    showStatus("Please login to GitHub first.", "error");
    return;
  }

  showStatus("Searching for patient on GitHub...", "success");

  try {
    // List files in patient_phenotype directory
    const dirUrl = `https://api.github.com/repos/${CONFIG.REPO_OWNER}/${CONFIG.REPO_NAME}/contents/${CONFIG.PHENOTYPE_DIR}?ref=${CONFIG.REPO_BRANCH}`;
    const dirResp = await fetch(dirUrl, {
      headers: { Authorization: `Bearer ${githubToken}` },
    });

    if (!dirResp.ok) {
      showStatus("Could not access patient_phenotype directory.", "error");
      return;
    }

    const files = await dirResp.json();
    const query = code || mrn;

    // Find matching file
    const match = files.find((f) => f.name.includes(query) && f.name.endsWith("_phenotype.txt"));

    if (!match) {
      showStatus(`No patient file found matching "${query}".`, "error");
      return;
    }

    // Fetch file content
    const fileResp = await fetch(match.url, {
      headers: { Authorization: `Bearer ${githubToken}`, Accept: "application/vnd.github.raw" },
    });

    if (!fileResp.ok) {
      showStatus("Could not read patient file.", "error");
      return;
    }

    const content = await fileResp.text();
    const lines = content.trim().split("\n");

    // Parse filename: code_mrn_phenotype.txt
    const nameParts = match.name.replace("_phenotype.txt", "").split("_");
    if (nameParts.length >= 2) {
      document.getElementById("patient-code").value = nameParts[0];
      document.getElementById("patient-mrn").value = nameParts[1];
    }

    // Clear existing rows and populate
    clearAllRows();
    clearAllPanelRows();

    let termCount = 0;
    let panelCount = 0;
    for (let i = 1; i < lines.length; i++) {
      const parts = lines[i].split("\t");
      if (parts.length < 2) continue;
      const col1 = parts[0];
      const col2 = parts[1];
      const weight = parts[2] || "1";

      if (col1.startsWith("HP:")) {
        // HPO term row
        const row = createRow();
        const rowNum = rowCount;
        const term = hpoById[col1];
        const genes = term ? (term.g || 0) : 0;

        row.dataset.hpId = col1;
        row.dataset.hpName = col2;
        row.querySelector(".search-input").value = `${col1} ${col2}`;
        row.querySelector(".search-input").classList.add("selected");
        row.querySelector(".weight-input").value = weight;
        document.getElementById(`selected-${rowNum}`).textContent = `${col1} ${col2} (${genes} genes)`;
        termCount++;
      } else if (col1 && !col2) {
        // Panel row (no hpo_name)
        const row = createPanelRow();
        const num = panelRowCount;
        row.dataset.panelName = col1;
        row.querySelector(".search-input").value = col1;
        row.querySelector(".search-input").classList.add("selected");
        row.querySelector(".weight-input").value = weight;
        document.getElementById(`panel-selected-${num}`).textContent = col1;
        panelCount++;
      }
    }

    // Add empty rows to reach minimums
    while (rowCount < 5) createRow();
    while (panelRowCount < 3) createPanelRow();

    renumberRows();
    if (githubToken) {
      document.getElementById("analysis-section").style.display = "block";
    }
    showStatus(`Loaded ${termCount} terms and ${panelCount} panels from ${match.name}`, "success");
  } catch (e) {
    showStatus(`Load error: ${e.message}`, "error");
  }
}

// ============================================================
// Generate Output
// ============================================================

function generateFile() {
  const code = document.getElementById("patient-code").value.trim();
  const mrn = document.getElementById("patient-mrn").value.trim();

  if (!code || !mrn) {
    showStatus("Please enter both patient code and MRN.", "error");
    return;
  }

  const lines = ["phenotype\thpo_name\tweight"];

  // HPO term rows
  document.querySelectorAll(".phenotype-row:not(.panel-row)").forEach((row) => {
    const hpId = row.dataset.hpId;
    const hpName = row.dataset.hpName;
    const weight = row.querySelector(".weight-input").value || "1";
    if (hpId && hpName) {
      lines.push(`${hpId}\t${hpName}\t${weight}`);
    }
  });

  // Panel rows
  document.querySelectorAll(".panel-row").forEach((row) => {
    const panelName = row.dataset.panelName;
    const weight = row.querySelector(".weight-input").value || "1";
    if (panelName) {
      lines.push(`${panelName}\t\t${weight}`);
    }
  });

  if (lines.length <= 1) {
    showStatus("No HPO terms or panels selected.", "error");
    return;
  }

  generatedContent = lines.join("\n") + "\n";
  const filename = `${code}_${mrn}_phenotype.txt`;
  const vipFilename = `${code}_${mrn}_phenotype_VIP.txt`;

  // Build VIP format: "hpo_name HP:XXXXXXX" for HPO terms, panel name for panels
  const vipLines = [];
  document.querySelectorAll(".phenotype-row:not(.panel-row)").forEach((row) => {
    const hpId = row.dataset.hpId;
    const hpName = row.dataset.hpName;
    if (hpId && hpName) {
      vipLines.push(`${hpName} ${hpId}`);
    }
  });
  document.querySelectorAll(".panel-row").forEach((row) => {
    const panelName = row.dataset.panelName;
    if (panelName) {
      vipLines.push(panelName);
    }
  });
  generatedVipContent = vipLines.join("\n") + "\n";

  document.getElementById("output-preview").style.display = "block";
  document.getElementById("output-content").textContent = generatedVipContent;
  document.getElementById("btn-download").style.display = "inline-block";
  document.getElementById("btn-download").onclick = () => downloadFile(filename);

  if (githubToken) {
    document.getElementById("btn-push").style.display = "inline-block";
    document.getElementById("btn-push").onclick = () => pushToGitHub(filename, vipFilename);
  }

  if (githubToken) {
    document.getElementById("analysis-section").style.display = "block";
  }

  showStatus(`Generated ${filename} + ${vipFilename} with ${lines.length - 1} terms.`, "success");
}

function downloadFile(filename) {
  const blob = new Blob([generatedContent], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename || "phenotype.txt";
  a.click();
  URL.revokeObjectURL(url);
}

// ============================================================
// GitHub OAuth
// ============================================================

function githubLogin() {
  if (!CONFIG.OAUTH_CLIENT_ID || !CONFIG.OAUTH_PROXY_URL) {
    const token = prompt(
      "GitHub OAuth not configured.\n\nEnter a GitHub Personal Access Token (with 'repo' scope):\n\nCreate one at: https://github.com/settings/tokens/new"
    );
    if (token) {
      githubToken = token;
      localStorage.setItem("github_token", token);
      checkGitHubAuth();
    }
    return;
  }

  const state = Math.random().toString(36).substring(2);
  sessionStorage.setItem("oauth_state", state);
  const authUrl = `https://github.com/login/oauth/authorize?client_id=${CONFIG.OAUTH_CLIENT_ID}&scope=repo,workflow&state=${state}`;
  window.location.href = authUrl;
}

function setManualToken() {
  const token = prompt(
    "Paste your GitHub Personal Access Token:\n\n" +
    "Create one at: https://github.com/settings/tokens/new\n" +
    "Required scopes: repo, workflow"
  );
  if (token) {
    githubToken = token.trim();
    localStorage.setItem("github_token", githubToken);
    checkGitHubAuth();
  }
}

function githubLogout() {
  githubToken = "";
  localStorage.removeItem("github_token");
  document.getElementById("github-user").style.display = "none";
  document.getElementById("btn-login").style.display = "inline-block";
  document.getElementById("btn-push").style.display = "none";
}

async function checkGitHubAuth() {
  const params = new URLSearchParams(window.location.search);
  if (params.has("code") && params.has("state")) {
    const state = sessionStorage.getItem("oauth_state");
    if (params.get("state") === state && CONFIG.OAUTH_PROXY_URL) {
      try {
        const resp = await fetch(CONFIG.OAUTH_PROXY_URL, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ code: params.get("code") }),
        });
        const data = await resp.json();
        if (data.access_token) {
          githubToken = data.access_token;
          localStorage.setItem("github_token", githubToken);
        }
      } catch (e) {
        console.error("OAuth token exchange failed:", e);
      }
      window.history.replaceState({}, "", window.location.pathname);
    }
  }

  if (!githubToken) return;

  try {
    const resp = await fetch("https://api.github.com/user", {
      headers: { Authorization: `Bearer ${githubToken}` },
    });
    if (resp.ok) {
      const user = await resp.json();
      document.getElementById("github-username").textContent = user.login;
      document.getElementById("github-user").style.display = "inline-flex";
      document.getElementById("btn-login").style.display = "none";
    } else {
      githubLogout();
    }
  } catch {
    githubLogout();
  }
}

// ============================================================
// GitHub Push
// ============================================================

async function pushOneFile(dir, fname, content) {
  const filePath = `${dir}/${fname}`;
  const apiUrl = `https://api.github.com/repos/${CONFIG.REPO_OWNER}/${CONFIG.REPO_NAME}/contents/${filePath}`;

  let sha = null;
  const checkResp = await fetch(`${apiUrl}?ref=${CONFIG.REPO_BRANCH}`, {
    headers: { Authorization: `Bearer ${githubToken}` },
  });
  if (checkResp.ok) {
    const existing = await checkResp.json();
    sha = existing.sha;
  }

  const body = {
    message: `Update ${fname}`,
    content: btoa(unescape(encodeURIComponent(content))),
    branch: CONFIG.REPO_BRANCH,
  };
  if (sha) body.sha = sha;

  const pushResp = await fetch(apiUrl, {
    method: "PUT",
    headers: {
      Authorization: `Bearer ${githubToken}`,
      "Content-Type": "application/json",
    },
    body: JSON.stringify(body),
  });

  if (!pushResp.ok) {
    const err = await pushResp.json();
    throw new Error(`${filePath}: ${err.message}`);
  }
  return filePath;
}

async function pushToGitHub(filename, vipFilename) {
  if (!githubToken) {
    showStatus("Not logged in to GitHub.", "error");
    return;
  }

  try {
    await pushOneFile(CONFIG.PHENOTYPE_DIR, filename, generatedContent);
    await pushOneFile(CONFIG.PHENOTYPE_VIP_DIR, vipFilename, generatedVipContent);

    showStatus(`Successfully pushed ${filename} + ${vipFilename} to GitHub!`, "success");
  } catch (e) {
    showStatus(`GitHub push error: ${e.message}`, "error");
  }
}

// ============================================================
// Status Messages
// ============================================================

function showStatus(msg, type) {
  document.querySelectorAll(".status-msg").forEach((el) => el.remove());

  const div = document.createElement("div");
  div.className = `status-msg ${type}`;
  div.textContent = msg;
  document.getElementById("actions").before(div);

  setTimeout(() => div.remove(), 5000);
}

// ============================================================
// Analysis Trigger (GitHub Actions workflow_dispatch)
// ============================================================

const WORKFLOW_FILE = "analyze.yml";

async function triggerAnalysis(genomeBuild) {
  const code = document.getElementById("patient-code").value.trim();
  const geneList = document.getElementById("gene-list-input").value.trim();
  if (!code) {
    showAnalysisStatus("No patient code.", "error");
    return;
  }
  if (!githubToken) {
    showAnalysisStatus("Not logged in to GitHub.", "error");
    return;
  }

  const btns = document.querySelectorAll(".btn-analysis");
  btns.forEach((b) => (b.disabled = true));
  const cmdPreview = geneList
    ? `run-vcf -${genomeBuild} ${code} '${geneList}'`
    : `run-vcf -${genomeBuild} ${code}`;
  showAnalysisStatus(`Dispatching ${cmdPreview}...`, "pending");

  try {
    const url = `https://api.github.com/repos/${CONFIG.REPO_OWNER}/${CONFIG.REPO_NAME}/actions/workflows/${WORKFLOW_FILE}/dispatches`;
    const resp = await fetch(url, {
      method: "POST",
      headers: {
        Authorization: `Bearer ${githubToken}`,
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        ref: CONFIG.REPO_BRANCH,
        inputs: {
          lis_id: code,
          genome_build: genomeBuild,
          gene_list: geneList,
        },
      }),
    });

    if (resp.status === 204) {
      showAnalysisStatus(`Triggered: ${cmdPreview}`, "success");
    } else {
      const err = await resp.json();
      showAnalysisStatus(`Failed: ${err.message}`, "error");
    }
  } catch (e) {
    showAnalysisStatus(`Error: ${e.message}`, "error");
  } finally {
    btns.forEach((b) => (b.disabled = false));
  }
}

function showAnalysisStatus(msg, type) {
  const el = document.getElementById("analysis-status");
  el.textContent = msg;
  el.style.color = type === "error" ? "#721c24" : type === "success" ? "#155724" : "#856404";
}

// ============================================================
// Terminal (Remote Command)
// ============================================================

const REMOTE_CMD_WORKFLOW = "remote-cmd.yml";

function toggleTerminal() {
  const section = document.getElementById("terminal-section");
  const btn = document.querySelector(".btn-terminal-toggle");
  if (section.style.display === "none") {
    section.style.display = "block";
    btn.innerHTML = "Terminal &#9660;";
    document.getElementById("terminal-input").focus();
  } else {
    section.style.display = "none";
    btn.innerHTML = "Terminal &#9654;";
  }
}

async function sendTerminalCmd() {
  const input = document.getElementById("terminal-input");
  const cmd = input.value.trim();
  if (!cmd) return;
  if (!githubToken) {
    showTerminalStatus("Not logged in to GitHub.", "error");
    return;
  }

  showTerminalStatus(`Sending: ${cmd}`, "pending");

  try {
    const url = `https://api.github.com/repos/${CONFIG.REPO_OWNER}/${CONFIG.REPO_NAME}/actions/workflows/${REMOTE_CMD_WORKFLOW}/dispatches`;
    const resp = await fetch(url, {
      method: "POST",
      headers: {
        Authorization: `Bearer ${githubToken}`,
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        ref: CONFIG.REPO_BRANCH,
        inputs: { command: cmd },
      }),
    });

    if (resp.status === 204) {
      showTerminalStatus(`Sent: ${cmd}`, "success");
      input.value = "";
    } else {
      const err = await resp.json();
      showTerminalStatus(`Failed: ${err.message}`, "error");
    }
  } catch (e) {
    showTerminalStatus(`Error: ${e.message}`, "error");
  }
}

function showTerminalStatus(msg, type) {
  const el = document.getElementById("terminal-status");
  el.textContent = msg;
  el.style.color = type === "error" ? "#721c24" : type === "success" ? "#155724" : "#856404";
}

// ============================================================
// Init
// ============================================================

loadHPOData();
