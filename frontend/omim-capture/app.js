// OMIM 擷取頁：token 管理、擴充功能下載、最近擷取列表（PoC）。
(() => {
  const $ = sel => document.querySelector(sel);
  const esc = value => String(value ?? "").replace(/[&<>"']/g, ch => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[ch]
  ));
  const STATUS = {
    not_in_workbook: ["工作簿沒有此 MIM", "warn"],
    missing_detail: ["只有名稱，缺說明", "warn"],
    outdated: ["OMIM 較新，需要更新", "warn"],
    no_done_date: ["沒有 Done 日期", "info"],
    up_to_date: ["已是最新", "ok"],
    gene_in_workbook: ["基因在工作簿中", "info"],
  };
  const SYNOPSIS = {
    present: ["已取得", "ok"],
    not_loaded: ["未展開", "warn"],
    absent: ["此 entry 沒有", "info"],
    not_applicable: ["不適用", "info"],
  };
  const badge = ([label, cls]) => `<span class="badge ${cls}">${esc(label)}</span>`;
  const fmtTime = iso => {
    if (!iso) return "—";
    const d = new Date(iso);
    return Number.isNaN(d.getTime()) ? iso : d.toLocaleString("zh-TW", { hour12: false });
  };

  async function api(path, options = {}) {
    const response = await fetch(path, { credentials: "same-origin", ...options });
    if (response.status === 401) {
      $("#login-needed").hidden = false;
      throw new Error("not logged in");
    }
    if (!response.ok) throw new Error(`${response.status} ${response.statusText}`);
    return response.json();
  }

  async function loadTokens() {
    const { tokens } = await api("/api/omim/capture-tokens");
    $("#token-table tbody").innerHTML = tokens.length ? tokens.map(t => `
      <tr><td>${esc(t.label || "（未命名）")}</td><td>${esc(fmtTime(t.created_at))}</td>
        <td>${esc(fmtTime(t.last_used_at))}</td>
        <td><button class="btn btn-secondary" data-revoke="${t.id}" type="button">停用</button></td></tr>`).join("")
      : `<tr><td colspan="4" class="muted-text">尚未產生 token。</td></tr>`;
  }

  let currentView = "pending";
  let pendingIds = [];

  async function loadCaptures() {
    if (currentView === "adopted") return loadAdoptions();
    const data = await api("/api/omim/captures?view=pending&limit=500");
    const captures = data.captures || [];
    pendingIds = captures.filter(c => !(c.blockers || []).length).map(c => c.id);
    $("#pending-count").textContent = `（${data.total || 0}）`;
    const adoptAll = $("#btn-adopt-all");
    adoptAll.hidden = !data.db_enabled || !captures.length;
    adoptAll.disabled = !pendingIds.length;
    adoptAll.textContent = pendingIds.length === captures.length
      ? `全部採用（${pendingIds.length}）` : `全部採用（${pendingIds.length}／${captures.length}）`;
    $("#capture-table tbody").innerHTML = captures.length ? captures.map(c => {
      const notes = (c.blockers || []).length ? c.blockers : (c.warnings || []);
      return `
      <tr class="clickable" data-capture="${esc(c.id)}">
        <td>${esc(c.prefix || "")}${esc(c.mim)}</td>
        <td>${esc(c.title)}</td>
        <td>${esc(c.description_source || "—")}</td>
        <td>${badge(SYNOPSIS[c.synopsis_status] || [c.synopsis_status || "—", "info"])}</td>
        <td>${esc(c.omim_edit_date || "—")}</td>
        <td>${esc(fmtTime(c.captured_at))}<br><span class="muted-text">${esc(c.user)}</span></td>
        <td>${(c.blockers || []).length ? badge(["需確認", "warn"]) : ""}${notes.length ? `<div class="muted-text">${notes.map(esc).join("<br>")}</div>` : ""}</td>
      </tr>`;
    }).join("") : `<tr><td colspan="7" class="muted-text">沒有待採用的擷取。</td></tr>`;
  }

  async function loadAdoptions() {
    const { adoptions } = await api("/api/omim/adoptions?limit=500");
    $("#adoption-table tbody").innerHTML = adoptions.length ? adoptions.map(a => `
      <tr class="clickable" data-capture="${esc(a.capture_id)}">
        <td>${esc(fmtTime(a.adopted_at))}</td>
        <td>${esc(a.mim)}</td>
        <td>${esc(a.title || "—")}</td>
        <td>${esc((a.genes || []).join("、") || "—")}${a.reviews ? `<br><span class="muted-text">新增 ${a.reviews} 項待確認</span>` : ""}</td>
        <td>${esc(a.adopted_by)}</td>
      </tr>`).join("") : `<tr><td colspan="5" class="muted-text">還沒有採用紀錄。</td></tr>`;
  }

  function setView(view) {
    currentView = view;
    document.querySelectorAll(".tab").forEach(t => t.classList.toggle("active", t.dataset.view === view));
    $("#capture-table").hidden = view !== "pending";
    $("#adoption-table").hidden = view !== "adopted";
    $("#pending-hint").hidden = view !== "pending";
    $("#btn-adopt-all").style.display = view === "pending" ? "" : "none";
    $("#capture-detail").hidden = true;
    $("#bulk-result").hidden = true;
    loadCaptures();
  }

  async function showCapture(id) {
    const box = $("#capture-detail");
    box.hidden = false;
    box.innerHTML = "載入中…";
    const [mim, name] = id.split("/");
    const { meta, summary } = await api(`/api/omim/captures/${encodeURIComponent(mim)}/${encodeURIComponent(name)}`);
    const cmp = summary.comparison || {};
    const links = (cmp.phenotype_links || []).map(l => `<li>${esc(l.gene)}：Disease${esc((l.slots || []).join(",") || "—")} · ${l.has_detail ? "有說明" : "只有名稱"} · Done ${esc(l.done || "—")} ${badge(STATUS[l.status] || [l.status, "info"])}</li>`).join("");
    const gene = cmp.gene && cmp.gene.in_workbook ? `
      <p><b>基因 ${esc(cmp.gene.gene)}</b>（Done ${esc(cmp.gene.done || "—")}）</p>
      ${list("OMIM 有、工作簿沒有", cmp.gene.added_in_omim, r => `${esc(r.phenotype)} (${esc(r.phenotype_mim)})(${esc(r.inheritance)})`)}
      ${list("工作簿有、OMIM 已不列出", cmp.gene.missing_from_omim, esc)}
      ${list("名稱不同（同 MIM）", cmp.gene.renamed, r => `工作簿：${esc(r.workbook)}<br>OMIM：${esc(r.omim)}`)}
      ${list("工作簿缺 MIM", cmp.gene.workbook_lines_without_mim, esc)}` : "";
    box.innerHTML = `
      <p><b>${esc(summary.prefix)}${esc(summary.mim)} ${esc(summary.title)}</b> ${badge(STATUS[cmp.status] || [cmp.status, "info"])}</p>
      <p class="muted-text">擷取：${esc(fmtTime(meta.captured_at))} · ${esc(meta.user)} · <a href="${esc(meta.url)}" target="_blank" rel="noopener">${esc(meta.url)}</a></p>
      ${list("關係表", summary.relationships, r => `${esc(r.phenotype)} (${esc(r.phenotype_mim || "無 MIM")})(${esc(r.inheritance)}) · ${esc(r.gene)} ${esc(r.gene_mim)}`)}
      ${links ? `<p><b>工作簿</b></p><ul>${links}</ul>` : ""}
      ${gene}
      ${list("注意", summary.warnings, esc)}
      ${(summary.disease_cells || []).map(c => `<p><b>${esc(c.gene)} Disease 欄預覽</b></p><pre>${esc(c.cell)}</pre>`).join("")}
      ${adoptSection(id, summary)}`;
  }

  function adoptSection(id, summary) {
    if (!summary.db_enabled) {
      return `<p class="muted-text">平台 OMIM 資料庫尚未啟用，無法採用。</p>`;
    }
    const preview = summary.adopt_preview || {};
    const lines = (preview.genes || []).map(g => {
      const current = (g.current_labels || []).join("、");
      const change = g.fills_mim_for ? `（補上 MIM，原為「${esc(g.fills_mim_for)}」）`
        : (current && current !== g.label ? `（原為「${esc(current)}」）` : (current ? "" : "（新增）"));
      return `<li>${esc(g.gene)}：${esc(g.label)}${change}</li>`;
    }).join("");
    const note = summary.entry_type === "gene"
      ? "依 OMIM 基因頁同步這個基因的疾病清單；OMIM 已不列出的舊疾病只會列入待確認。"
      : (preview.same_content ? "平台已採用過相同內容。"
        : (preview.text_exists ? "會取代目前的說明與 Clinical Synopsis。" : "會新增說明與 Clinical Synopsis。"));
    return `<p><b>採用</b></p>${lines ? `<ul>${lines}</ul>` : ""}
      <p class="muted-text">${esc(note)}</p>
      <button class="btn btn-primary" type="button" data-adopt="${esc(id)}" ${preview.same_content ? "disabled" : ""}>採用到平台</button>
      <span id="adopt-result" class="muted-text"></span>`;
  }

  async function loadDbStatus() {
    const status = await api("/api/omim/db/status");
    const box = $("#db-status");
    if (!status.enabled) {
      box.innerHTML = `尚未啟用，平台目前讀取 OMIM.xlsx。在 dev 機執行一次<br>
        <code>PYTHONPATH=backend python scripts/omim_db.py init</code><br>即可從 OMIM.xlsx 建立資料庫並啟用採用功能。`;
      $("#btn-export").hidden = true;
      return;
    }
    const c = status.counts || {};
    box.innerHTML = `已啟用（平台讀取資料庫）：基因 ${c.genes}、疾病 ${c.diseases}、已採用擷取的疾病說明 ${c.phenotype_texts}、待確認 ${c.open_reviews}、異動紀錄 ${c.changes}。
      <br>建立時間：${esc(fmtTime((status.meta || {}).created_at))}`;
    $("#btn-export").hidden = false;
  }

  const STEP_KIND = { gene_page: ["基因頁", "info"], missing: ["缺說明", "warn"], stale: ["逾一年", "info"] };
  const STATUS_BADGE = {
    in_progress: ["進行中", "warn"], gene_sync: ["基因頁待同步", "info"], phenotype: ["疾病待補", "info"],
  };
  let queueOffset = 0;

  function reviewHtml(r) {
    const d = r.detail || {};
    const label = r.current_label || d.label || d.disease_cell || d.omim_disease || "";
    const canRemove = r.disease_id && r.disease_active;
    return `<li data-review-item="${r.id}">${esc(r.kind_label)}：${esc(label)}
      ${d.suggestion ? `<br><span class="muted-text">可能是：${esc(d.suggestion)}</span>` : ""}
      <span class="review-actions"><button class="btn btn-secondary" type="button" data-review="${r.id}" data-action="keep">保留</button>${
        canRemove ? `<button class="btn btn-secondary" type="button" data-review="${r.id}" data-action="remove">移除疾病</button>` : ""}</span></li>`;
  }

  async function loadQueue(append = false) {
    if (!append) queueOffset = 0;
    const status = $("#queue-status").value;
    const scope = $("#queue-scope").value;
    const data = await api(`/api/omim/queue?status=${encodeURIComponent(status)}&scope=${encodeURIComponent(scope)}&offset=${queueOffset}&limit=100`);
    const c = data.counts || {};
    $("#queue-total").textContent = `（${data.total} 個基因）`;
    $("#queue-summary").textContent = `進行中 ${c.in_progress || 0}、基因頁待同步 ${c.gene_sync || 0}（其中疾病相關 ${c.panel_gene_sync || 0}）、只有疾病待補 ${c.phenotype || 0}；待補疾病頁 ${c.phenotype_steps || 0}、待確認 ${c.reviews || 0} 項。`;
    const rows = data.items.map(g => `
      <tr>
        <td>${g.url ? `<a href="${esc(g.url)}" target="_blank" rel="noopener"><b>${esc(g.gene || g.gene_mim)}</b></a>` : `<b>${esc(g.gene)}</b>`}
          <br><span class="muted-text">${esc(g.gene_mim)}${g.panel ? " · 疾病相關" : ""}</span>
          <br><button class="btn btn-secondary btn-mini" type="button" data-lookup="${esc(g.gene_mim || g.gene)}">查看</button></td>
        <td>${badge(STATUS_BADGE[g.status] || [g.status, "info"])}${g.synced_at ? `<br><span class="muted-text">基因頁同步 ${esc(fmtTime(g.synced_at))}</span>` : ""}${g.partially_updated ? `<br><span class="muted-text">已有疾病被更新（共用疾病）</span>` : ""}</td>
        <td><ol class="step-list">${g.steps.map(st => `<li>
          <a href="${esc(st.url)}" target="_blank" rel="noopener">${esc(st.mim)}</a>
          ${st.type === "gene" ? "" : esc((st.labels || [])[0] || "")}
          ${(st.kinds || []).map(k => badge(STEP_KIND[k] || [k, "info"])).join(" ")}
          ${st.captured ? badge(["已擷取待採用", "ok"]) : ""}${st.recently_adopted ? badge(["7 天內採用過", "info"]) : ""}</li>`).join("")}</ol></td>
        <td>${g.reviews.length ? `<ul class="review-list">${g.reviews.map(reviewHtml).join("")}</ul>` : ""}</td>
        <td>${g.reported ? esc(g.reported) : ""}</td>
      </tr>`).join("");
    const body = $("#queue-table tbody");
    if (append) body.insertAdjacentHTML("beforeend", rows);
    else body.innerHTML = rows || `<tr><td colspan="5" class="muted-text">沒有待補項目。</td></tr>`;
    queueOffset += data.items.length;
    $("#btn-queue-more").hidden = queueOffset >= data.total;
  }

  // ------------------------------------------------------- gene lookup
  const SOURCE_BADGE = {
    capture: ["擷取採用", "ok"], xlsx: ["OMIM.xlsx 原文", "info"], label: ["只有名稱 *", "warn"],
  };
  const ACTION_LABEL = {
    adopt_text: "採用說明", update_label: "更新疾病名稱", add_disease: "新增疾病",
    add_gene: "新增基因", remove_disease: "移除疾病",
  };
  let lookupTimer = null;

  let lookupItems = [];
  let lookupActive = -1;

  function hideLookupDropdown() {
    $("#lookup-dropdown").classList.add("hidden");
    lookupActive = -1;
  }

  function renderLookupDropdown() {
    const list = $("#lookup-dropdown");
    if (!lookupItems.length) { hideLookupDropdown(); return; }
    list.innerHTML = lookupItems.map((i, n) => `
      <li class="combobox-option${n === lookupActive ? " active" : ""}" role="option" data-index="${n}">
        <span class="opt-symbol">${esc(i.symbol)}</span>
        <span class="opt-name">${i.symbol !== i.gene ? `舊名稱／別名，現為 ${esc(i.gene)}` : ""}</span>
        <span class="opt-mim">${esc(i.gene_mim)}</span>
      </li>`).join("");
    list.classList.remove("hidden");
  }

  async function lookupSuggest() {
    const q = $("#lookup-q").value.trim();
    if (q.length < 1 || /^\d+$/.test(q)) { lookupItems = []; hideLookupDropdown(); return; }
    const { items } = await api(`/api/omim/genes/search?q=${encodeURIComponent(q)}&limit=15`);
    if ($("#lookup-q").value.trim() !== q) return;  // typed on meanwhile
    lookupItems = items;
    lookupActive = -1;
    renderLookupDropdown();
  }

  function pickLookup(index) {
    const item = lookupItems[index];
    if (!item) return;
    hideLookupDropdown();
    runLookup(item.symbol);
  }

  function geneButton(g) {
    return `<button class="btn btn-secondary" type="button" data-lookup="${esc(g.gene_mim || g.gene || g.symbol)}">${esc(g.gene || g.symbol)} ${esc(g.gene_mim || "")}${g.diseases != null ? `（${g.diseases} 個疾病）` : ""}</button>`;
  }

  function renderGene(g) {
    const q = g.queue;
    const statusBadge = q ? badge(STATUS_BADGE[q.status] || [q.status, "info"]) : badge(["已完成，不在待補清單", "ok"]);
    const diseases = g.diseases.map((d, i) => `
      <tr>
        <td>${i + 1}</td>
        <td>${d.url ? `<a href="${esc(d.url)}" target="_blank" rel="noopener">${esc(d.phenotype_mim)}</a> ` : '<span class="muted-text">無 MIM</span> '}${esc(d.label)}
          ${d.linked_genes.length ? `<br><span class="muted-text">也連到：${d.linked_genes.map(x => `<a href="#" data-lookup="${esc(x)}">${esc(x)}</a>`).join("、")}</span>` : ""}</td>
        <td>${badge(SOURCE_BADGE[d.source] || [d.source, "info"])}
          ${d.adopted_at ? `<br><span class="muted-text">${esc(fmtTime(d.adopted_at))} · ${esc(d.adopted_by)}${d.omim_edit_date ? `<br>OMIM 編輯 ${esc(d.omim_edit_date)}` : ""}</span>` : ""}</td>
        <td>${d.reported ? esc(d.reported) : ""}</td>
      </tr>
      <tr class="lookup-text-row"><td></td><td colspan="3"><details><summary>卡片／報告使用的內容</summary><pre>${esc(d.text)}</pre></details></td></tr>`).join("");
    const steps = q && q.steps.length ? `<p><b>待補步驟</b></p><ol class="step-list">${q.steps.map(st => `<li>
      <a href="${esc(st.url)}" target="_blank" rel="noopener">${esc(st.mim)}</a> ${st.type === "gene" ? "基因頁" : esc((st.labels || [])[0] || "")}
      ${(st.kinds || []).map(k => badge(STEP_KIND[k] || [k, "info"])).join(" ")}</li>`).join("")}</ol>` : "";
    const reviews = q && q.reviews.length ? `<p><b>待確認</b></p><ul class="review-list">${q.reviews.map(reviewHtml).join("")}</ul>` : "";
    const removed = g.removed.length ? `<p><b>已移除的疾病</b></p><ul>${g.removed.map(r =>
      `<li>${esc(r.label)} <span class="muted-text">${esc(fmtTime(r.updated_at))}</span></li>`).join("")}</ul>` : "";
    const changes = g.changes.length ? `<details class="lookup-changes"><summary>異動紀錄（${g.changes.length}）</summary>
      <table class="data-table"><thead><tr><th>時間</th><th>使用者</th><th>動作</th><th>內容</th></tr></thead><tbody>${
        g.changes.map(ch => `<tr><td>${esc(fmtTime(ch.at))}</td><td>${esc(ch.user)}</td><td>${esc(ACTION_LABEL[ch.action] || ch.action)}</td>
          <td>${esc(ch.phenotype_mim || "")} ${ch.before && ch.before.label ? `${esc(ch.before.label)} → ` : ""}${esc((ch.after && (ch.after.label || ch.after.title)) || (ch.before && ch.before.label) || "")}</td></tr>`).join("")
      }</tbody></table></details>` : "";
    return `
      <div class="lookup-head">
        <h3>${g.url ? `<a href="${esc(g.url)}" target="_blank" rel="noopener">${esc(g.gene)} ${esc(g.gene_mim)}</a>` : `${esc(g.gene)} ${esc(g.gene_mim)}`}</h3>
        ${statusBadge} ${badge([g.source === "db" ? "平台資料庫" : "OMIM.xlsx", "info"])}
      </div>
      <p class="muted-text">${g.aliases.length ? `別名：${esc(g.aliases.join("、"))} · ` : ""}Inheritance：${esc(g.inheritance || "—")} · Done：${esc(g.done || "—")}${g.legacy_done && g.legacy_done !== g.done ? `（OMIM.xlsx 原 Done：${esc(g.legacy_done)}）` : ""} · 基因頁同步：${esc(g.synced_at ? fmtTime(g.synced_at) : "尚未")}</p>
      ${g.diseases.length ? `<table class="data-table lookup-table"><thead><tr><th>#</th><th>疾病</th><th>內容來源</th><th>報告</th></tr></thead><tbody>${diseases}</tbody></table>`
        : '<p class="muted-text">平台沒有這個基因的疾病資料。</p>'}
      ${steps}${reviews}${removed}${changes}`;
  }

  async function runLookup(q) {
    q = String(q || "").trim();
    if (!q) return;
    $("#lookup-q").value = q;
    lookupItems = [];
    hideLookupDropdown();
    const box = $("#lookup-result");
    box.hidden = false;
    box.innerHTML = "查詢中…";
    try {
      const data = await api(`/api/omim/gene?q=${encodeURIComponent(q)}`);
      if (data.gene) {
        box.innerHTML = renderGene(data.gene);
      } else if (data.matches.length) {
        box.innerHTML = `<p>${esc(q)} 連到多個基因，請選一個：</p><div class="lookup-choices">${data.matches.map(geneButton).join("")}</div>`;
      } else {
        box.innerHTML = `<p>找不到「${esc(q)}」。</p>${(data.suggestions || []).length
          ? `<p class="muted-text">是不是要找：</p><div class="lookup-choices">${data.suggestions.map(geneButton).join("")}</div>` : ""}`;
      }
    } catch (err) {
      box.textContent = `查詢失敗：${err.message}`;
    }
  }

  function list(title, items, fmt) {
    return items && items.length
      ? `<p><b>${esc(title)}</b></p><ul>${items.map(i => `<li>${fmt(i)}</li>`).join("")}</ul>` : "";
  }

  $("#btn-token-create").addEventListener("click", async () => {
    const created = await api("/api/omim/capture-tokens", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ label: $("#token-label").value }),
    });
    $("#token-value").textContent = created.token;
    $("#token-new").hidden = false;
    $("#token-label").value = "";
    loadTokens();
  });
  $("#btn-token-copy").addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText($("#token-value").textContent);
      $("#btn-token-copy").textContent = "已複製";
    } catch (_) {
      const range = document.createRange();
      range.selectNodeContents($("#token-value"));
      getSelection().removeAllRanges();
      getSelection().addRange(range);
    }
  });
  $("#token-table").addEventListener("click", async ev => {
    const id = ev.target.closest("[data-revoke]")?.dataset.revoke;
    if (!id || !confirm("停用這個 token？使用它的擴充功能會無法再送出擷取。")) return;
    await api(`/api/omim/capture-tokens/${id}`, { method: "DELETE" });
    loadTokens();
  });
  $("#capture-table").addEventListener("click", ev => {
    const id = ev.target.closest("[data-capture]")?.dataset.capture;
    if (id) showCapture(id).catch(err => { $("#capture-detail").textContent = String(err); });
  });
  $("#btn-refresh").addEventListener("click", () => { loadCaptures(); loadDbStatus(); loadQueue(); });
  document.querySelectorAll(".tab").forEach(t => t.addEventListener("click", () => setView(t.dataset.view)));
  $("#adoption-table").addEventListener("click", ev => {
    const id = ev.target.closest("[data-capture]")?.dataset.capture;
    if (id) showCapture(id).catch(err => { $("#capture-detail").textContent = String(err); });
  });
  $("#btn-adopt-all").addEventListener("click", async () => {
    if (!pendingIds.length) return;
    if (!confirm(`採用 ${pendingIds.length} 筆擷取？有「需確認」標記的會被跳過。`)) return;
    const btn = $("#btn-adopt-all");
    btn.disabled = true;
    btn.textContent = "採用中…";
    const box = $("#bulk-result");
    try {
      const res = await api("/api/omim/captures/adopt-pending", {
        method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ ids: pendingIds }),
      });
      const reviews = res.adopted.reduce((n, a) => n + (a.reviews || 0), 0);
      box.innerHTML = `已採用 ${res.adopted.length} 筆${reviews ? `，新增 ${reviews} 項待確認` : ""}。`
        + (res.skipped.length ? `<br>跳過 ${res.skipped.length} 筆：<ul>${res.skipped.map(sk =>
          `<li>${esc(sk.mim || sk.id)} ${esc(sk.title || "")}：${esc((sk.reasons || []).join("、"))}</li>`).join("")}</ul>` : "");
    } catch (err) {
      box.textContent = `全部採用失敗：${err.message}`;
    }
    box.hidden = false;
    $("#capture-detail").hidden = true;
    await Promise.all([loadCaptures(), loadDbStatus(), loadQueue()]);
  });
  $("#capture-detail").addEventListener("click", async ev => {
    const btn = ev.target.closest("[data-adopt]");
    if (!btn) return;
    btn.disabled = true;
    const [mim, name] = btn.dataset.adopt.split("/");
    try {
      const res = await api(`/api/omim/captures/${encodeURIComponent(mim)}/${encodeURIComponent(name)}/adopt`, { method: "POST" });
      const genes = (res.result.genes || []).map(g => g.gene).filter(Boolean).join("、");
      $("#adopt-result").textContent = `已採用${genes ? `：${genes}` : ""}${res.result.reviews ? `；新增 ${res.result.reviews} 項待確認` : ""}${res.result.kept_previous_synopsis ? "；這次沒擷取到 Synopsis，保留原本的" : ""}`;
      loadDbStatus();
      loadCaptures();
    } catch (err) {
      btn.disabled = false;
      $("#adopt-result").textContent = `採用失敗：${err.message}`;
    }
  });
  $("#queue-status").addEventListener("change", () => loadQueue());
  $("#lookup-form").addEventListener("submit", ev => {
    ev.preventDefault();
    if (lookupActive >= 0) { pickLookup(lookupActive); return; }
    hideLookupDropdown();
    runLookup($("#lookup-q").value);
  });
  $("#lookup-q").addEventListener("input", () => {
    clearTimeout(lookupTimer);
    lookupTimer = setTimeout(() => lookupSuggest().catch(() => {}), 150);
  });
  $("#lookup-q").addEventListener("keydown", ev => {
    const open = !$("#lookup-dropdown").classList.contains("hidden");
    if (ev.key === "ArrowDown" && lookupItems.length) {
      ev.preventDefault();
      lookupActive = Math.min(lookupActive + 1, lookupItems.length - 1);
      renderLookupDropdown();
    } else if (ev.key === "ArrowUp" && open) {
      ev.preventDefault();
      lookupActive = Math.max(lookupActive - 1, 0);
      renderLookupDropdown();
    } else if (ev.key === "Escape" && open) {
      ev.preventDefault();
      hideLookupDropdown();
    }
  });
  $("#lookup-q").addEventListener("focus", () => { if (lookupItems.length) renderLookupDropdown(); });
  $("#lookup-q").addEventListener("blur", () => setTimeout(hideLookupDropdown, 150));
  $("#lookup-dropdown").addEventListener("mousedown", ev => {
    const li = ev.target.closest("[data-index]");
    if (!li) return;
    ev.preventDefault();
    pickLookup(Number(li.dataset.index));
  });
  document.addEventListener("click", ev => {
    const target = ev.target.closest("[data-lookup]");
    if (!target) return;
    ev.preventDefault();
    runLookup(target.dataset.lookup);
    $("#lookup-card").scrollIntoView({ behavior: "smooth", block: "start" });
  });
  $("#lookup-result").addEventListener("click", async ev => {
    const btn = ev.target.closest("[data-review]");
    if (!btn) return;
    const action = btn.dataset.action;
    if (action === "remove" && !confirm("移除這個疾病？卡片與報告將不再顯示它（已勾選的個案仍保留勾選時的名稱）。")) return;
    btn.disabled = true;
    await api(`/api/omim/reviews/${btn.dataset.review}`, {
      method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ action }),
    });
    runLookup($("#lookup-q").value);
    loadDbStatus();
  });
  $("#queue-scope").addEventListener("change", () => loadQueue());
  $("#btn-queue-more").addEventListener("click", () => loadQueue(true));
  $("#queue-table").addEventListener("click", async ev => {
    const btn = ev.target.closest("[data-review]");
    if (!btn) return;
    const action = btn.dataset.action;
    if (action === "remove" && !confirm("移除這個疾病？卡片與報告將不再顯示它（已勾選的個案仍保留勾選時的名稱）。")) return;
    btn.disabled = true;
    await api(`/api/omim/reviews/${btn.dataset.review}`, {
      method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ action }),
    });
    const li = btn.closest("[data-review-item]");
    if (li) li.innerHTML = `<span class="muted-text">${action === "remove" ? "已移除" : "已保留"}</span>`;
    loadDbStatus();
  });

  (async () => {
    try {
      const me = await api("/api/auth/me");
      $("#whoami").textContent = me.username;
      await Promise.all([loadTokens(), loadCaptures(), loadDbStatus(), loadQueue()]);
    } catch (_) { /* login hint already shown */ }
  })();
})();
