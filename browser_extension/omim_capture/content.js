// NGS OMIM Capture — content script for omim.org/entry/* pages.
//
// Shows a small panel with the platform's status for this MIM. The page is
// read only when the reviewer presses "Capture": we expand the Clinical
// Synopsis if needed, take the rendered page text (the same thing a manual
// Ctrl+A / Ctrl+C gives), and hand it to the background worker, which
// sends it to the NGS platform. No other OMIM pages are opened or crawled,
// except one same-origin /clinicalSynopsis/{MIM} request as a fallback when
// the synopsis cannot be expanded in place.

(() => {
  if (window.__ngsOmimCapture) return;
  window.__ngsOmimCapture = true;

  const MIM = (location.pathname.match(/\/entry\/(\d{6})/) || [])[1] || "";
  if (!MIM) return;

  const SYNOPSIS_CATEGORIES = [
    "INHERITANCE", "GROWTH", "HEAD & NECK", "CARDIOVASCULAR", "RESPIRATORY",
    "CHEST", "ABDOMEN", "GENITOURINARY", "SKELETAL", "SKIN, NAILS, & HAIR",
    "MUSCLE, SOFT TISSUES", "NEUROLOGIC", "VOICE", "METABOLIC FEATURES",
    "ENDOCRINE FEATURES", "HEMATOLOGY", "IMMUNOLOGY", "NEOPLASIA",
    "PRENATAL MANIFESTATIONS", "LABORATORY ABNORMALITIES", "MISCELLANEOUS",
    "MOLECULAR BASIS",
  ];
  const STATUS_TEXT = {
    not_in_workbook: ["平台工作簿沒有這個 MIM", "warn"],
    missing_detail: ["平台只有疾病名稱，缺說明（*）", "warn"],
    outdated: ["OMIM 編輯日期比平台 Done 新，需要更新", "warn"],
    no_done_date: ["平台有資料，但沒有 Done 日期", "info"],
    up_to_date: ["平台資料已是最新", "ok"],
    gene_in_workbook: ["基因在平台工作簿中", "info"],
  };
  const SYNOPSIS_STATUS = {
    present: "已取得",
    not_loaded: "未展開／未載入",
    absent: "此 entry 沒有 Clinical Synopsis",
    not_applicable: "基因 entry（不適用）",
  };

  const esc = value => String(value ?? "").replace(/[&<>"']/g, ch => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[ch]
  ));
  const sleep = ms => new Promise(resolve => setTimeout(resolve, ms));
  const send = message => new Promise(resolve => {
    try {
      chrome.runtime.sendMessage(message, response => {
        if (chrome.runtime.lastError) {
          resolve({ ok: false, error: chrome.runtime.lastError.message });
        } else {
          resolve(response || { ok: false, error: "no response" });
        }
      });
    } catch (err) {
      resolve({ ok: false, error: String(err) });
    }
  });

  // ---------------------------------------------------------------- panel
  const host = document.createElement("div");
  host.id = "ngs-omim-capture-host";
  host.style.cssText = "position:fixed;right:16px;bottom:16px;z-index:2147483647;";
  const root = host.attachShadow({ mode: "open" });
  root.innerHTML = `
    <style>
      :host { all: initial; }
      .panel { font: 13px/1.45 -apple-system, "Segoe UI", "Noto Sans TC", sans-serif;
        color: #1f2328; background: #fff; border: 1px solid #d0d7de; border-radius: 8px;
        box-shadow: 0 6px 24px rgba(31,35,40,.18); width: 360px; max-height: 70vh;
        display: flex; flex-direction: column; overflow: hidden; }
      .head { display: flex; align-items: center; gap: 8px; padding: 8px 10px;
        background: #24292f; color: #fff; }
      .head b { flex: 1; font-size: 13px; }
      .head button { background: none; border: 0; color: #fff; cursor: pointer; font-size: 14px; }
      .body { padding: 10px; overflow: auto; }
      .mini .body, .mini .foot { display: none; }
      .mini.panel { width: auto; }
      .foot { display: flex; gap: 6px; padding: 8px 10px; border-top: 1px solid #d0d7de; }
      .btn { flex: 1; border: 1px solid #1f6feb; background: #1f6feb; color: #fff;
        border-radius: 6px; padding: 6px 8px; cursor: pointer; font-size: 13px; }
      .btn.secondary { flex: 0 0 auto; background: #f6f8fa; color: #1f2328; border-color: #d0d7de; }
      .btn:disabled { opacity: .55; cursor: default; }
      .badge { display: inline-block; border-radius: 10px; padding: 1px 8px; font-size: 12px; }
      .ok { background: #dafbe1; color: #116329; }
      .warn { background: #fff8c5; color: #7d4e00; }
      .info { background: #ddf4ff; color: #0550ae; }
      .err { background: #ffebe9; color: #a40e26; }
      .muted { color: #59636e; font-size: 12px; }
      .row { margin: 6px 0; }
      ul { margin: 4px 0 0 18px; padding: 0; }
      li { margin: 2px 0; }
      details { margin-top: 6px; }
      pre { white-space: pre-wrap; font-size: 11px; background: #f6f8fa; padding: 6px;
        border-radius: 4px; max-height: 160px; overflow: auto; }
    </style>
    <div class="panel" id="panel">
      <div class="head"><b>NGS OMIM Capture · ${esc(MIM)}</b>
        <button id="min" title="縮小／展開">—</button></div>
      <div class="body" id="body"><span class="muted">讀取平台狀態中…</span></div>
      <div class="foot">
        <button class="btn" id="capture">Capture 送到平台</button>
        <button class="btn secondary" id="options" title="設定平台網址與 token">設定</button>
      </div>
    </div>`;
  document.documentElement.appendChild(host);
  const $ = sel => root.querySelector(sel);

  chrome.storage.local.get(["panelMinimized"], s => {
    if (s.panelMinimized) $("#panel").classList.add("mini");
  });
  $("#min").addEventListener("click", () => {
    const mini = $("#panel").classList.toggle("mini");
    chrome.storage.local.set({ panelMinimized: mini });
  });
  $("#options").addEventListener("click", () => send({ type: "open-options" }));
  $("#capture").addEventListener("click", () => runCapture());

  // ------------------------------------------------------- page reading
  function pageText() {
    const display = host.style.display;
    host.style.display = "none";
    try {
      return document.body ? document.body.innerText : "";
    } finally {
      host.style.display = display;
    }
  }

  function entryPrefix(text) {
    const match = text.match(new RegExp(`^\\s*([#%*+^])?\\s*${MIM}\\s*$`, "m"));
    return match ? (match[1] || "") : "";
  }

  function editDate(text) {
    const idx = text.indexOf("Edit History:");
    if (idx < 0) return "";
    const block = text.slice(idx, idx + 2000).split(/\n\s*\n/)[0];
    const dates = [...block.matchAll(/(\d{1,2})\/(\d{1,2})\/(\d{4})/g)]
      .map(m => `${m[3]}-${m[1].padStart(2, "0")}-${m[2].padStart(2, "0")}`);
    return dates.sort().pop() || "";
  }

  // Synopsis is visible when a category line sits between the relationships
  // table and the TEXT section.
  function synopsisVisible(text) {
    const start = text.search(/(Phenotype|Gene)-(Gene|Phenotype) Relationships/);
    const textIdx = text.search(/\n\s*▼\s*TEXT\b/);
    if (start < 0) return false;
    const slice = text.slice(start, textIdx > start ? textIdx : undefined);
    return SYNOPSIS_CATEGORIES.some(cat => new RegExp(`\\n\\s*${cat.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}\\s*\\n`).test(slice));
  }

  function synopsisToggles() {
    const out = [];
    document.querySelectorAll("a, button, span, div, h3, h4, h5, strong").forEach(el => {
      if (host.contains(el) || el.children.length > 2) return;
      const label = (el.textContent || "").replace(/\s+/g, " ").trim();
      if (!/^[^A-Za-z0-9]{0,4}\s*Clinical Synopsis\s*$/i.test(label)) return;
      if (!el.getClientRects().length) return;
      out.push(el);
    });
    // Prefer clickable elements outside the table of contents (TOC links
    // only scroll the page).
    const rank = el => (el.closest("[id*='toc' i], [class*='toc' i]") ? 2 : 0)
      + (el.matches("a, button, [onclick], [role='button'], [data-toggle]") ? 0 : 1);
    return out.sort((a, b) => rank(a) - rank(b));
  }

  async function expandSynopsis() {
    if (synopsisVisible(pageText())) return "already_visible";
    for (const el of synopsisToggles()) {
      el.click();
      for (let i = 0; i < 20; i += 1) {
        await sleep(250);
        if (synopsisVisible(pageText())) return "expanded";
      }
    }
    return "";
  }

  async function synopsisPageText() {
    try {
      const response = await fetch(`/clinicalSynopsis/${MIM}`, { credentials: "include" });
      if (!response.ok) return "";
      const doc = new DOMParser().parseFromString(await response.text(), "text/html");
      doc.querySelectorAll("script, style, link, img, iframe, noscript").forEach(n => n.remove());
      const box = document.createElement("div");
      box.style.cssText = "position:absolute;left:-100000px;top:0;width:900px;white-space:normal;";
      box.append(...Array.from(doc.body ? doc.body.childNodes : []));
      document.body.appendChild(box);
      try {
        return box.innerText || "";
      } finally {
        box.remove();
      }
    } catch (_) {
      return "";
    }
  }

  function relationshipRows() {
    for (const table of document.querySelectorAll("table")) {
      const head = (table.innerText || "").slice(0, 400);
      if (!/MIM\s*number/i.test(head) || !/Phenotype/i.test(head) || !/Inheritance/i.test(head)) continue;
      const rows = [];
      table.querySelectorAll("tr").forEach(tr => {
        if (tr.querySelector("th")) return;
        const cells = [...tr.querySelectorAll("td")].map(td => (td.innerText || "").trim());
        if (cells.some(Boolean)) rows.push(cells);
      });
      if (rows.length) return rows;
    }
    return [];
  }

  // ------------------------------------------------------------- render
  function badge(status) {
    const [label, cls] = STATUS_TEXT[status] || [status || "未知", "info"];
    return `<span class="badge ${cls}">${esc(label)}</span>`;
  }

  function renderLinks(comparison) {
    const links = (comparison && comparison.phenotype_links) || [];
    if (!links.length) return "";
    return `<div class="row"><b>平台工作簿</b><ul>${links.map(l => `
      <li>${esc(l.gene)}：Disease${esc((l.slots || []).join(",") || "—")}
        ${l.has_detail ? "有說明" : "只有名稱"} · Done ${esc(l.done || "—")} ${badge(l.status)}</li>`).join("")}
    </ul></div>`;
  }

  function renderGene(gene) {
    if (!gene) return "";
    if (!gene.in_workbook) return `<div class="row">${badge("not_in_workbook")}（基因）</div>`;
    const list = (title, items, fmt) => items && items.length
      ? `<div class="row"><b>${esc(title)}</b><ul>${items.map(i => `<li>${fmt(i)}</li>`).join("")}</ul></div>` : "";
    return `
      <div class="row"><b>基因 ${esc(gene.gene)}</b> · Done ${esc(gene.done || "—")}</div>
      ${list("OMIM 有、平台沒有", gene.added_in_omim, r => `${esc(r.phenotype)} (${esc(r.phenotype_mim)})(${esc(r.inheritance)})`)}
      ${list("平台有、OMIM 已不列出", gene.missing_from_omim, esc)}
      ${list("名稱不同（同 MIM）", gene.renamed, r => `平台：${esc(r.workbook)}<br>OMIM：${esc(r.omim)}`)}
      ${list("平台缺 MIM 的疾病", gene.workbook_lines_without_mim, esc)}`;
  }

  function renderStatus(result, localEditDate) {
    if (!result.ok) {
      $("#body").innerHTML = `<div class="badge err">無法連到平台</div><p class="muted">${esc(result.error)}</p>`;
      return;
    }
    const { comparison, latest_capture: latest } = result.data;
    $("#body").innerHTML = `
      <div class="row">${badge(comparison.status)}</div>
      <div class="row muted">此頁 OMIM 最後編輯：${esc(localEditDate || "—")}
        ${latest ? `<br>平台最後擷取：${esc(latest.captured_at)}（${esc(latest.user)}）` : ""}</div>
      ${renderLinks(comparison)}
      ${comparison.gene ? renderGene(comparison.gene) : ""}`;
  }

  function renderCapture(result) {
    if (!result.ok) {
      $("#body").innerHTML = `
        <div class="badge err">${result.queued ? "暫存在瀏覽器，稍後自動重送" : "擷取失敗"}</div>
        <p class="muted">${esc(result.error)}</p>
        ${result.queued ? `<p class="muted">待送筆數：${esc(result.outbox)}</p>` : ""}`;
      return;
    }
    const { capture, summary } = result.data;
    const own = (summary.relationships || []);
    $("#body").innerHTML = `
      <div class="row"><span class="badge ok">${capture.duplicate ? "已送出（內容與上次相同）" : "已保存到平台"}</span>
        ${badge(summary.comparison.status)}</div>
      <div class="row"><b>${esc(summary.prefix)}${esc(summary.mim)} ${esc(summary.title)}</b></div>
      ${own.length ? `<div class="row"><b>關係表</b><ul>${own.map(r => `
        <li>${esc(r.phenotype)} (${esc(r.phenotype_mim || "無 MIM")})(${esc(r.inheritance)}) · ${esc(r.gene)}</li>`).join("")}</ul></div>` : ""}
      <div class="row">Description：${summary.description_chars ? `${summary.description_chars} 字（${esc(summary.description_source)}）` : "—"}
        <br>Clinical Synopsis：${esc(SYNOPSIS_STATUS[summary.synopsis_status] || summary.synopsis_status)}
        ${summary.synopsis_lines ? `（${summary.synopsis_lines} 行${summary.synopsis_source === "clinical_synopsis_page" ? "，取自 Synopsis 頁" : ""}）` : ""}
        <br>OMIM 最後編輯：${esc(summary.omim_edit_date || "—")}</div>
      ${(summary.warnings || []).length ? `<div class="row"><span class="badge warn">注意</span><ul>${summary.warnings.map(w => `<li>${esc(w)}</li>`).join("")}</ul></div>` : ""}
      ${renderLinks(summary.comparison)}
      ${summary.comparison.gene ? renderGene(summary.comparison.gene) : ""}
      ${(summary.disease_cells || []).map(c => `
        <details><summary>預覽 ${esc(c.gene)} 的 Disease 欄內容</summary><pre>${esc(c.cell)}</pre></details>`).join("")}
      ${result.outboxSent ? `<p class="muted">另外補送了 ${result.outboxSent} 筆先前未送出的擷取。</p>` : ""}`;
  }

  // ------------------------------------------------------------ actions
  async function refreshStatus() {
    const localEdit = editDate(pageText());
    renderStatus(await send({ type: "status", mim: MIM, editDate: localEdit }), localEdit);
  }

  async function runCapture() {
    const button = $("#capture");
    button.disabled = true;
    button.textContent = "擷取中…";
    try {
      const prefix = entryPrefix(pageText());
      let synopsisMethod = "";
      let synopsisText = "";
      if (prefix !== "*") {
        synopsisMethod = await expandSynopsis();
        if (!synopsisMethod) {
          synopsisText = await synopsisPageText();
          synopsisMethod = synopsisText ? "synopsis_page" : "not_found";
        }
      }
      const payload = {
        mim: MIM,
        url: location.href,
        title: document.title,
        text: pageText(),
        extracted: {
          relationship_rows: relationshipRows(),
          synopsis_text: synopsisText,
          synopsis_method: synopsisMethod,
        },
        extension_version: chrome.runtime.getManifest().version,
        captured_at_client: new Date().toISOString(),
      };
      renderCapture(await send({ type: "capture", payload }));
    } catch (err) {
      renderCapture({ ok: false, error: String(err) });
    } finally {
      button.disabled = false;
      button.textContent = "Capture 送到平台";
    }
  }

  refreshStatus();
})();
