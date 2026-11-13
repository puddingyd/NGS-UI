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
      .nextbar { padding: 8px 10px; border-top: 1px solid #d0d7de; background: #f6f8fa; font-size: 12px; }
      .nextbar .next-label { margin: 2px 0 6px; }
      .nextbar .next-actions { display: flex; gap: 6px; }
      .nextbar .btn { flex: 0 0 auto; padding: 4px 8px; font-size: 12px; }
      .mini .nextbar { display: none; }
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
      <div class="nextbar" id="nextbar"><span class="muted">讀取待補清單中…</span></div>
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

  // Text between the relationships table and the TEXT section, where an
  // expanded synopsis shows up.
  function synopsisSlice(text) {
    let start = text.search(/(Phenotype|Gene)-(Gene|Phenotype) Relationships/);
    // Entries without a relationships table: start at the entry header.
    if (start < 0) start = text.search(new RegExp(`\\n\\s*[#%*+^]?\\s*${MIM}\\s*\\n`));
    if (start < 0) return "";
    const textIdx = text.search(/\n\s*▼\s*TEXT\b/);
    return text.slice(start, textIdx > start ? textIdx : undefined);
  }

  // Synopsis is visible when the fold's "▲ Close" link is shown (only there
  // when expanded; works for old-style headings such as "Neuro" / "Misc")
  // or a current-style category line sits in that slice.
  function synopsisVisible(text) {
    const slice = synopsisSlice(text);
    if (!slice) return false;
    if (/\n\s*▲\s*Close\s*(\n|$)/.test(slice)) return true;
    return SYNOPSIS_CATEGORIES.some(cat => new RegExp(`\\n\\s*${cat.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}\\s*\\n`).test(slice));
  }

  // True when following this element's link would leave the current page
  // (e.g. the orange "Clinical Synopsis" label links to /clinicalSynopsis/MIM).
  function navigates(el) {
    const a = el && el.closest ? el.closest("a[href]") : null;
    if (!a) return false;
    const raw = (a.getAttribute("href") || "").trim();
    if (!raw || raw.startsWith("#") || /^javascript:/i.test(raw)) return false;
    try {
      const url = new URL(a.href, location.href);
      return url.origin !== location.origin || url.pathname !== location.pathname
        || url.search !== location.search;
    } catch (_) {
      return true;
    }
  }

  // Safety net while we click: cancel any synthetic click that would
  // navigate away (that would kill this script before the capture is sent).
  // Real user clicks (isTrusted) are never touched.
  function blockSyntheticNavigation(ev) {
    if (ev.isTrusted) return;
    const target = ev.target && ev.target.closest ? ev.target.closest("a[href]") : null;
    if (target && navigates(target)) {
      ev.preventDefault();
      ev.stopImmediatePropagation();
    }
  }

  const TRIANGLE_RE = /^[►▼▶▸▾▴◀◂▲+\-−]$/;
  const TOGGLE_HINT_RE = /triangle|toggle|caret|collapse|fold|expand|chevron|arrow/i;

  // The same-page expander next to the "Clinical Synopsis" label: the small
  // triangle the reviewer normally clicks. Never the label link itself.
  function synopsisToggles() {
    const found = new Set();
    const add = el => {
      if (!el || host.contains(el) || navigates(el)) return;
      // Table-of-contents links only scroll the page.
      if (el.closest("[id*='toc' i], [class*='toc' i]")) return;
      found.add(el);
    };
    // 1. Explicit collapse controls that point at a synopsis fold.
    document.querySelectorAll([
      "[id*='ClinicalSynopsisToggle' i]",
      "[id*='ClinicalSynopsisTriangle' i]",
      "[data-toggle][href*='linicalSynopsis' i]",
      "[data-bs-toggle][href*='linicalSynopsis' i]",
      "[data-target*='linicalSynopsis' i]",
      "[data-bs-target*='linicalSynopsis' i]",
      "[aria-controls*='linicalSynopsis' i]",
      "a[href^='#'][href*='linicalSynopsis' i]",
    ].join(",")).forEach(add);
    // 2. Triangle / toggle-looking elements around a "Clinical Synopsis" label.
    const labels = [];
    document.querySelectorAll("a, span, div, h2, h3, h4, h5, strong, b").forEach(el => {
      if (host.contains(el) || el.children.length > 3) return;
      const label = (el.textContent || "").replace(/\s+/g, " ").trim();
      if (/^[^A-Za-z0-9]{0,4}\s*Clinical Synopsis\s*$/i.test(label)) labels.push(el);
    });
    for (const label of labels) {
      if (label.closest("[id*='toc' i], [class*='toc' i]")) continue;
      let scope = label;
      for (let depth = 0; depth < 3 && scope.parentElement; depth += 1) {
        scope = scope.parentElement;
        scope.querySelectorAll("*").forEach(el => {
          if (el === label || el.contains(label)) return;
          const text = (el.textContent || "").trim();
          const hint = `${el.id || ""} ${typeof el.className === "string" ? el.className : ""}`;
          if (TRIANGLE_RE.test(text) || (TOGGLE_HINT_RE.test(hint) && text.length <= 3)) add(el);
        });
        if (found.size) break;
      }
      // A label that is itself a same-page toggle (e.g. "► Clinical Synopsis").
      if (!navigates(label)) add(label);
    }
    const rank = el => (el.matches("[data-toggle], [data-bs-toggle], [aria-controls], a[href^='#']") ? 0 : 1)
      + (el.getClientRects().length ? 0 : 2);
    return [...found].sort((x, y) => rank(x) - rank(y)).slice(0, 6);
  }

  async function expandSynopsis() {
    if (synopsisVisible(pageText())) return "already_visible";
    document.addEventListener("click", blockSyntheticNavigation, true);
    try {
      for (const el of synopsisToggles()) {
        const target = el.closest("a, button, [role='button'], [onclick], [data-toggle], [data-bs-toggle]") || el;
        if (navigates(target)) continue;
        const before = synopsisSlice(pageText()).length;
        target.click();
        for (let i = 0; i < 12; i += 1) {
          await sleep(i < 4 ? 100 : 250);
          const text = pageText();
          if (synopsisVisible(text)) return "expanded";
          // Something opened even if we do not recognise its headings: stop
          // here instead of clicking another toggle (that could close it).
          if (synopsisSlice(text).length > before + 30) return "expanded";
        }
      }
    } finally {
      document.removeEventListener("click", blockSyntheticNavigation, true);
    }
    return "";
  }

  // The synopsis is often already in the page but collapsed (display:none).
  // Copy the smallest hidden block holding several synopsis categories into
  // an off-screen box to read it with its line structure.
  function hiddenSynopsisText() {
    let best = null;
    document.querySelectorAll("div, section, td, table").forEach(el => {
      if (host.contains(el) || el.getClientRects().length) return;
      const text = el.textContent || "";
      const hits = SYNOPSIS_CATEGORIES.filter(cat => text.includes(cat)).length;
      if (hits < 2) return;
      if (!best || text.length < (best.textContent || "").length) best = el;
    });
    if (!best) return "";
    const box = document.createElement("div");
    box.style.cssText = "position:absolute;left:-100000px;top:0;width:900px;";
    const clone = best.cloneNode(true);
    clone.style.display = "block";
    clone.querySelectorAll("*").forEach(n => {
      if (n.style && n.style.display === "none") n.style.display = "";
    });
    clone.classList.remove("collapse");
    box.appendChild(clone);
    document.body.appendChild(box);
    try {
      return box.innerText || "";
    } finally {
      box.remove();
    }
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
      ${result.outboxSent ? `<p class="muted">另外補送了 ${result.outboxSent} 筆先前未送出的擷取。</p>` : ""}
      <div id="adopt-box" class="row">${renderAdoptBox(summary)}</div>`;
    const adoptBtn = root.querySelector("#adopt");
    if (adoptBtn) adoptBtn.addEventListener("click", () => runAdopt(capture.id, adoptBtn, false));
    const adoptNext = root.querySelector("#adopt-next");
    if (adoptNext) adoptNext.addEventListener("click", () => runAdopt(capture.id, adoptNext, true));
  }

  function renderAdoptBox(summary) {
    if (!summary.db_enabled) {
      return `<p class="muted">平台尚未啟用 OMIM 資料庫，暫時只保存擷取內容，無法採用。</p>`;
    }
    const preview = summary.adopt_preview || {};
    const isGene = summary.entry_type === "gene";
    const lines = (preview.genes || []).map(g => {
      const current = (g.current_labels || []).join("、");
      const change = g.fills_mim_for ? `補上 MIM，原為「${esc(g.fills_mim_for)}」`
        : (!current ? "新增" : (current === g.label ? "沿用" : `原為「${esc(current)}」`));
      return `<li>${esc(g.gene)}：${esc(g.label)}（${change}）</li>`;
    }).join("");
    const unlisted = preview.unlisted || [];
    if (summary.entry_type === "moved_removed") {
      const moved = summary.moved_to || [];
      return `<p class="muted">${moved.length
        ? `OMIM 已把這個 entry 搬到 ${esc(moved.join("、"))}。平台已記下，待補清單會改到新的 MIM；採用只是把這筆擷取標成已處理。`
        : "OMIM 已移除這個 entry。採用後，如果它是平台上的基因，會列入「待確認」讓你決定。"}</p>
        <div class="next-actions" style="display:flex;gap:6px"><button class="btn" id="adopt">採用到平台</button></div>`;
    }
    let note = preview.omim_lists_no_disease
      ? (unlisted.length
        ? `OMIM 沒有列出這個基因的疾病。採用後記錄為已核對；平台現有的 ${unlisted.length} 個疾病（${esc(unlisted.join("、"))}）會列入待確認，不會刪除。`
        : "OMIM 沒有列出這個基因的疾病，平台也沒有；採用後記錄為已核對。")
      : isGene ? "採用後會依 OMIM 基因頁同步這個基因的疾病清單（OMIM 已不列出的舊疾病只會列入待確認，不會刪除）。"
      : (preview.same_content ? "平台已採用過相同內容。"
        : (preview.text_exists ? "會取代平台上這個疾病目前的說明與 Clinical Synopsis。"
          : "會新增這個疾病的說明與 Clinical Synopsis。"));
    if (!isGene && summary.entry_type !== "gene_phenotype" && summary.relationship_section === "absent") {
      note += "<br>OMIM 頁面沒有 Gene-Phenotype Relationships（例如染色體異常），採用後只保存說明，不會連到任何基因，也會從待補清單移除。";
    }
    if (!isGene && summary.synopsis_status === "not_loaded") {
      note += "<br><b>注意：Clinical Synopsis 沒有擷取到，建議展開後重新 Capture 再採用。</b>";
    }
    const blocked = (summary.blockers || []).length > 0;
    return `
      ${lines ? `<b>採用後</b><ul>${lines}</ul>` : ""}
      <p class="muted">${note}</p>
      <div class="next-actions" style="display:flex;gap:6px">
        <button class="btn" id="adopt" ${preview.same_content ? "disabled" : ""}>採用到平台</button>
        ${blocked || preview.same_content ? "" : `<button class="btn" id="adopt-next">採用並前往下一筆</button>`}
      </div>`;
  }

  async function runAdopt(captureId, button, goNext) {
    const label = button.textContent;
    root.querySelectorAll("#adopt, #adopt-next").forEach(b => { b.disabled = true; });
    button.textContent = "採用中…";
    const result = await send({ type: "adopt", captureId });
    const box = root.querySelector("#adopt-box");
    if (!result.ok) {
      root.querySelectorAll("#adopt, #adopt-next").forEach(b => { b.disabled = false; });
      button.textContent = label;
      box.insertAdjacentHTML("beforeend", `<p class="badge err">${esc(result.error)}</p>`);
      return;
    }
    const r = result.data.result || {};
    const genes = (r.genes || []).map(g => esc(g.gene)).filter(Boolean).join("、");
    box.innerHTML = `
      <span class="badge ok">已採用到平台</span>
      ${genes ? `<p class="muted">已更新：${genes}${r.text ? `（說明${r.text === "added" ? "新增" : "已取代"}）` : ""}</p>` : ""}

      ${r.kept_previous_synopsis ? `<p class="muted">這次沒擷取到 Synopsis，已保留平台原本的內容。</p>` : ""}
      <p class="muted">平台上的個案重新載入後會顯示新內容。</p>`;
    const geneReviews = r.gene_reviews || [];
    if (geneReviews.length) {
      box.insertAdjacentHTML("beforeend", renderGeneReviews(geneReviews)
        + (goNext ? `<p class="muted">處理完上面的項目後，按下方「前往」繼續。</p>` : ""));
      box.addEventListener("click", handleReviewClick);
    }
    if (goNext && !geneReviews.length) {
      await goToNext(null);
    } else {
      loadNext();
    }
  }

  // ---------------------------------------------------------- work queue
  let nextData = null;
  const KIND_TEXT = {
    gene_page: "同步這個基因的疾病清單", gene_updated: "OMIM 之後更新了基因頁，重新同步",
    gene_moved: "OMIM 已把這個基因搬到新的 MIM，同步新的基因頁",
    missing: "缺說明", updated: "OMIM 之後有更新",
    new_disease: "平台沒有的疾病（OMIM 有更新），採用後會加到所屬基因",
    new_gene: "平台沒有的基因（OMIM 有更新），採用後會建立這個基因",
  };

  async function loadNext() {
    const bar = root.querySelector("#nextbar");
    const result = await send({ type: "next", current: MIM });
    if (!result.ok) {
      bar.innerHTML = `<span class="muted">待補清單：${esc(result.error)}</span>`;
      nextData = null;
      return null;
    }
    nextData = result.data;
    if (!nextData.step) {
      bar.innerHTML = `<span class="muted">待補清單已完成 🎉</span>`;
      return nextData;
    }
    const g = nextData.gene;
    const step = nextData.step;
    const kinds = (step.kinds || []).map(k => KIND_TEXT[k] || k).join("、")
      + (step.update_date ? `（OMIM ${step.update_date}）` : "");
    const pseudo = String(g.gene_mim || "").startsWith("new");
    const done = g.total_steps - g.remaining_steps;
    bar.innerHTML = `
      <b>${nextData.continuing ? "同一個基因的下一步" : (nextData.related ? "下一個基因（共用剛才的疾病）" : "下一個基因")}</b>
      <span class="muted">（還有 ${esc(nextData.remaining_genes)} 個基因）</span>
      <div class="next-label"><b>${esc(g.gene)}</b>${g.panel ? " · 疾病相關" : ""} · 第 ${done + 1}／${g.total_steps} 步${g.reviews ? ` · 待確認 ${esc(g.reviews)} 項` : ""}
        <br>${step.type === "gene" ? `基因頁 ${esc(step.mim)}` : `${esc(step.mim)} ${esc((step.labels || [])[0] || "")}`}
        <br><span class="muted">${esc(kinds)}${g.reported ? ` · 報告用過 ${esc(g.reported)} 次` : ""}</span></div>
      <div class="next-actions">
        <button class="btn" id="go-next">前往</button>
        <button class="btn" id="capture-next" title="把目前這一頁 Capture 送到平台（不採用），再前往下一步；之後到平台的 OMIM 擷取頁「全部採用」">Capture 並前往</button>
        <button class="btn secondary" id="skip-current" title="這一頁先不處理，之後不再推薦">略過這頁</button>
        ${pseudo ? "" : `<button class="btn secondary" id="skip-gene" title="${esc(g.gene)} 整個基因先不處理">略過 ${esc(g.gene)}</button>`}
      </div>`;
    root.querySelector("#go-next").addEventListener("click", () => goToNext(null));
    root.querySelector("#capture-next").addEventListener("click", () => captureAndNext());
    root.querySelector("#skip-current").addEventListener("click", () => goToNext("page"));
    if (!pseudo) root.querySelector("#skip-gene").addEventListener("click", () => goToNext("gene", g.gene_mim));
    return nextData;
  }

  async function goToNext(skipMode, geneMim) {
    if (skipMode === "page") await send({ type: "skip", mim: MIM });
    if (skipMode === "gene") await send({ type: "skip", mim: `g${geneMim}` });
    // Always ask again: an adoption may have just changed this gene's steps.
    const data = await loadNext();
    if (data && data.step) {
      await send({ type: "set-gene", gene: data.gene.gene_mim });
      location.href = `${location.origin}/entry/${data.step.mim}`;
      return true;
    }
    return false;
  }

  function renderGeneReviews(reviews) {
    if (!reviews || !reviews.length) return "";
    return `<div class="row"><span class="badge warn">這個基因有 ${reviews.length} 項對不上，請確認</span><ul id="gene-reviews">${
      reviews.map(r => {
        const d = r.detail || {};
        const label = r.current_label || d.label || d.disease_cell || d.omim_disease || "";
        return `<li data-review="${esc(r.id)}">${esc(r.kind_label)}：<b>${esc(label)}</b>
          ${d.suggestion ? `<br><span class="muted">可能是：${esc(d.suggestion)}</span>` : ""}
          <br><button class="btn secondary" data-action="keep">保留</button>
          ${r.disease_id && r.disease_active ? `<button class="btn secondary" data-action="remove">移除疾病</button>` : ""}</li>`;
      }).join("")
    }</ul></div>`;
  }

  async function handleReviewClick(ev) {
    const btn = ev.target.closest("[data-action]");
    const li = ev.target.closest("[data-review]");
    if (!btn || !li) return;
    btn.disabled = true;
    const result = await send({ type: "review", id: li.dataset.review, action: btn.dataset.action });
    if (result.ok) {
      li.innerHTML = `<span class="muted">${btn.dataset.action === "remove" ? "已移除" : "已保留"}</span>`;
      if (!root.querySelector("#gene-reviews [data-action]")) loadNext();
    } else {
      btn.disabled = false;
      li.insertAdjacentHTML("beforeend", `<br><span class="badge err">${esc(result.error)}</span>`);
    }
  }

  // ------------------------------------------------------------ actions
  async function refreshStatus() {
    const localEdit = editDate(pageText());
    renderStatus(await send({ type: "status", mim: MIM, editDate: localEdit }), localEdit);
  }

  // Capture this page (without adopting) and move on; stops on this page
  // when the capture failed or needs a look (e.g. Synopsis not captured), so
  // it can be fixed right here. Adopt everything later on /omim/.
  async function captureAndNext() {
    const buttons = root.querySelectorAll("#go-next, #capture-next, #skip-current, #skip-gene");
    buttons.forEach(b => { b.disabled = true; });
    const button = root.querySelector("#capture-next");
    button.textContent = "擷取中…";
    const result = await runCapture();
    const blockers = result && result.ok ? (result.data.summary.blockers || []) : null;
    if (blockers && !blockers.length) {
      button.textContent = "前往中…";
      const moved = await goToNext(null);
      if (moved) return;
      button.textContent = "Capture 並前往";
      return;
    }
    buttons.forEach(b => { b.disabled = false; });
    button.textContent = "Capture 並前往";
    const msg = blockers ? `這一頁需要確認（${blockers.join("、")}），沒有自動前往；可重新 Capture，或按「前往」先跳過。`
      : "Capture 沒有成功送到平台，沒有自動前往。";
    $("#body").insertAdjacentHTML("afterbegin", `<p class="badge warn">${esc(msg)}</p>`);
  }

  async function runCapture() {
    const button = $("#capture");
    button.disabled = true;
    button.textContent = "擷取中…";
    let result = null;
    try {
      const prefix = entryPrefix(pageText());
      let synopsisMethod = "";
      let synopsisText = "";
      if (prefix !== "*") {
        synopsisMethod = await expandSynopsis();
        if (!synopsisMethod) {
          synopsisText = hiddenSynopsisText();
          synopsisMethod = synopsisText ? "hidden_fold" : "";
        }
        if (!synopsisMethod) {
          // Background fetch only — the page itself never navigates.
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
      result = await send({ type: "capture", payload });
      renderCapture(result);
    } catch (err) {
      result = { ok: false, error: String(err) };
      renderCapture(result);
    } finally {
      button.disabled = false;
      button.textContent = "Capture 送到平台";
    }
    return result;
  }

  refreshStatus();
  loadNext();
})();
