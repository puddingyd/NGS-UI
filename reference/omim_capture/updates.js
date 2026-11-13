// NGS OMIM Capture — content script for the OMIM Update List
// (overview omim.org/statistics/update, months /statistics/updates/{year}/{month}).
//
// On a month page: "擷取這個月份" reads the listed entries (date, New /
// Updated Entries or Clinical Synopses, MIM) and sends them to the platform.
// On the overview page: "擷取缺少的月份" reads the months the platform still
// needs, one same-origin request every few seconds, stopping on any error.
// Nothing is read until a button is pressed.

(() => {
  if (window.__ngsOmimUpdates) return;
  window.__ngsOmimUpdates = true;

  const PAGE = location.pathname.match(/\/statistics\/updates\/(\d{4})\/(\d{1,2})\/?$/);
  const DELAY_MS = 3000;
  const DATE_RE = /^(January|February|March|April|May|June|July|August|September|October|November|December)\s+\d{1,2}(st|nd|rd|th)?,?\s+\d{4}$/i;
  const SECTION_RE = /^(New|Updated)\s+(Entries|Clinical Synopses)\s*:?$/i;
  const HREF_RE = /\/(entry|clinicalSynopsis)\/(\d{6})(?:[/?#]|$)/;

  const esc = value => String(value ?? "").replace(/[&<>"']/g, ch => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[ch]
  ));
  const sleep = ms => new Promise(resolve => setTimeout(resolve, ms));
  const send = message => new Promise(resolve => {
    try {
      chrome.runtime.sendMessage(message, response => {
        resolve(response || { ok: false, error: chrome.runtime.lastError?.message || "no response" });
      });
    } catch (err) {
      resolve({ ok: false, error: String(err) });
    }
  });
  const ym = (y, m) => `${y}-${String(m).padStart(2, "0")}`;

  const host = document.createElement("div");
  host.id = "ngs-omim-updates-host";
  host.style.cssText = "position:fixed;right:16px;bottom:16px;z-index:2147483647;";
  const root = host.attachShadow({ mode: "open" });
  root.innerHTML = `
    <style>
      :host { all: initial; }
      .panel { font: 13px/1.45 -apple-system, "Segoe UI", "Noto Sans TC", sans-serif;
        color: #1f2328; background: #fff; border: 1px solid #d0d7de; border-radius: 8px;
        box-shadow: 0 6px 24px rgba(31,35,40,.18); width: 340px; max-height: 70vh;
        display: flex; flex-direction: column; overflow: hidden; }
      .head { display: flex; align-items: center; gap: 8px; padding: 8px 10px; background: #24292f; color: #fff; }
      .head b { flex: 1; }
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
      .err { background: #ffebe9; color: #a40e26; }
      .muted { color: #59636e; font-size: 12px; }
      .row { margin: 6px 0; }
      progress { width: 100%; }
    </style>
    <div class="panel" id="panel">
      <div class="head"><b>NGS OMIM 更新清單</b><button id="min" title="縮小／展開">—</button></div>
      <div class="body" id="body"><span class="muted">讀取平台狀態中…</span></div>
      <div class="foot" id="foot"></div>
    </div>`;
  document.documentElement.appendChild(host);
  const $ = sel => root.querySelector(sel);
  $("#min").addEventListener("click", () => $("#panel").classList.toggle("mini"));

  // ------------------------------------------------------------ reading
  // Walk the page in document order: date headings and "New Entries:" style
  // labels set the context for the entry links that follow them.
  function extractRows(doc) {
    const rows = [];
    let day = "";
    let section = "";
    const body = doc.body || doc.documentElement;
    for (const el of body.querySelectorAll("*")) {
      if (el === host || host.contains(el)) continue;
      if (el.tagName === "A") {
        const match = (el.getAttribute("href") || "").match(HREF_RE);
        if (!match || !day || !section) continue;
        const own = (el.textContent || "").trim();
        let prefix = (own.match(/^([#%*+^])/) || [])[1] || "";
        if (!prefix) {
          const prev = el.previousSibling;
          const text = prev ? (prev.textContent || "").trim() : "";
          prefix = /[#%*+^]$/.test(text) ? text.slice(-1) : "";
        }
        rows.push({ date: day, section, mim: match[2], prefix, title: titleFor(el), link: match[1] });
        continue;
      }
      if (el.children.length > 3) continue;
      const text = (el.textContent || "").replace(/\s+/g, " ").trim();
      if (text.length > 60) continue;
      if (DATE_RE.test(text)) { day = text; section = ""; }
      else if (SECTION_RE.test(text)) section = text;
    }
    return rows;
  }

  // The entry title sits in the same table row, or in the element right
  // after the link (or after the link's cell), depending on the layout.
  function titleFor(link) {
    const clean = node => (node && node.textContent || "").replace(/\s+/g, " ").trim();
    const strip = text => text.replace(/^[#%*+^]?\s*\d{6}\s*/, "").trim();
    const tr = link.closest("tr");
    if (tr) {
      const text = strip(clean(tr));
      if (text && text.length <= 400) return text;
    }
    let node = link;
    for (let depth = 0; depth < 3 && node; depth += 1) {
      const own = strip(clean(node));
      if (own && own.length <= 400) return own;
      const next = node.nextElementSibling;
      if (next && !next.querySelector("a[href*='/entry/'], a[href*='/clinicalSynopsis/']")) {
        const text = clean(next);
        if (text && text.length <= 400) return text;
      }
      node = node.parentElement;
    }
    return "";
  }

  async function sendMonth(year, month, doc) {
    const rows = extractRows(doc);
    if (!rows.length) return { ok: false, error: "頁面讀不到任何更新項目" };
    return send({ type: "updates-month", payload: { year, month, rows } });
  }

  // ------------------------------------------------------------- panels
  let status = null;

  function coverageHtml() {
    if (!status) return "";
    const missing = status.missing || [];
    return `
      <div class="row">已連續涵蓋：${status.coverage_start ? `<b>${esc(status.coverage_start)}</b> 至上個月` : "尚無"}${status.history_complete ? "（歷史已補齊）" : ""}
        <br><span class="muted">最新一筆：${esc(status.last_listed_date || "—")}</span></div>
      <div class="row">${missing.length
        ? `還要擷取：<b>${missing.length}</b> 個月份<br><span class="muted">${esc(missing.slice(0, 6).join("、"))}${missing.length > 6 ? "…" : ""}（只擷取已結束的月份）</span>`
        : `已擷取到 <b>${esc(status.last_full_month || "上個月")}</b>，本月結束後再來擷取。`}</div>`;
  }

  async function loadStatus() {
    const result = await send({ type: "updates-status" });
    if (!result.ok) {
      $("#body").innerHTML = `<span class="badge err">${esc(result.error)}</span>`;
      return null;
    }
    status = result.data;
    if (!status.db_enabled) {
      $("#body").innerHTML = `<span class="muted">平台 OMIM 資料庫尚未啟用。</span>`;
      return null;
    }
    return status;
  }

  async function monthPage() {
    const year = Number(PAGE[1]);
    const month = Number(PAGE[2]);
    await loadStatus();
    const target = ym(year, month);
    if (!status) return;
    const stored = (status.months || []).find(m => m.month === target);
    $("#body").innerHTML = `<div class="row"><b>${esc(target)}</b>：${
      stored ? `平台已有 ${esc(stored.rows)} 筆${stored.complete ? "" : "（尚未完整）"}` : "平台還沒有這個月份"}</div>`
      + coverageHtml();
    $("#foot").innerHTML = `<button class="btn" id="go">擷取這個月份</button>`;
    $("#go").addEventListener("click", async () => {
      $("#go").disabled = true;
      $("#go").textContent = "擷取中…";
      const result = await sendMonth(year, month, document);
      $("#go").disabled = false;
      $("#go").textContent = "擷取這個月份";
      $("#body").insertAdjacentHTML("afterbegin", result.ok
        ? `<div class="row"><span class="badge ok">已送到平台：${esc(result.data.rows)} 筆</span></div>`
        : `<div class="row"><span class="badge err">${esc(result.error)}</span></div>`);
    });
  }

  let stopRequested = false;

  async function overviewPage() {
    if (!await loadStatus()) return;
    $("#body").innerHTML = coverageHtml();
    const missing = status.missing || [];
    // Months stored before titles were read: fetch them again to fill titles.
    const untitled = (status.untitled || []).filter(m => !missing.includes(m));
    if (!missing.length && !untitled.length) {
      $("#foot").innerHTML = `<span class="muted">需要的月份都擷取過了。</span>`;
      return;
    }
    $("#foot").innerHTML = `
      ${missing.length ? `<button class="btn" id="go">擷取缺少的月份（${missing.length}）</button>` : ""}
      ${untitled.length ? `<button class="btn${missing.length ? " secondary" : ""}" id="titles" title="這些月份是在讀取標題之前擷取的，重新擷取以補上標題">補標題（${untitled.length} 個月份）</button>` : ""}
      <button class="btn secondary" id="stop" hidden>停止</button>`;
    if (missing.length) $("#go").addEventListener("click", () => fetchMissing(missing));
    if (untitled.length) $("#titles").addEventListener("click", () => fetchMissing(untitled));
    $("#stop").addEventListener("click", () => { stopRequested = true; $("#stop").disabled = true; });
  }

  async function fetchMissing(months) {
    stopRequested = false;
    root.querySelectorAll("#go, #titles").forEach(b => { b.disabled = true; });
    $("#stop").hidden = false;
    $("#stop").disabled = false;
    $("#body").insertAdjacentHTML("beforeend", `
      <div class="row"><progress id="prog" max="${months.length}" value="0"></progress>
        <span class="muted" id="prog-text"></span></div>`);
    let done = 0;
    let stopReason = "";
    for (const month of months) {
      if (stopRequested) { stopReason = "已停止，下次會從剩下的月份繼續。"; break; }
      const [y, m] = month.split("-").map(Number);
      $("#prog-text").textContent = `讀取 ${month}…（${done}/${months.length}）`;
      let response;
      try {
        response = await fetch(`/statistics/updates/${y}/${m}`, { credentials: "include" });
      } catch (err) {
        stopReason = `讀取 ${month} 失敗（${err.message || err}），可能被 OMIM 暫時擋下，請稍後再試。`;
        break;
      }
      if (!response.ok) {
        stopReason = `讀取 ${month} 時 OMIM 回應 ${response.status}，可能被暫時擋下，請稍後再試。`;
        break;
      }
      const doc = new DOMParser().parseFromString(await response.text(), "text/html");
      const result = await sendMonth(y, m, doc);
      if (!result.ok) {
        stopReason = `${month}：${result.error}（可能被 OMIM 暫時擋下，或頁面格式改了）`;
        break;
      }
      done += 1;
      $("#prog").value = done;
      $("#prog-text").textContent = `${month}：${result.data.rows} 筆（${done}/${months.length}）`;
      if (done < months.length) await sleep(DELAY_MS);
    }
    $("#stop").hidden = true;
    await loadStatus();
    const message = stopReason ? `<div class="row"><span class="badge warn">${esc(stopReason)}</span></div>`
      : `<div class="row"><span class="badge ok">完成 ${done} 個月份</span></div>`;
    await overviewPage();
    $("#body").insertAdjacentHTML("afterbegin", message);
  }

  if (PAGE) monthPage();
  else overviewPage();
})();
