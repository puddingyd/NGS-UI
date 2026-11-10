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

  async function loadCaptures() {
    const { captures } = await api("/api/omim/captures?limit=100");
    $("#capture-table tbody").innerHTML = captures.length ? captures.map(c => `
      <tr class="clickable" data-capture="${esc(c.id)}">
        <td>${esc(c.prefix || "")}${esc(c.mim)}</td>
        <td>${esc(c.title)}</td>
        <td>${esc(c.description_source || "—")}</td>
        <td>${badge(SYNOPSIS[c.synopsis_status] || [c.synopsis_status || "—", "info"])}</td>
        <td>${esc(c.omim_edit_date || "—")}</td>
        <td>${esc(fmtTime(c.captured_at))}<br><span class="muted-text">${esc(c.user)}${c.duplicate_of ? "（內容相同）" : ""}</span></td>
        <td>${(c.warnings || []).length ? badge([`${c.warnings.length} 項`, "warn"]) : ""}</td>
      </tr>`).join("")
      : `<tr><td colspan="7" class="muted-text">還沒有擷取紀錄。</td></tr>`;
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
      ${(summary.disease_cells || []).map(c => `<p><b>${esc(c.gene)} Disease 欄預覽</b></p><pre>${esc(c.cell)}</pre>`).join("")}`;
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
  $("#btn-refresh").addEventListener("click", () => loadCaptures());

  (async () => {
    try {
      const me = await api("/api/auth/me");
      $("#whoami").textContent = me.username;
      await Promise.all([loadTokens(), loadCaptures()]);
    } catch (_) { /* login hint already shown */ }
  })();
})();
