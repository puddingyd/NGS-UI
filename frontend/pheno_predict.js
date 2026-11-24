// Phenotype prediction popup (plan §4.3–4.8 engine), shared by the main page
// and /phenotype/.  PhenoPredict.open([{phenotype, weight}]) posts the HPO
// terms to /api/pheno-engine/predict and shows the top 10 genes and diseases.
(function () {
  "use strict";

  const API = "/api/pheno-engine/predict";
  let modal = null;
  let request = 0;

  function esc(value) {
    return String(value == null ? "" : value)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;").replace(/'/g, "&#39;");
  }

  function num(value, digits) {
    const x = Number(value);
    return Number.isFinite(x) ? x.toFixed(digits) : "—";
  }

  function pct(value) {
    const x = Number(value);
    return Number.isFinite(x) ? `${Math.round(x * 100)}%` : "—";
  }

  function diseaseLink(id) {
    const text = String(id || "");
    let url = "";
    if (text.startsWith("OMIM:")) url = `https://omim.org/entry/${encodeURIComponent(text.slice(5))}`;
    else if (text.startsWith("ORPHA:")) url = `https://www.orpha.net/en/disease/detail/${encodeURIComponent(text.slice(6))}`;
    return url ? `<a href="${url}" target="_blank" rel="noopener">${esc(text)}</a>` : esc(text);
  }

  function ensureModal() {
    if (modal) return modal;
    modal = document.createElement("div");
    modal.className = "ppred-modal";
    modal.hidden = true;
    modal.setAttribute("role", "dialog");
    modal.setAttribute("aria-modal", "true");
    modal.setAttribute("aria-labelledby", "ppred-title");
    modal.innerHTML = `
      <div class="ppred-card">
        <div class="ppred-head">
          <h3 id="ppred-title">依 HPO 預測致病基因與疾病</h3>
          <button type="button" class="ppred-close" aria-label="關閉">×</button>
        </div>
        <div class="ppred-body"></div>
      </div>`;
    document.body.appendChild(modal);
    modal.addEventListener("click", (event) => {
      if (event.target === modal || event.target.closest(".ppred-close")) close();
      const row = event.target.closest(".ppred-disease-row");
      if (row && !event.target.closest("a")) {
        const detail = modal.querySelector(`.ppred-explain[data-for="${row.dataset.id}"]`);
        if (detail) {
          detail.hidden = !detail.hidden;
          row.classList.toggle("open", !detail.hidden);
        }
      }
    });
    document.addEventListener("keydown", (event) => {
      if (event.key === "Escape" && modal && !modal.hidden) close();
    });
    return modal;
  }

  function close() {
    if (!modal) return;
    modal.hidden = true;
    request += 1;
  }

  function renderTerms(data) {
    const kept = (data.terms || []).map((t) => {
      const w = Number(t.weight) !== 1 ? `（權重 ${esc(t.weight)}）` : "";
      const via = t.how !== "current" ? `（由 ${esc(t.input_id)} 換算）` : "";
      return `<span class="ppred-chip">${esc(t.name)} <span class="ppred-id">${esc(t.hpo_id)}</span>${via}${w}</span>`;
    }).join("");
    const removed = (data.removed || []).map((r) =>
      `<li><span class="ppred-id">${esc(r.hpo_id)}</span> ${esc(r.reason_zh)}</li>`).join("");
    return `
      <div class="ppred-section">
        <div class="ppred-label">計分的 HPO（${(data.terms || []).length} 個）</div>
        <div class="ppred-chips">${kept || "—"}</div>
        ${removed ? `<details class="ppred-removed"><summary>未計分 ${(data.removed || []).length} 個</summary><ul>${removed}</ul></details>` : ""}
      </div>`;
  }

  function renderGenes(data) {
    if (!(data.genes || []).length) {
      return `<p class="ppred-muted">這個知識庫沒有疾病–基因對應，只能排疾病。</p>`;
    }
    const rows = data.genes.map((g) => {
      const tags = [];
      if (g.lead) tags.push(`<span class="ppred-tag ppred-tag-lead" title="第 1 名且領先第 2 名 ≥3 nats">明顯領先</span>`);
      else if (g.near_tie && g.rank > 1) tags.push(`<span class="ppred-tag" title="與第 1 名相差 ≤1 nat">接近同分</span>`);
      return `<tr>
        <td class="ppred-num">${esc(g.rank)}</td>
        <td class="ppred-gene">${esc(g.gene)} ${tags.join("")}</td>
        <td class="ppred-num">${num(g.G, 2)}</td>
        <td class="ppred-num">${Number(g.P) > 0 ? "+" : ""}${num(g.P, 1)}</td>
        <td>${diseaseLink(g.best_disease)} ${esc(g.best_disease_name)}</td>
        <td class="ppred-num">${pct(g.explained)}</td>
      </tr>`;
    }).join("");
    return `<table class="ppred-table">
      <thead><tr><th>#</th><th>基因</th><th title="最佳疾病校正後分數 S*">分數</th>
        <th title="表型點數 P（暫定校準）">P</th><th>最吻合的疾病</th><th title="疾病解釋了病人多少表型資訊">解釋</th></tr></thead>
      <tbody>${rows}</tbody></table>`;
  }

  function renderDiseases(data) {
    const rows = (data.diseases || []).map((d) => {
      const id = esc(d.disease_id);
      const explain = (d.explain || []).map((r) => `<tr>
          <td>${esc(r.hpo_name)} <span class="ppred-id">${esc(r.hpo_id)}</span></td>
          <td>${esc(r.annotation_name)} <span class="ppred-id">${esc(r.annotation)}</span></td>
          <td>${esc(r.relation_zh)}</td>
          <td class="ppred-num">${num(r.s, 2)}</td>
          <td class="ppred-num">${num(r.net, 2)}</td>
        </tr>`).join("");
      return `<tr class="ppred-disease-row" data-id="${id}" title="點一下看每個 HPO 怎麼配對">
          <td class="ppred-num">${esc(d.rank)}</td>
          <td>${diseaseLink(d.disease_id)}</td>
          <td>${esc(d.name)}</td>
          <td class="ppred-num">${num(d.S_star, 2)}</td>
          <td class="ppred-num">${pct(d.explained)}</td>
          <td>${esc((d.genes || []).slice(0, 4).join("、")) || "—"}</td>
        </tr>
        <tr class="ppred-explain" data-for="${id}" hidden><td colspan="6">
          <table class="ppred-subtable"><thead><tr><th>病人 HPO</th><th>疾病最吻合的註記</th><th>關係</th>
            <th>配對分</th><th title="配對分扣掉這個疾病的基準分">淨貢獻</th></tr></thead><tbody>${explain}</tbody></table>
        </td></tr>`;
    }).join("");
    return `<table class="ppred-table">
      <thead><tr><th>#</th><th>疾病</th><th>名稱</th><th title="校正後分數 S*">分數</th><th>解釋</th><th>基因</th></tr></thead>
      <tbody>${rows}</tbody></table>`;
  }

  function render(data) {
    const body = modal.querySelector(".ppred-body");
    body.innerHTML = `
      ${renderTerms(data)}
      <div class="ppred-section">
        <div class="ppred-label">基因排名前 ${(data.genes || []).length} 名</div>
        ${renderGenes(data)}
      </div>
      <div class="ppred-section">
        <div class="ppred-label">疾病排名前 ${(data.diseases || []).length} 名 <span class="ppred-muted">（點一列可看每個 HPO 的配對）</span></div>
        ${renderDiseases(data)}
      </div>
      <p class="ppred-note">知識庫 ${esc(data.snapshot)}（${Number(data.n_diseases || 0).toLocaleString()} 個疾病、${Number(data.n_genes || 0).toLocaleString()} 個基因）。
        分數為研究用途的表型排序，不是診斷；${esc(data.calibration || "")}。目前不影響變異卡片的排序。</p>`;
  }

  async function open(items) {
    ensureModal();
    const terms = (items || []).filter((t) => t && t.phenotype);
    const body = modal.querySelector(".ppred-body");
    modal.hidden = false;
    if (!terms.length) {
      body.innerHTML = `<p class="ppred-muted">還沒有選任何 HPO term。</p>`;
      return;
    }
    const token = ++request;
    body.innerHTML = `<p class="ppred-muted">計算中…（第一次使用需載入知識庫，約數秒）</p>`;
    try {
      const resp = await fetch(API, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        credentials: "same-origin",
        body: JSON.stringify({ hpo: terms, top: 10 }),
      });
      const data = await resp.json().catch(() => ({}));
      if (token !== request) return;
      if (!resp.ok) throw new Error(data.detail || resp.statusText || "預測失敗");
      render(data);
    } catch (err) {
      if (token !== request) return;
      body.innerHTML = `<p class="ppred-error">${esc(err && err.message ? err.message : err)}</p>`;
    }
  }

  window.PhenoPredict = { open, close };
})();
