// HPO tree picture for one candidate disease: patient terms ↔ disease annotations.
//
// Classic script (no build step, no external library): window.HpoGraph =
//   { open({title, subtitle, load}), render(el, graph), layout(graph), collapse(graph), symbol(row) }.
// `load(allAnnotations)` returns a Promise of the /prioritization/graph payload.
// Layout: organ systems on top, each term one layer below its lowest parent in
// the picture (so is_a edges always point up), barycentre ordering; chains of
// ≥2 plain intermediate terms collapse into one "⋯ n 層 ⋯" node.
(function () {
  "use strict";

  const NODE_W = 168, CHAIN_W = 96, NODE_H = 40, GAP_X = 16, GAP_Y = 64, PAD = 16;
  const ROLE_ZH = { patient: "病人的 HPO", annotation: "疾病註記（有配對）", other_annotation: "疾病其他註記",
    meet: "交會的 term（兩邊往上相遇處）", organ: "器官系統" };

  const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));

  function symbol(r) {
    switch (r && r.relation) {
      case "same": return "＝";
      case "patient_specific": return `↑${r.steps_patient}`;
      case "patient_broader": return `↓${r.steps_annotation}`;
      case "collateral": return `⋀${r.steps_patient}+${r.steps_annotation}`;
      default: return "✕";
    }
  }

  function symbolText(r) {
    switch (r && r.relation) {
      case "same": return "相同";
      case "patient_specific": return `病人較具體：往上 ${r.steps_patient} 層就是註記`;
      case "patient_broader": return `病人較籠統：註記在下面 ${r.steps_annotation} 層`;
      case "collateral": return `旁系：病人往上 ${r.steps_patient} 層、註記往上 ${r.steps_annotation} 層後相遇`;
      default: return "沒有相關註記";
    }
  }

  // ── collapse plain chains ────────────────────────────────────────────────
  function collapse(graph) {
    const nodes = new Map((graph.nodes || []).map((n) => [n.id, { ...n, roles: n.roles || [] }]));
    const parents = new Map(), children = new Map();
    for (const id of nodes.keys()) { parents.set(id, []); children.set(id, []); }
    for (const [c, p] of graph.edges || []) {
      if (!nodes.has(c) || !nodes.has(p)) continue;
      parents.get(c).push(p); children.get(p).push(c);
    }
    const plain = (id) => nodes.get(id).roles.length === 0 && parents.get(id).length === 1 && children.get(id).length === 1;
    const alias = new Map();
    const out = new Map();
    const done = new Set();
    for (const id of [...nodes.keys()].sort()) {
      if (done.has(id) || !plain(id)) continue;
      // walk down to the lowest plain node, then collect upwards
      let low = id;
      while (plain(children.get(low)[0]) && !done.has(children.get(low)[0])) low = children.get(low)[0];
      const chain = [low];
      while (plain(parents.get(chain[chain.length - 1])[0]) && !done.has(parents.get(chain[chain.length - 1])[0])) {
        chain.push(parents.get(chain[chain.length - 1])[0]);
      }
      chain.forEach((x) => done.add(x));
      if (chain.length < 2) continue;
      const vid = `chain:${chain[0]}`;
      out.set(vid, { id: vid, chain: chain.map((x) => nodes.get(x)), roles: ["chain"], name: `⋯ ${chain.length} 層 ⋯` });
      chain.forEach((x) => alias.set(x, vid));
    }
    for (const [id, n] of nodes) if (!alias.has(id)) out.set(id, n);
    const map = (id) => alias.get(id) || id;
    const edgeSet = new Set();
    const edges = [];
    for (const [c, p] of graph.edges || []) {
      if (!nodes.has(c) || !nodes.has(p)) continue;
      const a = map(c), b = map(p);
      const key = `${a}>${b}`;
      if (a === b || edgeSet.has(key)) continue;
      edgeSet.add(key); edges.push([a, b]);
    }
    return { nodes: out, edges, alias: map };
  }

  // ── layered layout ───────────────────────────────────────────────────────
  function layout(graph) {
    const g = collapse(graph);
    const parents = new Map([...g.nodes.keys()].map((k) => [k, []]));
    const children = new Map([...g.nodes.keys()].map((k) => [k, []]));
    for (const [c, p] of g.edges) { parents.get(c).push(p); children.get(p).push(c); }
    const layer = new Map();
    const depth = (id, seen = new Set()) => {
      if (layer.has(id)) return layer.get(id);
      if (seen.has(id)) return 0;                       // defensive: the HPO is a DAG
      seen.add(id);
      const ps = parents.get(id);
      const v = ps.length ? 1 + Math.max(...ps.map((p) => depth(p, seen))) : 0;
      layer.set(id, v);
      return v;
    };
    [...g.nodes.keys()].forEach((id) => depth(id));
    const layers = [];
    for (const [id, l] of layer) (layers[l] ||= []).push(id);
    layers.forEach((ids) => ids.sort((a, b) => String(g.nodes.get(a).name).localeCompare(String(g.nodes.get(b).name))));
    const pos = new Map();
    const setPos = () => layers.forEach((ids) => ids.forEach((id, i) => pos.set(id, i)));
    setPos();
    const bary = (id, nbrs) => {
      const xs = nbrs.get(id).map((n) => pos.get(n)).filter((x) => x !== undefined);
      return xs.length ? xs.reduce((a, b) => a + b, 0) / xs.length : pos.get(id);
    };
    for (let it = 0; it < 6; it++) {
      for (let l = 1; l < layers.length; l++) {
        layers[l].sort((a, b) => bary(a, parents) - bary(b, parents) || a.localeCompare(b)); setPos();
      }
      for (let l = layers.length - 2; l >= 0; l--) {
        layers[l].sort((a, b) => bary(a, children) - bary(b, children) || a.localeCompare(b)); setPos();
      }
    }
    const width = (id) => (g.nodes.get(id).roles.includes("chain") ? CHAIN_W : NODE_W);
    const rowW = layers.map((ids) => ids.reduce((s, id) => s + width(id), 0) + GAP_X * Math.max(ids.length - 1, 0));
    const total = Math.max(...rowW, NODE_W);
    const xy = new Map();
    layers.forEach((ids, l) => {
      let x = PAD + (total - rowW[l]) / 2;
      ids.forEach((id) => {
        xy.set(id, { x, y: PAD + l * (NODE_H + GAP_Y), w: width(id), h: NODE_H });
        x += width(id) + GAP_X;
      });
    });
    return { ...g, xy, layers, width: total + PAD * 2, height: PAD * 2 + layers.length * (NODE_H + GAP_Y) - GAP_Y };
  }

  // ── rendering ────────────────────────────────────────────────────────────
  function nodeClass(n) {
    const r = new Set(n.roles);
    if (r.has("chain")) return "hpog-chain";
    if (r.has("patient") && r.has("annotation")) return "hpog-both";
    if (r.has("patient")) return "hpog-patient";
    if (r.has("annotation")) return "hpog-annotation";
    if (r.has("other_annotation")) return "hpog-other";
    if (r.has("meet")) return "hpog-meet";
    if (r.has("organ")) return "hpog-organ";
    return "hpog-plain";
  }

  const LINE_CHARS = 25, MAX_LINES = 2;

  // Organ systems are all "Abnormality of (the) X": show X so the meaningful word fits.
  function organName(name) {
    const short = String(name || "").replace(/^Abnormality of (the )?/i, "");
    return short ? short[0].toUpperCase() + short.slice(1) : String(name || "");
  }

  // Up to two lines broken at spaces; only a third line's worth is cut with "…".
  function wrap(text, width = LINE_CHARS, maxLines = MAX_LINES) {
    const words = String(text).split(/\s+/).filter(Boolean);
    const lines = [];
    let cur = "";
    for (const w of words) {
      if (!cur) cur = w;
      else if ((cur + " " + w).length <= width) cur += " " + w;
      else { lines.push(cur); cur = w; }
    }
    if (cur) lines.push(cur);
    const out = lines.slice(0, maxLines).map((l) => (l.length > width ? `${l.slice(0, width - 1)}…` : l));
    if (lines.length > maxLines) {
      const last = out[maxLines - 1];
      out[maxLines - 1] = `${last.slice(0, width - 1).replace(/\s+$/, "")}…`;
    }
    return out;
  }

  function label(n) {
    if (n.roles.includes("chain")) return [n.name];
    const r = n.roles;
    const text = r.includes("organ") && !r.includes("patient") && !r.includes("annotation")
      ? organName(n.name || n.id) : String(n.name || n.id);
    const star = r.includes("meet") && !r.includes("patient") ? "★ " : "";
    return wrap(star + text);
  }

  function textSvg(lines, cx, cy) {
    const lh = 13, top = cy + 4 - ((lines.length - 1) * lh) / 2;
    return `<text text-anchor="middle">${lines.map((l, i) =>
      `<tspan x="${cx}" y="${top + i * lh}">${esc(l)}</tspan>`).join("")}</text>`;
  }

  function tooltip(n) {
    if (n.roles.includes("chain")) return n.chain.map((c) => `${c.name}（${c.id}）`).join("\n");
    const roles = n.roles.map((r) => ROLE_ZH[r]).filter(Boolean).join("、");
    return [`${n.name}（${n.id}）`, `IC ${Number(n.ic).toFixed(2)}；${n.n_diseases} 個疾病帶有此 term 或其下層`, roles]
      .filter(Boolean).join("\n");
  }

  function pairPaths(graph, L) {
    // node and edge keys (after collapsing) of each patient term's pairing
    const out = new Map();
    for (const p of graph.pairs || []) {
      const nodes = new Set(), edges = new Set();
      for (const path of [p.path_patient || [], p.path_annotation || []]) {
        const ids = path.map(L.alias);
        ids.forEach((id) => nodes.add(id));
        for (let i = 0; i + 1 < ids.length; i++) if (ids[i] !== ids[i + 1]) edges.add(`${ids[i]}>${ids[i + 1]}`);
      }
      out.set(p.hpo_id, { nodes, edges });
    }
    return out;
  }

  function render(el, graph) {
    const L = layout(graph);
    const paths = pairPaths(graph, L);
    const edgeSvg = L.edges.map(([c, p]) => {
      const a = L.xy.get(c), b = L.xy.get(p);
      const x1 = a.x + a.w / 2, y1 = a.y, x2 = b.x + b.w / 2, y2 = b.y + b.h, my = (y1 + y2) / 2;
      return `<path class="hpog-edge" data-edge="${esc(`${c}>${p}`)}" d="M${x1},${y1} C${x1},${my} ${x2},${my} ${x2},${y2}"/>`;
    }).join("");
    const nodeSvg = [...L.nodes.values()].map((n) => {
      const { x, y, w, h } = L.xy.get(n.id);
      const unmatched = n.roles.includes("patient") && (graph.pairs || []).some((p) => p.hpo_id === n.id && p.relation === "none");
      return `<g class="hpog-node ${nodeClass(n)}${unmatched ? " hpog-unmatched" : ""}" data-id="${esc(n.id)}">
        <title>${esc(tooltip(n))}</title>
        <rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${n.roles.includes("chain") ? h / 2 : 6}"/>
        ${textSvg(label(n), x + w / 2, y + h / 2)}</g>`;
    }).join("");
    const order = { match: 0, partial: 1, none: 2 };
    const rows = [...(graph.pairs || [])].sort((a, b) => (order[a.match] ?? 3) - (order[b.match] ?? 3));
    const names = Object.fromEntries((graph.terms || []).map((t) => [t.hpo_id, t.name]));
    const table = rows.map((p) => `<tr data-hpo="${esc(p.hpo_id)}" class="hpog-row-${esc(p.match)}">
        <td>${esc(names[p.hpo_id] || p.hpo_id)} <span class="hpog-id">${esc(p.hpo_id)}</span></td>
        <td>${p.relation === "none" ? "—" : `${esc(p.annotation_name)} <span class="hpog-id">${esc(p.annotation)}</span>`}</td>
        <td><span class="hpog-sym" title="${esc(symbolText(p))}">${esc(symbol(p))}</span> ${esc(p.relation_zh || "")}</td>
        <td>${p.relation === "none" ? "—" : `${esc(p.mica_name)} <span class="hpog-id">IC ${Number(p.ic_mica).toFixed(2)}</span>`}</td>
        <td class="hpog-num">${p.s == null ? "—" : Number(p.s).toFixed(2)}</td>
        <td class="hpog-num">${p.net == null ? "—" : Number(p.net).toFixed(2)}</td></tr>`).join("");
    el.innerHTML = `
      <div class="hpog-canvas"><svg width="${L.width}" height="${L.height}" viewBox="0 0 ${L.width} ${L.height}"
        role="img" aria-label="HPO 樹狀圖">${edgeSvg}${nodeSvg}</svg></div>
      <table class="hpog-table"><thead><tr><th>病人 HPO</th><th>配對的疾病註記</th><th>關係（層數）</th>
        <th>交會的 term</th><th class="hpog-num">配對分</th><th class="hpog-num">淨貢獻</th></tr></thead>
        <tbody>${table}</tbody></table>`;
    const svg = el.querySelector("svg");
    const highlight = (hpo) => {
      const hit = paths.get(hpo);
      svg.classList.toggle("hpog-dim", !!hit);
      svg.querySelectorAll(".hpog-hl").forEach((x) => x.classList.remove("hpog-hl"));
      el.querySelectorAll(".hpog-table tr.hpog-active").forEach((x) => x.classList.remove("hpog-active"));
      if (!hit) return;
      hit.nodes.forEach((id) => svg.querySelector(`.hpog-node[data-id="${CSS.escape(id)}"]`)?.classList.add("hpog-hl"));
      hit.edges.forEach((k) => svg.querySelector(`.hpog-edge[data-edge="${CSS.escape(k)}"]`)?.classList.add("hpog-hl"));
      el.querySelector(`.hpog-table tr[data-hpo="${CSS.escape(hpo)}"]`)?.classList.add("hpog-active");
    };
    el._hpog = { highlight, paths };
    if (!el._hpogBound) {                         // the container is reused across loads
      el._hpogBound = true;
      el.addEventListener("mouseover", (ev) => {
        const st = el._hpog;
        const node = ev.target.closest(".hpog-node");
        const row = ev.target.closest(".hpog-table tr[data-hpo]");
        if (row) st.highlight(row.dataset.hpo);
        else if (node && st.paths.has(node.dataset.id)) st.highlight(node.dataset.id);
      });
      el.addEventListener("mouseleave", () => el._hpog.highlight(null));
    }
    return L;
  }

  // ── modal ────────────────────────────────────────────────────────────────
  let modal = null, reqId = 0, current = null;

  function close() {
    if (!modal) return;
    modal.hidden = true;
    current = null;
    reqId++;
  }

  function build() {
    modal = document.createElement("div");
    modal.className = "hpog-modal";
    modal.hidden = true;
    modal.innerHTML = `<div class="hpog-card" role="dialog" aria-modal="true">
      <div class="hpog-head"><div><div class="hpog-title"></div><div class="hpog-sub muted"></div></div>
        <button type="button" class="btn btn-ghost hpog-close">關閉</button></div>
      <div class="hpog-tools">
        <span class="hpog-legend">
          <span class="hpog-key hpog-patient">病人的 HPO</span>
          <span class="hpog-key hpog-annotation">疾病註記</span>
          <span class="hpog-key hpog-both">兩者相同</span>
          <span class="hpog-key hpog-meet">★ 交會的 term</span>
          <span class="hpog-key hpog-organ">器官系統</span>
          <span class="hpog-key hpog-patient hpog-unmatched">沒有相關註記</span>
        </span>
        <label class="hpog-others"><input type="checkbox"> 顯示疾病其他註記 <span class="hpog-others-n"></span></label>
      </div>
      <div class="hpog-note muted">越上面越籠統、越下面越具體；直系是一條直線，旁系是兩條線在交會的 term 相遇。滑鼠移到病人的 HPO 或下方表格，只顯示那一組配對的路徑。</div>
      <div class="hpog-body"></div></div>`;
    document.body.appendChild(modal);
    modal.addEventListener("click", (ev) => {
      if (ev.target === modal || ev.target.closest(".hpog-close")) close();
    });
    modal.querySelector(".hpog-others input").addEventListener("change", (ev) => {
      if (current) load(ev.target.checked);
    });
    document.addEventListener("keydown", (ev) => {
      if (ev.key === "Escape" && modal && !modal.hidden) {
        ev.stopImmediatePropagation(); ev.preventDefault(); close();
      }
    }, true);
  }

  async function load(all) {
    const body = modal.querySelector(".hpog-body");
    const id = ++reqId;
    body.innerHTML = `<div class="muted hpog-msg">載入中…</div>`;
    try {
      const graph = await current.load(!!all);
      if (id !== reqId) return;
      if (!graph || !graph.nodes) throw new Error("沒有資料");
      render(body, graph);
      const n = modal.querySelector(".hpog-others-n");
      n.textContent = graph.others_total ? `（${all ? `${graph.others_shown}／` : ""}${graph.others_total}）` : "（0）";
    } catch (err) {
      if (id === reqId) body.innerHTML = `<div class="hpog-msg">無法載入：${esc(err && err.message || err)}</div>`;
    }
  }

  function open(opts) {
    if (!modal) build();
    current = opts;
    modal.querySelector(".hpog-title").textContent = opts.title || "HPO 樹狀圖";
    modal.querySelector(".hpog-sub").textContent = opts.subtitle || "";
    modal.querySelector(".hpog-others input").checked = false;
    modal.querySelector(".hpog-others-n").textContent = "";
    modal.hidden = false;
    load(false);
  }

  window.HpoGraph = { open, close, render, layout, collapse, symbol, symbolText, label, wrap };
})();
