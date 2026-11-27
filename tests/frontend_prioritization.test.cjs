const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const test = require('node:test');
const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
function fn(name) {
  const start = source.indexOf(`function ${name}(`);
  assert.ok(start >= 0, name);
  return source.slice(start, source.indexOf('\n}\n', start) + 3);
}
function setup(data) {
  const context = vm.createContext({ state: { data } });
  const start = source.indexOf('const PRIO_MODE_ZH =');
  const end = source.indexOf('};', start) + 2;
  const k = source.indexOf('const PRIO_KIND_LABEL =');
  vm.runInContext(source.slice(start, end) + source.slice(k, source.indexOf(';', k) + 1) + fn('_prioInfo') +
    fn('_prioNum') + fn('_prioSigned') + fn('_prioScoreText') +
    fn('_prioHintLines') + fn('_prioMode') + fn('_snvSortScore') +
    '\nconst displaySnvHgvs = (v, id) => v.hgvs || id;\n', context);
  return context;
}

test('score text: total (V + P + C) (rank N), negatives and missing P', () => {
  const c = setup({});
  assert.equal(c._prioScoreText({ prio: { total: 17.84, v: 10, p: 7.84, c: 0, rank: 1 } }), '17.8 (10 + 7.8 + 0) (rank 1)');
  assert.equal(c._prioScoreText({ prio: { total: 6.4, v: 10, p: 0.4, c: -4, rank: 3 } }), '6.4 (10 + 0.4 − 4) (rank 3)');
  assert.equal(c._prioScoreText({ prio: { total: 4, v: 4, p: null, c: 0 } }), '4 (4 + — + 0)');
  assert.equal(c._prioScoreText({ prio: { total: -1.2, v: -2, p: -3.1, c: 3.9, rank: 9 } }), '-1.2 (-2 − 3.1 + 3.9) (rank 9)');
  assert.equal(c._prioScoreText({}), '—');
});

test('hint lines explain pair partner, disease and inheritance', () => {
  const c = setup({ prioritization: { mode: 'prio' }, variants: { b: { hgvs: 'c.1A>G' } } });
  const lines = c._prioHintLines({ prio: { total: 8, v: 3.9, p: 5.1, c: -1, mode: 'pair', partner: 'b',
    v_single: 4, partner_v: 10, pheno_counted: true, gene: 'GENEA', gene_rank: 4, gene_points: 7.5,
    disease_id: 'OMIM:1', disease_name: 'Disease', disease_rank: 3, moi: ['ar'], moi_source: 'OMIM' } });
  const text = lines.join('\n');
  assert.match(text, /c\.1A>G/);
  assert.match(text, /表型排名第 4/);
  assert.match(text, /Disease（OMIM:1，疾病第 3 名）/);
  assert.match(text, /phase 未知/);
  assert.match(text, /tier 內依此總分排序/);
});

test('hint lines explain the variant-type prior', () => {
  const c = setup({ prioritization: { mode: 'prio' } });
  const text = c._prioHintLines({ prio: { total: 10.5, v: 0.3, p: 10.2, c: 0, mode: 'dominant', pheno_counted: true,
    gene: 'KCNQ2', moi: ['dominant'], type_prior: { class: 'noncoding', label: '非編碼區（UTR、內含子 ≥+7／≤−21）', points: -3.3, splice: 0.02,
    source: 'ClinVar 2026-11-22' } } }).join('\n');
  assert.match(text, /含變異類型先驗/);
  assert.match(text, /變異類型先驗 -3\.3：非編碼區（UTR、內含子 ≥\+7／≤−21），剪接預測 0\.02/);
  const none = c._prioHintLines({ prio: { v: 2, p: 1, c: 0, mode: 'dominant', pheno_counted: true, moi: [] } }).join('\n');
  assert.doesNotMatch(none, /變異類型先驗/);
});

test('panel-only cases sort by Score (old)', () => {
  const old = setup({ prioritization: { mode: 'old', reason: '沒有 HPO' } });
  assert.equal(old._snvSortScore({ total_score: 80, prio: {}, prio_total: 12 }), 80);
  const text = old._prioHintLines({ prio: { v: 4, c: 0, mode: 'dominant', pheno_counted: false, moi: [] } }).join('\n');
  assert.match(text, /表型 P —：沒有 HPO/);
  assert.match(text, /排序沿用 Score \(old\)/);
  const prio = setup({ prioritization: { mode: 'prio' } });
  assert.equal(prio._snvSortScore({ total_score: 80, prio: {}, prio_total: 12 }), 12);
  assert.equal(prio._snvSortScore({ total_score: 80 }), 80);
});

test('candidate header: score, disease with rank, inheritance mode; button disabled without HPO', () => {
  const start = source.indexOf('const PRIO_MODE_ZH =');
  const end = source.indexOf('};', start) + 2;
  const context = vm.createContext({
    state: { data: { prioritization: { mode: 'old', reason: '沒有 HPO' } } },
    document: { getElementById: () => button },
  });
  const button = { disabled: false, title: '' };
  const s2 = source.indexOf('const PRIO_MODE_SHORT =');
  const e2 = source.indexOf('};', s2) + 2;
  vm.runInContext(source.slice(start, end) + source.slice(s2, e2) + fn('_prioModeShort') + fn('_prioNum') + fn('_prioSigned') + fn('_prioCandidateSummary') +
    fn('updatePrioCandidatesButton') + fn('_prioDiseaseLink'), context);
  const s = context._prioCandidateSummary({ total: 13.54, v: 6, p: 8.5, c: -1, mode: 'pair', moi_class: 'ar',
    disease_id: 'OMIM:100003', disease_name: 'Broad syndrome C', disease_rank: 2 });
  assert.equal(s.score, '13.5 (6 + 8.5 − 1)');
  assert.equal(s.disease, 'Broad syndrome C（OMIM:100003，疾病第 2 名）');
  assert.equal(s.mode, 'AR，兩個變異');
  assert.equal(context._prioModeShort({ mode: 'single', moi_class: 'ar' }), 'AR，一個變異');
  assert.equal(context._prioModeShort({ mode: 'dominant', moi_class: 'dominant' }), 'AD');
  assert.equal(context._prioModeShort({ mode: 'unknown', moi_class: 'unknown' }), '遺傳模式未知');
  assert.equal(context._prioCandidateSummary({ total: 4.6, v: 10, p: -5.4, c: 0, mode: 'dominant', moi_class: 'unknown', no_disease: true }).disease, '沒有已知疾病');
  context.updatePrioCandidatesButton();
  assert.equal(button.disabled, true);
  assert.equal(button.title, '沒有 HPO');
  context.state.data.prioritization = { mode: 'prio', candidates: [] };
  context.updatePrioCandidatesButton();
  assert.equal(button.disabled, false);
  assert.equal(context._prioDiseaseLink('ORPHA:5'), 'https://www.orpha.net/en/disease/detail/5');
  assert.equal(context._prioDiseaseLink(''), '');
});


test('CNV and mito hints, rank override', () => {
  const c = setup({ cnv_prioritization: { mode: 'prio' }, mito_prioritization: { mode: 'old', reason: '沒有 HPO' } });
  const ev = { prio: { total: 15.2, v: 10, p: 5.2, c: 0, rank: 4, annotsv_score: 0.99, n_genes: 3, gene: 'GENEA',
    mode: 'dominant', pheno_counted: true, moi: ['dominant'] } };
  assert.equal(c._prioScoreText(ev, 2), '15.2 (10 + 5.2 + 0) (rank 2)');
  const text = c._prioHintLines(ev, 'cnv').join('\n');
  assert.match(text, /AnnotSV ranking score 0\.99/);
  assert.match(text, /受影響基因 3 個/);
  assert.match(text, /全部 CNV 中的名次/);
  assert.match(text, /與 SNV 配成隱性一對的缺失會列入/);
  const mito = c._prioHintLines({ prio: { total: 7, v: 7, p: null, c: 0, mode: 'mito', v_basis: 'ClinVar P/LP',
    pheno_counted: false, moi: [] } }, 'mito').join('\n');
  assert.match(mito, /ClinVar P\/LP（ClinVar P\/LP＋MITOMAP/);
  assert.match(mito, /表型 P —：沒有 HPO/);
  assert.match(mito, /沿用原本排序/);
});

test('candidate list: CNV deletion pairs with the strongest het SNV of a recessive gene', () => {
  const data = {
    prioritization: { mode: 'prio', candidate_count: 2, candidates: [
      { gene: 'GENEA', total: 9, v: 4, p: 5, c: 0, mode: 'dominant', rank: 1, variants: ['a'], others: [], n_combinations: 1 },
      { gene: 'GENED', total: 2, v: 6, p: 0, c: -4, mode: 'single', rank: 2, variants: ['d'], others: [], n_combinations: 2 },
    ] },
    variants: {
      a: { zygosity: 'het', prio: { gene: 'GENEA', v_single: 4, v: 4 } },
      d: { zygosity: 'heterozygous', prio: { gene: 'GENED', v_single: 6, v: 6 } },
      d2: { zygosity: 'het', prio: { gene: 'GENED', v_single: 2, v: 2 } },
      h: { zygosity: 'hom', prio: { gene: 'GENEQ', v_single: 9, v: 9 } },
    },
    cnv_variants: {
      e1: { prio: { v: 10, ar_genes: [
        { gene: 'GENED', gene_rank: 7, diseases: [{ disease_id: 'OMIM:4', disease_name: 'Eye', p: 5, moi_class: 'ar' }] },
        { gene: 'GENEQ', diseases: [{ disease_id: 'OMIM:9', p: 9, moi_class: 'ar' }] },
      ] } },
    },
  };
  const context = vm.createContext({ state: { data } });
  vm.runInContext(fn('_prioPairPoints') + fn('_prioCnvPairs') + fn('_prioSameGeneHets') + fn('_prioCandidateList'), context);
  const { list, total } = context._prioCandidateList();
  assert.equal(total, 2);                                      // GENEQ has no het SNV → no row
  const pairV = context._prioPairPoints(6, 10);
  assert.ok(Math.abs(pairV - 6.0) < 0.2);
  assert.equal(list[0].gene, 'GENED');                          // 6.x + 5 > 9 → new top row
  assert.deepEqual(Array.from(list[0].variants), ['d', 'e1']);
  assert.equal(list[0].cnv_pair, true);
  assert.equal(list[0].c, 0);
  assert.equal(list[0].others[0].mode, 'single');              // the old best moves to 其他組合
  assert.equal(list[0].n_combinations, 3);
  assert.deepEqual(Array.from(list[0].same_gene), ['d2']);         // the other het SNV stays visible
  assert.deepEqual(list.map(c => c.rank), [1, 2]);
  assert.equal(data.prioritization.candidates[1].gene, 'GENED'); // payload not mutated
  data.cnv_variants = {};
  assert.equal(context._prioCandidateList().list, data.prioritization.candidates);
});

test('candidate HPO chips: matched first, partial, then unexplained terms', () => {
  const context = vm.createContext({
    escapeHtml: s => String(s).replace(/</g, '&lt;'), escapeAttr: s => String(s).replace(/"/g, '&quot;'),
  });
  vm.runInContext(fn('_prioHpoSymbol') + fn('_prioHpoChips'), context);
  const html = context._prioHpoChips({ matched: 1, partial: 1, rows: [
    { hpo_id: 'HP:3', match: 'none', relation_zh: '無共同上層', net: -0.2 },
    { hpo_id: 'HP:2', match: 'partial', annotation: 'HP:9', annotation_name: 'B', relation: 'collateral', relation_zh: '旁系',
      net: 0.5, steps_patient: 2, steps_annotation: 1, mica_name: 'Root-ish' },
    { hpo_id: 'HP:1', match: 'match', annotation: 'HP:1', annotation_name: 'A', relation: 'same', relation_zh: '相同', net: 2 },
  ] }, [{ hpo_id: 'HP:1', name: 'Seizure' }, { hpo_id: 'HP:2', name: 'Ataxia' }]);
  assert.match(html, /HPO 符合 1／3、部分 1/);
  const order = [...html.matchAll(/prio-hpo-(match|partial|none)"/g)].map(m => m[1]);
  assert.deepEqual(order, ['match', 'partial', 'none']);
  assert.match(html, /<\/b> HP:3</);                    // name missing → ID
  assert.match(html, /疾病沒有相關註記/);
  assert.match(html, />⋀2\+1<\/b> Ataxia/);
  assert.match(html, />＝<\/b> Seizure/);
  assert.match(html, /在「Root-ish」相遇/);
  assert.doesNotMatch(html, /基準分）/);
  assert.equal(context._prioHpoSymbol({ relation: 'patient_specific', steps_patient: 3 }), '↑3');
  assert.equal(context._prioHpoSymbol({ relation: 'patient_broader', steps_annotation: 1 }), '↓1');
});

test('HPO graph: chains collapse, parents always drawn above children', () => {
  const src = fs.readFileSync(path.join(__dirname, '../frontend/hpo_graph.js'), 'utf8');
  const context = vm.createContext({ window: {}, document: {}, CSS: { escape: s => s } });
  vm.runInContext(src, context);
  const G = context.window.HpoGraph;
  const graph = {
    nodes: [
      { id: 'O', name: 'Organ', roles: ['organ'] }, { id: 'A', name: 'a', roles: [] }, { id: 'B', name: 'b', roles: [] },
      { id: 'C', name: 'c', roles: [] }, { id: 'M', name: 'Meet', roles: ['meet'] },
      { id: 'P', name: 'Patient', roles: ['patient'] }, { id: 'D', name: 'Disease', roles: ['annotation'] },
      { id: 'Q', name: 'Patient2', roles: ['patient'] },
    ],
    // O ← A ← B ← C ← M ;  M ← P ; M ← D ; O ← Q
    edges: [['A', 'O'], ['B', 'A'], ['C', 'B'], ['M', 'C'], ['P', 'M'], ['D', 'M'], ['Q', 'O']],
    pairs: [{ hpo_id: 'P', path_patient: ['P', 'M'], path_annotation: ['D', 'M'], relation: 'collateral' }],
  };
  const L = G.layout(graph);
  const chain = [...L.nodes.values()].find(n => n.roles.includes('chain'));
  assert.ok(chain && chain.chain.length === 3, 'A-B-C collapse into one chain node');
  assert.equal(L.alias('B'), chain.id);
  for (const [c, p] of L.edges) assert.ok(L.xy.get(p).y < L.xy.get(c).y, `${p} above ${c}`);
  assert.equal(L.xy.get('O').y, Math.min(...[...L.xy.values()].map(v => v.y)));
  assert.equal(G.symbol({ relation: 'collateral', steps_patient: 1, steps_annotation: 2 }), '⋀1+2');
  assert.equal(G.symbol({ relation: 'none' }), '✕');
  // a single plain node is kept (only chains of ≥2 collapse)
  const L2 = G.layout({ nodes: [{ id: 'O', name: 'o', roles: ['organ'] }, { id: 'X', name: 'x', roles: [] },
    { id: 'P', name: 'p', roles: ['patient'] }], edges: [['X', 'O'], ['P', 'X']] });
  assert.ok(L2.nodes.has('X'));
  // organ systems drop "Abnormality of (the)" and wrap instead of cutting the key words
  const lab = n => [...G.label(n)];
  assert.deepEqual(lab({ name: 'Abnormality of the musculoskeletal system', roles: ['organ'] }), ['Musculoskeletal system']);
  assert.deepEqual(lab({ name: 'Abnormality of blood and blood-forming tissues', roles: ['organ'] }),
    ['Blood and blood-forming', 'tissues']);
  assert.deepEqual(lab({ name: 'Growth abnormality', roles: ['organ'] }), ['Growth abnormality']);
  assert.deepEqual(lab({ name: 'Multiple cafe-au-lait spots', roles: ['annotation'] }), ['Multiple cafe-au-lait', 'spots']);
  assert.deepEqual(lab({ name: 'Scoliosis', roles: ['annotation', 'meet'] }), ['★ Scoliosis']);
  assert.equal(lab({ name: 'a b c d e f g h i j k l m n o p q r s t u v w x y z a b c d e f', roles: [] }).length, 2);
});
