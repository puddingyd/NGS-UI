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
  vm.runInContext(source.slice(start, end) + fn('_prioNum') + fn('_prioSigned') + fn('_prioScoreText') +
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
