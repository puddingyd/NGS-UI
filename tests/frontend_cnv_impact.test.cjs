const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
function fixture(render = false) {
  const state = { currentLIS: 'sample', data: { sample_id: 'sample', cnv_variants: {}, sv_variants: {},
    cnv_categories: {}, sv_categories: {} }, reports: { status: {}, edits: {} } };
  const element = () => ({ innerHTML: '', children: [], classList: { add() {}, toggle() {} },
    appendChild(child) { this.children.push(child); },
    addEventListener(event, callback) { this[event] = callback; } });
  const panels = Object.fromEntries(['CNV-1A', 'CNV-1B', 'SV-2A', 'SV-2B'].map(tier => [tier, { ...element(), dataset: { tier } }]));
  const bar = element();
  const document = { createElement: element,
    getElementById(id) { return id === 'cnv-sv-tab-bar' ? bar : { querySelectorAll() { return Object.values(panels); } }; },
    querySelector(selector) { return panels[selector.match(/data-tier="([^"]+)"/)[1]]; } };
  const c = vm.createContext({ state, document,
    _statusValues: raw => String(raw || '').split(','),
    _cnvSvAcmgClassValue: (id, v) => Number(state.reports.edits[id]?.ACMG_class_sv ?? v.acmg_class),
    renderCnvSvCard: (v, id) => ({ ...element(), variantId: id }),
    panels, bar,
    escapeHtml: value => String(value).replaceAll('<', '&lt;').replaceAll('>', '&gt;'),
  });
  vm.runInContext(source.slice(source.indexOf('const CNV_SV_TIER_ORDER'), source.indexOf('// ---------- Mitochondria tier tabs')), c);
  if (!render) c.renderCnvSvTabBar = () => {};
  return c;
}
function variant(id, category = 'noncoding', changes = {}) {
  const impact = { category, reasons: [{ gene: 'TEST', impact: '純內含子' }], hpo_score: 0, mechanism: 0 };
  return { id, source: 'cnv', CHROM: '1', POS: 100, END: 200, sv_type: 'DEL', in_panel: true,
    acmg_class: 3, ranking_score: 0, cnv_sv_sort_score: 50, genes: [],
    impact_clinical: impact, impact_all: impact, ...changes };
}

for (const tier of ['CNV-1A', 'CNV-1B', 'SV-2A', 'SV-2B']) {
  test(`${tier}: defaults, show all, protected P/LP and reviewer status`, () => {
    const c = fixture();
    assert.equal(c._cnvSvPassesImpact(variant('n'), tier), false);
    assert.equal(c._cnvSvPassesImpact(variant('f', 'functional'), tier), true);
    assert.equal(c._cnvSvPassesImpact(variant('u', 'unknown'), tier), true);
    assert.equal(c._cnvSvPassesImpact(variant('p', 'noncoding', { acmg_class: 4 }), tier), tier.endsWith('A'));
    c.state.reports.status.n = 'C,0';
    assert.equal(c._cnvSvPassesImpact(variant('n'), tier), true);
    c.state.reports.status.n = '0';
    assert.equal(c._cnvSvPassesImpact(variant('n'), tier), false);
    c.state.reports.edits.n = { ACMG_class_sv: '5' };
    assert.equal(c._cnvSvPassesImpact(variant('n'), tier), tier.endsWith('A'));
    delete c.state.reports.edits.n;
    c._cnvSvImpactFilter(tier).noncoding = true;
    assert.equal(c._cnvSvPassesImpact(variant('n'), tier), true);
  });
}

test('Clinical uses matched genes; Pathogenic uses every involved gene', () => {
  const c = fixture();
  const v = variant('a', 'noncoding', { impact_all: { category: 'functional' } });
  assert.equal(c._cnvSvPassesImpact(v, 'CNV-1A'), false);
  assert.equal(c._cnvSvPassesImpact(v, 'CNV-1B'), true);
});

test('all callers and WGS/WES receive the same display rules; no source eligibility changes', () => {
  const c = fixture();
  for (const caller of ['dragen', 'nckuh']) for (const testType of ['WGS', 'WES']) {
    c.state.data.pipeline_type = caller; c.state.data.test_type = testType;
    assert.equal(c._cnvSvPassesImpact(variant('a'), 'CNV-1A'), false);
    assert.equal(c._cnvSvPassesImpact(variant('a', 'functional'), 'SV-2A'), true);
  }
});

test('fixed ordering is unaffected by unrelated scaled ranking scores or rescue provenance', () => {
  const c = fixture(), variants = c.state.data.cnv_variants;
  variants.a = variant('a', 'functional', { ranking_score: 0.1, impact_clinical: { category: 'functional', hpo_score: 60 } });
  variants.b = variant('b', 'functional', { ranking_score: 1, impact_clinical: { category: 'functional', hpo_score: 10 } });
  variants.p = variant('p', 'noncoding', { acmg_class: 5 });
  assert.ok(c._cnvSvCompareImpact('a', 'b', 'CNV-1A') < 0);
  assert.ok(c._cnvSvCompareImpact('p', 'a', 'CNV-1A') < 0);
  variants.a.cnv_sv_sort_score = -100; variants.b.cnv_sv_sort_score = 200;
  variants.b.cnv_rescue = { rule: 'A' };
  assert.ok(c._cnvSvCompareImpact('a', 'b', 'CNV-1A') < 0);
});

test('merged parent uses actual segments, preserves clinical scope and member protection', () => {
  const c = fixture(), variants = c.state.data.cnv_variants;
  variants.a = variant('a', 'noncoding', {
    p_loss: { diseases: ['Disease A'], sources: ['CLN:1'], coords: ['1:100-200'] },
  });
  variants.b = variant('b', 'functional', { POS: 300, END: 400, in_panel: false,
    p_loss: { diseases: ['Disease B'], sources: ['dbVar:2'], coords: ['1:300-400'] } });
  c.state.data.cnv_categories['CNV-1A'] = ['a'];
  const parent = c._cnvSvBuildParent({ member_ids: ['a', 'b'] });
  assert.equal(parent.impact_clinical.category, 'noncoding');
  assert.equal(parent.impact_all.category, 'functional');
  assert.deepEqual(Array.from(parent.p_loss.diseases), ['Disease A', 'Disease B']);
  assert.equal(c._cnvSvIdsForTier('CNV-1A').length, 0);
  assert.equal(c._cnvSvIdsForTier('CNV-1A', false).length, 1);
  c.state.reports.status.a = '1';
  assert.equal(c._cnvSvIdsForTier('CNV-1A').length, 1);
});

test('adjacent merge accepts the small chr7 boundary overlap', () => {
  const c = fixture(), variants = c.state.data.cnv_variants;
  variants.a = variant('a', 'functional', { POS: 73303743, END: 74416985 });
  variants.b = variant('b', 'functional', { POS: 74416499, END: 74727986 });
  const groups = c._cnvSvAdjacentMergeGroups(['a', 'b']);
  assert.equal(groups.length, 1);
  assert.deepEqual(Array.from(groups[0], v => v.id), ['a', 'b']);
});

test('adjacent merge protects 10 kb and 10% overlap boundaries', () => {
  const c = fixture(), variants = c.state.data.cnv_variants;
  const grouped = (first, second) => {
    variants.a = variant('a', 'functional', { POS: first[0], END: first[1] });
    variants.b = variant('b', 'functional', { POS: second[0], END: second[1] });
    return c._cnvSvAdjacentMergeGroups(['a', 'b']).length === 1;
  };
  assert.equal(grouped([100, 100100], [90100, 190100]), true);
  assert.equal(grouped([100, 100100], [90099, 190100]), false);
  assert.equal(grouped([100, 50100], [44100, 94100]), false);
  assert.equal(grouped([100, 100100], [10100, 90100]), false);
  assert.equal(grouped([100, 200], [250200, 250300]), true);
  assert.equal(grouped([100, 200], [250201, 250301]), false);
});

test('toolbar remains usable with zero visible results; hide unknown only when absent', () => {
  const c = fixture();
  c.state.data.cnv_variants.a = variant('a');
  let box = c._cnvSvImpactToolbar('CNV-1A', ['a']);
  assert.match(box.innerHTML, /Only UTR \/ intronic <span>\(1 \/ 1\)<\/span>/);
  assert.doesNotMatch(box.innerHTML, /<legend|<button|顯示的影響類型|臨床優先排序|近似位點/);
  assert.doesNotMatch(box.innerHTML, /data-impact="unknown"/);
  box.change({ target: { dataset: { impact: 'noncoding' }, checked: true } });
  assert.equal(c._cnvSvPassesImpact(variant('a'), 'CNV-1A'), true);
  box.change({ target: { dataset: { impact: 'noncoding' }, checked: false } });
  assert.equal(c._cnvSvPassesImpact(variant('a'), 'CNV-1A'), false);
  c.state.data.cnv_variants.u = variant('u', 'unknown');
  box = c._cnvSvImpactToolbar('CNV-1A', ['a', 'u']);
  assert.match(box.innerHTML, /data-impact="unknown" checked/);
  box.change({ target: { dataset: { impact: 'unknown' }, checked: false } });
  assert.equal(c._cnvSvPassesImpact(variant('u', 'unknown'), 'CNV-1A'), false);
  c.state.data.sample_id = 'new-sample';
  assert.equal(c._cnvSvPassesImpact(variant('u', 'unknown'), 'CNV-1A'), true);
});

test('protected near-duplicates remain visible and card evidence is escaped', () => {
  const c = fixture();
  c.state.data.cnv_variants.a = variant('a', 'functional');
  c.state.data.cnv_variants.b = variant('b', 'functional');
  assert.equal(c._cnvSvClusterIds(['a', 'b']).reps.length, 1);
  c.state.reports.status.b = '2';
  assert.equal(c._cnvSvClusterIds(['a', 'b']).reps.length, 2);
  const v = variant('x', 'unknown', { impact_clinical: { category: 'unknown', reasons: [{ gene: '<script>', impact: '不足' }] } });
  assert.match(c._renderCnvSvImpactReason(v, 'CNV-1A'), /&lt;script&gt;/);
  assert.equal(c._renderCnvSvImpactReason(v, undefined), '');
  v.impact_clinical = { category: 'functional', hpo_score: 42.9,
    reasons: [{ gene: 'SLC25A24', impact: '編碼外顯子' }] };
  assert.match(c._renderCnvSvImpactReason(v, 'CNV-1A'),
    /Exonic \/ splicing · 臨床相關基因：SLC25A24 \(coding exon\) · HPO match 42\.9 \/ 100/);
});


test('actual four-panel rendering includes filters even when all events are hidden', () => {
  const c = fixture(true);
  c.state.data.has_phenotype = true;
  for (const [source, prefix] of [['cnv', 'CNV-1'], ['sv', 'SV-2']]) {
    const id = source + '-n';
    c.state.data[source + '_variants'][id] = variant(id, 'noncoding', { source });
    c.state.data[source + '_categories'][prefix + 'A'] = [id];
    c.state.data[source + '_categories'][prefix + 'B'] = [];
  }
  c.renderCnvSvTabBar();
  assert.match(c.bar.innerHTML, /0 \/ 1/);
  for (const tier of ['CNV-1A', 'SV-2A']) {
    const panel = c.panels[tier];
    assert.match(panel.children[0].innerHTML, /\(1 \/ 1\)/);
    assert.match(panel.children[1].innerHTML, /可勾選其他影響類型/);
    panel.children[0].change({ target: { dataset: { impact: 'noncoding' }, checked: true } });
    const body = panel.children.at(-1);
    assert.equal(body.children.length, 1);
    assert.match(body.children[0].variantId, /-n$/);
  }
});
