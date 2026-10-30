const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const test = require('node:test');
const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
function fn(name) {
  const start = source.indexOf(`function ${name}(`);
  return source.slice(start, source.indexOf('\n}\n', start) + 3);
}
function setup() {
  const context = vm.createContext({});
  const start = source.indexOf('const IN_SILICO_REFERENCES =');
  const end = source.indexOf('function _hasNum', start);
  vm.runInContext(source.slice(start, end) + fn('_hasNum') + fn('_populatedInSilicoTools') + fn('_inSilicoMinimumCount') +
    '\nfunction _evidence(cls, label) { return {cls, label}; }\n' + fn('_toolEvidence') + fn('_gpnStarScoreText') + fn('_gpnStarPairText') + fn('_renderInSilicoCell') + fn('_annotationHint') + `
const escapeHtml = x => String(x);
const escapeAttr = x => String(x);
const fmtNum = x => String(x);`, context);
  return context;
}
test('Star follows MSA and both remain visible, including a zero score', () => {
  const context = setup();
  for (const score of [-9, 0, 2.5]) {
    const tools = context._populatedInSilicoTools({ SpliceAI_score: .1,
      GPN_STAR_score: score, GPN_MSA_score: -10, ESM1b_score: -2 });
    assert.deepEqual(Array.from(tools, t => t.label), ['SpliceAI', 'GPN-MSA', 'GPN-Star-V/M', 'ESM1b']);
    const evidence = context._toolEvidence({ GPN_STAR_score: score }, tools[2]);
    assert.equal(evidence.cls, '');
    assert.equal(evidence.label, 'V100 無分數；M447 不套色');
  }
});
test('legacy and Somatic MSA stay named MSA when Star is absent', () => {
  const context = setup();
  for (const score of [null, '', undefined, NaN]) {
    const tools = context._populatedInSilicoTools({ GPN_STAR_score: score, GPN_MSA_score: -8 });
    assert.deepEqual(Array.from(tools, t => t.label), ['GPN-MSA']);
    assert.equal(context._toolEvidence({ GPN_MSA_score: -8 }, tools[0]).cls, 'sig-lp');
  }
});

test('the More boundary never separates Star from a visible MSA', () => {
  const context = setup();
  const keys = ['pknn', 'alphamissense', 'pangolin', 'revel', 'spliceai', 'gpnmsa', 'gpnstar', 'esm1b'];
  assert.equal(context._inSilicoMinimumCount(keys), 7);
  assert.equal(context._inSilicoMinimumCount(keys.filter(k => k !== 'gpnstar')), 6);
});


test('V/M pair renders one row with only V coloured, fixed order and missing placeholders', () => {
  const context = setup();
  for (const [v, expected, cls] of [
    [{ GPN_STAR_V_score: -7, GPN_STAR_score: -8 }, '-7.0/-8.0', 'sig-vus'],
    [{ GPN_STAR_V_score: 0, GPN_STAR_score: 0 }, '0.0/0.0', 'sig-b'],
    [{ GPN_STAR_V_score: -7.14 }, '-7.1/—', 'sig-vus'],
    [{ GPN_STAR_score: -8.23 }, '—/-8.2', ''],
    [{ GPN_STAR_V_score: NaN, GPN_STAR_score: 2 }, '—/2.0', ''],
  ]) {
    const tools = context._populatedInSilicoTools(v);
    assert.equal(tools.length, 1);
    const html = context._renderInSilicoCell(v, tools[0]);
    assert.equal(context._gpnStarPairText(v), expected);
    assert.equal((html.match(/class="in-silico-row"/g) || []).length, 1);
    assert.ok(html.includes('GPN-Star-V/M'));
    const [vText, mText] = expected.split('/');
    assert.ok(html.includes(`class="gpn-star-v-score ${cls}">${vText}</span>/${mText}`));
    assert.ok(html.includes('class="v ">'));
    assert.ok(html.includes('>紀錄</a>'));
    assert.ok(!html.includes('gpn-star-record-link'));
  }
  assert.equal(context._populatedInSilicoTools({GPN_STAR_V_score: '', GPN_STAR_score: null}).length, 0);
});

test('V100 research colours follow full-precision bin boundaries for every V score', () => {
  const context = setup();
  const tool = Array.from(context._populatedInSilicoTools({ GPN_STAR_V_score: -8 }))[0];
  for (const [score, cls] of [
    [-11.0001, 'sig-p'], [-11, 'sig-p'], [-9.0001, 'sig-p'],
    [-9, 'sig-lp'], [-8.0001, 'sig-lp'], [-8, 'sig-vus'],
    [-6, 'sig-vus'], [-5.0001, 'sig-vus'], [-5, 'sig-b'],
    [-2.0001, 'sig-b'], [-2, 'sig-b'], [0, 'sig-b'],
  ]) {
    const v = { GPN_STAR_V_score: score, GPN_STAR_score: -9, CHROM: 'X', Consequence: 'intron_variant' };
    assert.equal(context._toolEvidence(v, tool).cls, cls);
    assert.ok(context._renderInSilicoCell(v, tool).includes(`class="gpn-star-v-score ${cls}"`));
  }
});

test('research note links to the checked-in Markdown record', () => {
  const context = setup();
  const recordUrl = vm.runInContext('IN_SILICO_CALIBRATIONS.gpnstar.recordUrl', context);
  assert.match(recordUrl, /github\.com\/puddingyd\/NGS-UI\/blob\/.*\/docs\/research\/GPN-Star_.*\.md$/);
  assert.ok(fs.existsSync(path.join(__dirname, '../docs/research', path.basename(recordUrl))));
  const html = context._renderInSilicoCell({ GPN_STAR_V_score: -9 }, Array.from(context._populatedInSilicoTools({ GPN_STAR_V_score: -9 }))[0]);
  assert.ok(html.includes('ClinVar (2026-10-25 版本) ≥2 stars'));
  assert.ok(html.includes('AF <0.01 (gnomAD v4.1 joint)'));
  assert.ok(html.includes('非正式官方公開之標準。'));
  assert.ok(html.includes(`href="${recordUrl}"`));
});

test('tail labels omit the research-bin parenthetical', () => {
  const context = setup();
  const tool = Array.from(context._populatedInSilicoTools({ GPN_STAR_V_score: -12 }))[0];
  assert.equal(context._toolEvidence({ GPN_STAR_V_score: -12 }, tool).label, 'V100: S>11');
  assert.equal(context._toolEvidence({ GPN_STAR_V_score: -1 }, tool).label, 'V100: S≤2');
});
