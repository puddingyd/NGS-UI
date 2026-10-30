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
    [-11.6001, 'sig-p'], [-11.6, 'sig-p'], [-11.5999, 'sig-p'],
    [-9.8, 'sig-p'], [-9.7999, 'sig-lp'], [-8.4, 'sig-lp'],
    [-8.3999, 'sig-vus'], [-5.1001, 'sig-vus'], [-5.1, 'sig-lb'],
    [-3.8001, 'sig-lb'], [-3.8, 'sig-b'], [-1.7001, 'sig-b'],
    [-1.7, 'sig-b'], [0, 'sig-b'],
  ]) {
    const v = { GPN_STAR_V_score: score, GPN_STAR_score: -9, CHROM: 'X', Consequence: 'intron_variant' };
    assert.equal(context._toolEvidence(v, tool).cls, cls);
    assert.ok(context._renderInSilicoCell(v, tool).includes(`class="gpn-star-v-score ${cls}"`));
  }
});

test('research note links to the checked-in Markdown record', () => {
  const context = setup();
  const recordUrl = vm.runInContext('IN_SILICO_CALIBRATIONS.gpnstar.recordUrl', context);
  assert.match(recordUrl, /github\.com\/puddingyd\/NGS-UI\/blob\/.*\/docs\/research\/GPN-Star_V100_PP3_BP4_Final_Methods_Results_Conclusion_20261030\.md$/);
  assert.ok(fs.existsSync(path.join(__dirname, '../docs/research', path.basename(recordUrl))));
  const html = context._renderInSilicoCell({ GPN_STAR_V_score: -9 }, Array.from(context._populatedInSilicoTools({ GPN_STAR_V_score: -9 }))[0]);
  assert.ok(html.includes('ClinVar (2026-10-25 版本) ≥2 stars'));
  assert.ok(html.includes('AF <0.01 (gnomAD v4.1 joint)'));
  assert.ok(html.includes('非正式官方公開之標準。'));
  assert.ok(html.includes(`href="${recordUrl}"`));
  assert.ok(!html.includes('本位點：'));
  for (const line of [
    'Threshold (V model)',
    '- BP4_Strong：≥ -1.7 Strong（sparse data; need further validation）',
    '- BP4_Moderate：-1.7 to -3.8',
    '- BP4_Supporting：-3.8 to -5.1',
    '- Indeterminate：-5.1 to -8.4',
    '- PP3_Supporting：-8.4 to -9.8',
    '- PP3_Moderate：-9.8 to -11.6',
    '- PP3_Strong：≤ -11.6 Strong（sparse data; need further validation）',
  ]) {
    assert.ok(html.includes(line));
  }
});

test('tail labels omit the research-bin parenthetical', () => {
  const context = setup();
  const tool = Array.from(context._populatedInSilicoTools({ GPN_STAR_V_score: -12 }))[0];
  for (const [score, expected] of [
    [-11.6, 'V100: S≥11.6'],
    [-11.5999, 'V100: PP3 Moderate 候選（9.8≤S<11.6）'],
    [-1.7001, 'V100: BP4 Moderate 候選（1.7<S≤3.8）'],
    [-1.7, 'V100: S≤1.7'],
  ]) {
    assert.equal(context._toolEvidence({ GPN_STAR_V_score: score }, tool).label, expected);
  }
});
