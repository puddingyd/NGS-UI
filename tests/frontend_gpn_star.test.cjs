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
  const start = source.indexOf('const IN_SILICO_TOOLS =');
  const end = source.indexOf('function _hasNum', start);
  vm.runInContext(source.slice(start, end) + fn('_hasNum') + fn('_populatedInSilicoTools') + fn('_inSilicoMinimumCount') +
    '\nfunction _evidence(cls, label) { return {cls, label}; }\n' + fn('_toolEvidence') + fn('_gpnStarPairText') + fn('_renderInSilicoCell') + fn('_annotationHint') + `
const IN_SILICO_CALIBRATIONS = {};
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
    assert.equal(evidence.label, 'Research only; no calibrated evidence');
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


test('V/M pair renders one neutral row in fixed order with one decimal and missing placeholders', () => {
  const context = setup();
  for (const [v, expected] of [
    [{ GPN_STAR_V_score: -7, GPN_STAR_score: -8 }, '-7.0/-8.0'],
    [{ GPN_STAR_V_score: 0, GPN_STAR_score: 0 }, '0.0/0.0'],
    [{ GPN_STAR_V_score: -7.14 }, '-7.1/—'],
    [{ GPN_STAR_score: -8.23 }, '—/-8.2'],
    [{ GPN_STAR_V_score: NaN, GPN_STAR_score: 2 }, '—/2.0'],
  ]) {
    const tools = context._populatedInSilicoTools(v);
    assert.equal(tools.length, 1);
    const html = context._renderInSilicoCell(v, tools[0]);
    assert.equal(context._gpnStarPairText(v), expected);
    assert.equal((html.match(/class="in-silico-row"/g) || []).length, 1);
    assert.ok(html.includes('GPN-Star-V/M'));
    assert.ok(html.includes(`>${expected}</span>`));
    assert.ok(html.includes('Research only'));
    assert.ok(!html.includes('sig-'));
  }
  assert.equal(context._populatedInSilicoTools({GPN_STAR_V_score: '', GPN_STAR_score: null}).length, 0);
});
