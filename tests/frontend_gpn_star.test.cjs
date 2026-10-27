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
    '\nfunction _evidence(cls, label) { return {cls, label}; }\n' + fn('_toolEvidence'), context);
  return context;
}
test('Star follows MSA and both remain visible, including a zero score', () => {
  const context = setup();
  for (const score of [-9, 0, 2.5]) {
    const tools = context._populatedInSilicoTools({ SpliceAI_score: .1,
      GPN_STAR_score: score, GPN_MSA_score: -10, ESM1b_score: -2 });
    assert.deepEqual(Array.from(tools, t => t.label), ['SpliceAI', 'GPN-MSA', 'GPN-Star', 'ESM1b']);
    const evidence = context._toolEvidence({ GPN_STAR_score: score }, tools[2]);
    assert.equal(evidence.cls, '');
    assert.equal(evidence.label, 'Contextual score; no calibrated evidence');
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
