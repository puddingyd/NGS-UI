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
  const elements = new Map();
  function element(id) {
    if (!elements.has(id)) elements.set(id, { checked: false, hidden: false,
      classList: { toggle(_name, value) { elements.get(id).hidden = value; } } });
    return elements.get(id);
  }
  const context = vm.createContext({ document: { getElementById: element },
    state: { currentLIS: 'S1', data: { somatic: { completed: false } } },
    somaticUi: { sid: 'S1' }, _passesNckuhCommonFilter: () => true,
    _numericValue: x => x, _isReferenceZygosity: () => false, _isClinvarPlp: () => false });
  vm.runInContext('function somaticSampleId() { return state.data?.sample_id || state.currentLIS; }\n' +
    fn('_passesMainSnvDisplayFilters') + fn('renderSomaticControls'), context);
  return { context, element };
}
test('Somatic checkbox is absent before a completed run, including zero-result completion', () => {
  const { context, element } = setup();
  context.renderSomaticControls();
  assert.equal(element('filter-somatic-label').hidden, true);
  context.state.data.somatic.completed = true;
  context.renderSomaticControls();
  assert.equal(element('filter-somatic-label').hidden, false);
});
test('Somatic toggle only changes additional calls, not germline low-VAF behavior', () => {
  const { context, element } = setup();
  const germline = { alt_af: 0.03 };
  const original = JSON.stringify(germline);
  element('filter-somatic').checked = true;
  assert.equal(context._passesMainSnvDisplayFilters({ ...germline, somatic: true }), true);
  assert.equal(context._passesMainSnvDisplayFilters(germline), false);
  element('filter-somatic').checked = false;
  element('filter-vaf').checked = true;
  assert.equal(context._passesMainSnvDisplayFilters({ ...germline, somatic: true }), false);
  assert.equal(context._passesMainSnvDisplayFilters(germline), true);
  assert.equal(JSON.stringify(germline), original);
});
test('explicitly targeted somatic variants survive germline gene-scope filters', () => {
  const { context, element } = setup();
  element('filter-somatic').checked = true;
  element('filter-disease-associated').checked = true;
  element('filter-in-panel-only').checked = true;
  assert.equal(context._passesMainSnvDisplayFilters({ somatic: true, in_panel: false }), true);
  assert.equal(context._passesMainSnvDisplayFilters({ in_panel: false }), false);
});

test('Somatic progress uses the same stable percentage model as the tertiary panel', () => {
  const context = vm.createContext({ SOMATIC_PROGRESS: { queued: 1, mutect2: 15,
    'annotation:vep': 72, publishing: 98, completed: 100 } });
  vm.runInContext(fn('somaticProgressPercent'), context);
  assert.equal(context.somaticProgressPercent({ status: 'running', step: 'mutect2' }), 15);
  assert.equal(context.somaticProgressPercent({ status: 'running', step: 'annotation:vep' }), 72);
  assert.equal(context.somaticProgressPercent({ status: 'completed', step: 'completed' }), 100);
  assert.equal(context.somaticProgressPercent({ status: 'failed', step: 'publishing' }), 98);
});

test('Somatic history uses deletion and card quality metrics are always expanded', () => {
  const html = fs.readFileSync(path.join(__dirname, '../frontend/index.html'), 'utf8');
  assert.match(source, /class="btn btn-danger somatic-delete"/);
  assert.doesNotMatch(source, /somatic-rerun/);
  assert.match(source, /class="somatic-qc-metrics">品質：/);
  assert.doesNotMatch(source, /<details><summary>品質資訊<\/summary>/);
  assert.match(html, /<span>VAF &lt; 0\.2<\/span>/);
  assert.doesNotMatch(html, /VAF &lt; 0\.2 \/ zygosity=ref/);
});

test('Somatic modal uses concise coordinate labels and no standard-only predictor notice', () => {
  const html = fs.readFileSync(path.join(__dirname, '../frontend/index.html'), 'utf8');
  assert.match(html, /Genomic position（GRCh38）/);
  assert.doesNotMatch(html, /Genomic position（GRCh38，1-based）/);
  assert.doesNotMatch(html, /基因和座標取聯集。僅新增 germline 未有的點位/);
  assert.match(source, /<summary>座標<\/summary>/);
  assert.doesNotMatch(source, /<summary>座標（1-based）<\/summary>/);
  assert.doesNotMatch(source, /標準 Somatic run 不會填入/);
});

test('Somatic candidates use a separate card modal and refresh is automatic', () => {
  const html = fs.readFileSync(path.join(__dirname, '../frontend/index.html'), 'utf8');
  assert.match(html, /id="somatic-candidates-modal"/);
  assert.match(source, /renderSomaticCandidates/);
  assert.match(source, /renderVariantCard\(entry\.variant, entry\.id/);
  assert.match(source, /顯示原始 Log/);
  assert.doesNotMatch(html, /id="somatic-refresh-btn"/);
  assert.doesNotMatch(source, /button\.id === "somatic-refresh-btn"/);
});
