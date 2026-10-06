const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
const start = source.indexOf('function renderPloidySexStatus(');
const end = source.indexOf('function ploidySexSummary(', start);

test('isolated NCKUH sex suspect is yellow while an autosomal signal remains red', () => {
  const classes = new Set();
  const sexControl = {
    classList: { toggle(name, enabled) {
      if (enabled) classes.add(name);
      else classes.delete(name);
    } },
  };
  const label = { textContent: '', hidden: true };
  const context = vm.createContext({
    state: { data: { ploidy: {}, meta: {} } },
    document: { getElementById(id) {
      return id === 'm-sex-control' ? sexControl : id === 'm-ploidy-call' ? label : null;
    } },
  });
  vm.runInContext(source.slice(start, end), context);

  context.state.data.ploidy = {
    exists: true, karyotype: 'XY', alert_level: 'review', aneuploidy_suspected: true,
    abnormal_chromosomes: [{ chrom: 'chrY' }],
  };
  context.renderPloidySexStatus('M');
  assert.ok(classes.has('ploidy-review'));
  assert.ok(!classes.has('ploidy-aneuploid'));
  assert.match(label.textContent, /1 筆待複核訊號/);

  context.state.data.ploidy.alert_level = 'high';
  context.state.data.ploidy.abnormal_chromosomes.push({ chrom: 'chr21' });
  context.renderPloidySexStatus('M');
  assert.ok(classes.has('ploidy-aneuploid'));
  assert.ok(!classes.has('ploidy-review'));
  assert.match(label.textContent, /2 條染色體劑量訊號/);
});

test('review modal keeps suspect depth visible with an amber status', () => {
  const elements = Object.fromEntries([
    'ploidy-modal', 'ploidy-summary', 'ploidy-alert-table-body',
    'ploidy-all-table-body', 'ploidy-raw-table-body', 'ploidy-alert-wrap',
    'ploidy-alert-empty', 'ploidy-alert-count',
  ].map(id => [id, {
    innerHTML: '', textContent: '', classList: { toggle() {}, remove() {} },
  }]));
  const row = {
    chrom: 'chrY', alt: '.', filter: 'SUSPECT', qual: null, DC: 20,
    NDC: 1.771, observed_ratio: 0.8855, ratio_source: 'native', end: 57227415,
    is_abnormal: true, alert_level: 'review', confidence: 'suspect',
    dosage_call: 'gain', call_label: 'Gain signal',
    interpretation: 'possible sex-chromosome dosage abnormality',
  };
  const context = vm.createContext({
    state: { data: { meta: { Sex: 'M' }, ploidy: {
      exists: true, karyotype: 'XY', alert_level: 'review', pipeline_kind: 'nckuh',
      seq_type: 'WGS', aneuploidy_suspected: true, chromosomes: [row],
      abnormal_chromosomes: [row], qc_warnings: [row],
    } } },
    document: { getElementById: id => elements[id] },
    fmtNum: value => String(value),
    fmtTxt: value => value == null ? '—' : String(value),
    escapeHtml: value => String(value ?? ''),
    escapeAttr: value => String(value ?? ''),
  });
  const modalStart = source.indexOf('function ploidyAlertLevel(');
  const modalEnd = source.indexOf('document.addEventListener("click"', modalStart);
  vm.runInContext(source.slice(modalStart, modalEnd), context);
  context.openPloidyModal();

  assert.match(elements['ploidy-summary'].innerHTML, /ploidy-status-review/);
  assert.match(elements['ploidy-summary'].innerHTML, /ploidy-fact-warn/);
  assert.match(elements['ploidy-alert-table-body'].innerHTML, /ploidy-row-review/);
  assert.match(elements['ploidy-raw-table-body'].innerHTML, /1\.771/);
  assert.match(elements['ploidy-raw-table-body'].innerHTML, /0\.8855/);
});
