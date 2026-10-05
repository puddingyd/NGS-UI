const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');

const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
const start = source.indexOf('function renderPloidyFindingCard(');
const end = source.indexOf('// Reads cnv_variants', start);
const escape = value => String(value ?? '').replace(/[&<>"']/g, character => ({
  '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
})[character]);
const context = vm.createContext({
  document: { createElement: () => ({ dataset: {}, innerHTML: '' }) },
  state: { reports: { status: { 'PLOIDY-chr21-GAIN-test': '1' } } },
  ROH_GRCH38_CHROM_LENGTHS: { chr21: 46709983 },
  _fmtPos: value => Number(value).toLocaleString('en-US'),
  fmtNum: value => String(value),
  escapeHtml: escape,
  escapeAttr: escape,
  _renderStatusRadio: (_id, status) => `<span class="status-radio">${status}</span>`,
  statusOptions: () => ['1', '2', 'C', '0'],
  getEdit: () => '',
});
vm.runInContext(source.slice(start, end), context);

test('chromosome dosage card separates reference extent from measured breakpoints', () => {
  const card = context.renderPloidyFindingCard({
    CHROM: 'chr21', dosage_call: 'gain', interpretation: 'possible trisomy 21',
    NDC: 1.346, filter: 'SUSPECT', pipeline_source: 'NCKUH_PLOIDY_MOSDEPTH',
  }, 'PLOIDY-chr21-GAIN-test');
  assert.match(card.innerHTML, /possible trisomy 21/);
  assert.match(card.innerHTML, /NDC 1\.346 ↑ · SUSPECT/);
  assert.match(card.innerHTML, /GRCh38 參考染色體範圍：chr21:1–46,709,983（非實測斷點）/);
  assert.match(card.innerHTML, /NDC 不代表確定的拷貝數或鑲嵌比例/);
  assert.match(card.innerHTML, /class="status-radio">1</);
});
