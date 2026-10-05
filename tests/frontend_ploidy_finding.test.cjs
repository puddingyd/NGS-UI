const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');

const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
const start = source.indexOf('function renderPloidyFindingCard(');
const end = source.indexOf('// Reads cnv_variants', start);
const acmgStart = source.indexOf('const SV_ACMG_LABELS =');
const acmgEnd = source.indexOf('function _fmtPos(', acmgStart);
const classStart = source.indexOf('function _cnvSvAcmgClassValue(');
const classEnd = source.indexOf('function _renderCnvSvHeader(', classStart);
const editorStart = source.indexOf('function _renderCnvSvComment(');
const editorEnd = source.indexOf('function renderCnvSvCard(', editorStart);
const escape = value => String(value ?? '').replace(/[&<>"']/g, character => ({
  '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
})[character]);
const context = vm.createContext({
  document: { createElement: () => ({ dataset: {}, innerHTML: '' }) },
  state: { reports: { status: { 'PLOIDY-chr21-GAIN-test': '1' }, edits: {
    'PLOIDY-chr21-GAIN-test': { disease: '唐氏症', comment: '人工複核', ACMG_class_sv: '5' },
  } } },
  ROH_GRCH38_CHROM_LENGTHS: { chr21: 46709983 },
  _fmtPos: value => Number(value).toLocaleString('en-US'),
  fmtNum: value => String(value),
  escapeHtml: escape,
  escapeAttr: escape,
  _renderStatusRadio: (_id, status) => `<span class="status-radio">${status}</span>`,
  statusOptions: () => ['1', '2', 'C', '0'],
  getEdit: (id, field) => context.state.reports.edits[id]?.[field] || '',
});
vm.runInContext(source.slice(acmgStart, acmgEnd), context);
vm.runInContext(source.slice(classStart, classEnd), context);
vm.runInContext(source.slice(editorStart, editorEnd), context);
vm.runInContext(source.slice(start, end), context);

test('ploidy card uses CNV-style title and manual ACMG and disease fields', () => {
  const card = context.renderPloidyFindingCard({
    CHROM: 'chr21', dosage_call: 'gain', interpretation: 'possible trisomy 21',
    NDC: 1.346, filter: 'SUSPECT', pipeline_source: 'NCKUH_PLOIDY_MOSDEPTH',
  }, 'PLOIDY-chr21-GAIN-test');
  assert.match(card.innerHTML, /Ploidy VCF<\/span>\s*<span class="sv-type-pill sv-type-ploidy">Trisomy/);
  assert.match(card.innerHTML, /NDC 1\.346 ↑ · SUSPECT/);
  assert.match(card.innerHTML, /chr21（GRCh38 參考染色體範圍：chr21:1–46,709,983）/);
  assert.match(card.innerHTML, /class="cnv-sv-acmg-select sig-p"[^>]*>\s*<option/);
  assert.match(card.innerHTML, /option value="5" selected>Pathogenic/);
  assert.match(card.innerHTML, /cnv-sv-disease-text[^>]*>唐氏症/);
  assert.match(card.innerHTML, /cnv-sv-comment-text[^>]*>人工複核/);
  assert.ok(card.innerHTML.indexOf('cnv-sv-disease-text') < card.innerHTML.indexOf('cnv-sv-comment-text'));
  assert.doesNotMatch(card.innerHTML, /NDC 不代表確定的拷貝數/);
  assert.match(card.innerHTML, /class="status-radio">1</);
});
