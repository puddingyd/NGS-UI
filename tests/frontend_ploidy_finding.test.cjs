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

test('analysis and report ploidy controls stay in sync and schedule autosave', async () => {
  const id = 'PLOIDY-chr21-GAIN-test';
  const listeners = {};
  const controls = [];
  let scheduled;
  let delay;
  let saved;
  function control(kind) {
    const classes = new Set([kind]);
    const el = {
      dataset: { id }, value: '',
      classList: {
        add: name => classes.add(name),
        remove: (...names) => names.forEach(name => classes.delete(name)),
        contains: name => classes.has(name),
      },
      matches: selector => selector.split(',').some(s => s.trim() === `.${kind}`),
      closest: selector => selector.includes('.ploidy-finding-card') ? { dataset: { id } } : null,
    };
    controls.push(el);
    return el;
  }
  const analysisAcmg = control('cnv-sv-acmg-select');
  const reportAcmg = control('cnv-sv-acmg-select');
  const analysisDisease = control('cnv-sv-disease-text');
  const reportDisease = control('cnv-sv-disease-text');
  const analysisComment = control('cnv-sv-comment-text');
  const reportComment = control('cnv-sv-comment-text');
  const c = vm.createContext({
    state: { currentLIS: 'LIS-1', dirty: false, reports: { edits: {} } },
    CSS: { escape: value => value },
    document: {
      addEventListener: (event, handler) => { (listeners[event] ||= []).push(handler); },
      querySelectorAll: selector => {
        if (selector === '.js-save-hint') return [];
        const names = [...selector.matchAll(/\.([\w-]+)\[data-id="[^"]+"\]/g)].map(m => m[1]);
        return controls.filter(el => names.some(name => el.classList.contains(name)));
      },
    },
    setTimeout: (callback, ms) => { scheduled = callback; delay = ms; return 1; },
    clearTimeout: () => {},
    renderCnvSvTabBar: () => {},
    renderReportSections: () => {},
    saveChanges: async options => { saved = options; },
  });
  vm.runInContext(source.slice(acmgStart, acmgEnd), c);
  const syncStart = source.indexOf('function getEdit(');
  const syncEnd = source.indexOf('function _syncVariantCheckboxes(', syncStart);
  vm.runInContext(source.slice(syncStart, syncEnd), c);
  const saveStart = source.indexOf('let _autoSaveTimer =');
  const saveEnd = source.indexOf('// Native browser confirmation', saveStart);
  vm.runInContext(source.slice(saveStart, saveEnd), c);
  const hookStart = source.indexOf('// CNV/SV, Ploidy, and Mito edit hooks.');
  const hookEnd = source.indexOf('// Click on a truncated cell', hookStart);
  vm.runInContext(source.slice(hookStart, hookEnd), c);

  analysisAcmg.value = '5';
  listeners.change[0]({ target: analysisAcmg });
  assert.equal(c.state.reports.edits[id].ACMG_class_sv, '5');
  assert.equal(reportAcmg.value, '5');
  assert.ok(analysisAcmg.classList.contains('sig-p'));
  assert.ok(reportAcmg.classList.contains('sig-p'));

  reportAcmg.value = '3';
  listeners.change[0]({ target: reportAcmg });
  assert.equal(analysisAcmg.value, '3');
  assert.ok(analysisAcmg.classList.contains('sig-vus'));
  assert.ok(!analysisAcmg.classList.contains('sig-p'));

  analysisDisease.value = '唐氏症';
  listeners.input[0]({ target: analysisDisease });
  assert.equal(reportDisease.value, '唐氏症');
  reportComment.value = '人工複核';
  listeners.input[0]({ target: reportComment });
  assert.equal(analysisComment.value, '人工複核');
  assert.deepEqual(JSON.parse(JSON.stringify(c.state.reports.edits[id])), {
    ACMG_class_sv: '3', disease: '唐氏症', comment: '人工複核',
  });
  assert.equal(delay, 1500);
  await scheduled();
  assert.equal(saved.silent, true);
});

test('ploidy ACMG dropdown has a distinct color for each class', () => {
  const style = fs.readFileSync(path.join(__dirname, '../frontend/style.css'), 'utf8');
  for (const [value, className] of [[5, 'sig-p'], [4, 'sig-lp'], [3, 'sig-vus'],
    [2, 'sig-lb'], [1, 'sig-b']]) {
    assert.match(style, new RegExp(`\\.ploidy-finding-card \\.cnv-sv-acmg-select\\.${className} \\{ background-color:`));
    assert.match(source, new RegExp(`${value}: "${className}"`));
  }
});
