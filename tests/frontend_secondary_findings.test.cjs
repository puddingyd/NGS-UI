const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const test = require('node:test');

const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
const definitions = source.slice(source.indexOf('const SECONDARY_PANEL_DEFS = ['), source.indexOf('const REPORT_SECTION_DEFS = ['));
const selection = source.slice(source.indexOf('function _isSecondaryEligible(id)'), source.indexOf('function getEdit(id, field)'));
const panels = ['acmg_sf', 'hereditary_cancer', 'stroke', 'carrier'];

function setup() {
  const state = {
    data: {
      variants: {
        shared: { tier: '1C', CLNSIG: 'Uncertain_significance' },
        pathogenic: { tier: '1A', CLNSIG: 'Pathogenic' },
      },
      categories: Object.fromEntries(panels.map(key => [key, ['shared', 'pathogenic']])),
    },
    reports: { secondary_findings: {} },
  };
  const context = vm.createContext({
    state,
    _isClinvarPlp: variant => ['Pathogenic', 'Likely_pathogenic'].includes(variant?.CLNSIG),
    renderReportSections() {},
    renderCandidateSections() {},
    _syncStatusRadios() {},
    updateSaveHint() {},
  });
  vm.runInContext(definitions + selection, context);
  return { state, context };
}

test('all four panels share ClinVar-only default selection', () => {
  const { context } = setup();
  assert.deepEqual(Array.from(vm.runInContext('SECONDARY_PANEL_DEFS.map(def => def.key)', context)), panels);
  for (const panel of panels) {
    assert.equal(context.isSecondarySelected('pathogenic', panel), true);
    assert.equal(context.isSecondarySelected('shared', panel), false);
  }
  assert.equal(context._isSecondaryEligible('shared'), true);
});

test('selecting or dismissing cancer findings synchronizes every matching panel', () => {
  const { state, context } = setup();
  context.setPanelStatus('shared', 'hereditary_cancer', true);
  for (const panel of panels) {
    assert.equal(context.getPanelStatus('shared', panel), '✓');
    assert.deepEqual(Array.from(state.reports.secondary_findings[panel].selected), ['shared']);
  }
  context.setPanelStatus('shared', 'acmg_sf', false);
  for (const panel of panels) {
    assert.equal(context.getPanelStatus('shared', panel), '');
    assert.deepEqual(Array.from(state.reports.secondary_findings[panel].dismissed), ['shared']);
  }
  assert.equal(state.dirty, true);
});

test('an existing cancer-panel dismissal overrides selection until explicitly selected again', () => {
  const { state, context } = setup();
  state.reports.secondary_findings.acmg_sf = { selected: ['pathogenic'] };
  state.reports.secondary_findings.hereditary_cancer = { dismissed: ['pathogenic'] };
  for (const panel of panels) assert.equal(context.isSecondarySelected('pathogenic', panel), false);
  context.setPanelStatus('pathogenic', 'hereditary_cancer', true);
  for (const panel of panels) assert.equal(context.isSecondarySelected('pathogenic', panel), true);
});

test('a cancer-only variant is not added to another panel by selection', () => {
  const { state, context } = setup();
  state.data.variants.cancerOnly = { tier: '1C', CLNSIG: 'Uncertain_significance' };
  state.data.categories.hereditary_cancer.push('cancerOnly');
  context.setPanelStatus('cancerOnly', 'hereditary_cancer', true);
  assert.deepEqual(Object.keys(state.reports.secondary_findings), ['hereditary_cancer']);
  assert.equal(context.isSecondarySelected('cancerOnly', 'hereditary_cancer'), true);
});
