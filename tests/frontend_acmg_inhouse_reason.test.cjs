const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const test = require('node:test');

const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
const escape = value => String(value).replace(/[&<>"']/g, c => ({
  '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
})[c]);
const start = source.indexOf('const ACMG_INHOUSE_STATUS =');
const end = source.indexOf('function renderAcmgCriteriaEditor()');
const context = vm.createContext({ escapeHtml: escape, escapeAttr: escape, _acmgEditor: null,
  _acmgCatalog: { strengths: [{ value: 'strong', label: 'Strong' }, { value: 'supporting', label: 'Supporting' }] } });
vm.runInContext(source.slice(start, end) + '\nthis.setEditor = e => { _acmgEditor = e; };', context);

test('In-house reasons render under each criterion, met first, with strength and points', () => {
  context.setEditor({ inhouseEvidence: {
    status: 'stale', version: 'ngsui-acmg-v1.9', classification: 'Likely benign', points: -2,
    rank_classification: 'Uncertain significance', clinvar_action: 'adopt', clinvar_source: 'ClinVar latest',
    clinvar_category: 'VUS', clinvar_stars: 2, vcep: { vcep: 'FBN1 VCEP', svi: 'GN022' },
    evidence: [
      { code: 'BS1', status: 'met', strength: 'strong', points: -4, reason: '最大 AF 0.00017 > 5e-05' },
      { code: 'PM1', status: 'not_met', strength: null, points: 0, reason: '不在 <熱點>' },
      { code: 'PS4', status: 'not_evaluated', strength: null, points: 0, reason: '' },
    ] } });
  const bs1 = context._acmgInhouseReasonHtml('BS1');
  assert.match(bs1, /acmg-inhouse-met/);
  assert.match(bs1, /In-house：符合 Strong（-4）<\/b>：最大 AF 0\.00017/);
  assert.match(context._acmgInhouseReasonHtml('PM1'), /不符合<\/b>：不在 &lt;熱點&gt;/);
  assert.match(context._acmgInhouseReasonHtml('PS4'), /未評估<\/b><\/div>$/);
  assert.equal(context._acmgInhouseReasonHtml('PP5'), '');
  const header = context._acmgInhouseHeaderHtml();
  assert.doesNotMatch(header, /ngsui-acmg|下方各項/);
  assert.match(header, /舊版規則/);
  assert.match(header, /排序改用 ClinVar latest 的 VUS（2★）/);
  context.setEditor({ inhouseEvidence: { status: 'current', evidence: [], clinvar_action: 'overridden',
    flags: [{ code: 'clinvar_overridden', message: 'ClinVar P（2★）但最大 AF 0.34 > 0.05' }] } });
  assert.match(context._acmgInhouseHeaderHtml(), /ClinVar P（2★）但最大 AF 0\.34/);
  context.setEditor({ inhouseEvidence: { status: 'current', evidence: [] } });
  assert.equal(context._acmgInhouseHeaderHtml(), '');                 // nothing to say → no header
  context.setEditor({ inhouseEvidence: null });
  assert.equal(context._acmgInhouseReasonHtml('BS1'), '');
  assert.equal(context._acmgInhouseHeaderHtml(), '');
});
