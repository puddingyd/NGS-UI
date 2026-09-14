const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const test = require('node:test');

const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
const escape = value => String(value).replace(/[&<>"']/g, character => ({
  '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
})[character]);
const context = vm.createContext({ escapeHtml: escape, escapeAttr: escape });
vm.runInContext(source.slice(source.indexOf('const CLINVAR_ABBREV ='),
  source.indexOf('function vusSubclassForScore')), context);
const { formatClinvar, renderClinvarValue } = context;
const visibleText = html => html.replace(/<[^>]*>/g, '');
const tokens = html => [...html.matchAll(/class="clinvar-token ([^"]*)">([^<]*)/g)]
  .map(match => [match[1], match[2]]);

test('pipeline underscores before counts are removed; counts and review stars survive', () => {
  const sig = 'Pathogenic_(22)|Uncertain_significance_(1)';
  assert.equal(formatClinvar(sig, '', 1), 'P(22)|VUS(1)(1★)');
  const html = renderClinvarValue(sig, '', 1);
  assert.equal(visibleText(html), 'P(22)|VUS(1)(1★)');
  assert.deepEqual(tokens(html), [['sig-p', 'P(22)'], ['sig-vus', 'VUS(1)']]);
  assert.match(html, /title="Pathogenic \(22\)\|Uncertain significance \(1\) \(1★\)"/);
});

test('case, spaces, single counted assertions, and combined classifications normalize', () => {
  for (const [raw, expected] of [
    [' likely_PATHOGENIC_ ( 12 ) ', 'LP(12)'],
    ['Benign_(4)', 'B(4)'],
    ['Likely benign(2)', 'LB(2)'],
    ['Variant of uncertain significance', 'VUS'],
    ['Pathogenic/Likely_pathogenic', 'P/LP'],
    ['Likely_benign/Benign', 'LB/B'],
    ['P(22)|VUS(1)', 'P(22)|VUS(1)'],
  ]) assert.equal(formatClinvar(raw, '', null), expected);
});

test('conflicting calls use assertion counts from CLNSIGCONF, including zero-star reviews', () => {
  const sig = 'conflicting_classifications_of_pathogenicity';
  const conf = 'Benign_(12),_Likely_benign_(3),_Uncertain_significance_(1)';
  const html = renderClinvarValue(sig, conf, 0);
  assert.equal(visibleText(html), 'B(12)|LB(3)|VUS(1)(0★)');
  assert.deepEqual(tokens(html), [['sig-b', 'B(12)'], ['sig-lb', 'LB(3)'], ['sig-vus', 'VUS(1)']]);
  assert.match(html, /conflicting classifications of pathogenicity/);
  assert.equal(formatClinvar('Conflicting_interpretations_of_pathogenicity', conf, 1),
    'B(12)|LB(3)|VUS(1)(1★)');
});

test('opposing P and B assertions remain distinct; VUS is yellow regardless of order', () => {
  for (const sig of [
    'Pathogenic_(22)|Benign_(1)|Uncertain_significance_(2)',
    'Uncertain_significance_(2)&Benign_(1)&Pathogenic_(22)',
  ]) {
    const values = tokens(renderClinvarValue(sig, '', 2));
    assert.deepEqual(values.find(item => item[1] === 'P(22)'), ['sig-p', 'P(22)']);
    assert.deepEqual(values.find(item => item[1] === 'B(1)'), ['sig-b', 'B(1)']);
    assert.deepEqual(values.find(item => item[1] === 'VUS(2)'), ['sig-vus', 'VUS(2)']);
  }
  assert.deepEqual(tokens(renderClinvarValue('Pathogenic/Benign', '', null)),
    [['sig-p', 'P'], ['sig-b', 'B']]);
});

test('missing values never acquire a classification or a standalone star rating', () => {
  for (const value of [null, undefined, '', '.', 'NA', 'N/A', 'na']) {
    assert.equal(formatClinvar(value, 'Pathogenic_(1)', 0), '—');
    assert.equal(renderClinvarValue(value, 'Pathogenic_(1)', 0), '—');
  }
  for (const conf of ['', '.', 'NA', 'N/A']) {
    assert.equal(formatClinvar('Conflicting_classifications_of_pathogenicity', conf, 1), 'Conflict(1★)');
  }
});

test('nonstandard clinical terms remain visible and uncolored, with escaped source text', () => {
  assert.equal(formatClinvar('risk_factor_(3)|drug_response_(1)', '', null), 'risk factor(3)|drug response(1)');
  assert.deepEqual(tokens(renderClinvarValue('Pathogenic,_low_penetrance|risk_factor_(3)', '', null)),
    [['sig-p', 'P'], ['', 'low penetrance'], ['', 'risk factor(3)']]);
  const html = renderClinvarValue('<img src=x onerror="alert(1)">_(3)', '', null);
  assert.ok(!html.includes('<img'));
  assert.ok(html.includes('&lt;img'));
});
