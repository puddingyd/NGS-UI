const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');

const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');

function context(reportDiseases = {}) {
  const c = vm.createContext({
    state: { reports: { edits: { v1: { report_diseases: reportDiseases } } } },
  });
  const start = source.indexOf('function pickedDiseaseSlots(');
  const end = source.indexOf('// HGVS', start);
  vm.runInContext(`
    const OMIM_DISEASE_SLOT_COUNT = 16;
    const INHERITANCE_LABELS = {
      AD: '體染色體顯性遺傳', AR: '體染色體隱性遺傳',
      XLD: '性染色體顯性遺傳', XLR: '性染色體隱性遺傳',
      XL: '性聯遺傳', YL: 'Y 染色體遺傳', MT: '粒線體遺傳',
      DD: '雙等位基因顯性遺傳', IC: '細胞質遺傳',
    };
    ${source.slice(start, end)}
  `, c);
  return c;
}

test('SNV report summary keeps every checked disease and phenotype MIM in slot order', () => {
  const c = context({ 2: true, 1: true });
  const summary = c.pickedDiseaseSummary('v1', {
    Disease1: 'Disease A (600001)(AD)',
    Disease2: 'Disease B (600001)(AR)',
    Disease3: 'Disease C (600003)(XLR)',
  });

  assert.equal(summary.names, 'Disease A、Disease B');
  assert.equal(summary.inheritance, '體染色體顯性遺傳、體染色體隱性遺傳');
  assert.equal(summary.phenotypeMims, '600001、600001');
});

test('SNV report summary falls back to the first disease when none is checked', () => {
  const c = context();
  const summary = c.pickedDiseaseSummary('v1', {
    Disease1: 'Disease A (600001)(AD)',
    Disease2: 'Disease B (600002)(AR)',
  });

  assert.equal(summary.names, 'Disease A');
  assert.equal(summary.phenotypeMims, '600001');
});
