const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');

const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');

function context(reportDiseases = {}, extraEdits = {}) {
  const c = vm.createContext({
    state: { reports: { edits: { v1: { report_diseases: reportDiseases, ...extraEdits } } } },
  });
  const start = source.indexOf('function pickedDiseaseSlots(');
  const end = source.indexOf('// HGVS', start);
  const labelsStart = source.indexOf('const INHERITANCE_LABELS =');
  const labelsEnd = source.indexOf('const ZYG_CH', labelsStart);
  vm.runInContext(`
    const OMIM_DISEASE_SLOT_COUNT = 16;
    ${source.slice(labelsStart, labelsEnd)}
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

test('SNV report summary describes SMu as a somatic variant', () => {
  const c = context({ 1: true });
  const summary = c.pickedDiseaseSummary('v1', {
    Disease1: '?RAS-associated autoimmune lymphoproliferative syndrome type 1, somatic (614470)(SMu)\n\nCurator detail.',
  });

  assert.equal(summary.names, '?RAS-associated autoimmune lymphoproliferative syndrome type 1, somatic');
  assert.equal(summary.inheritance, '體細胞變異');
  assert.equal(summary.phenotypeMims, '614470');
});

test('content-keyed picks survive a reordered OMIM workbook', () => {
  const c = context({ 3: true }, {
    report_disease_picks: {
      '613986|pituitary hormone deficiency combined 6': {
        label: 'Pituitary hormone deficiency, combined, 6 (613986)(AD)', slot: 3,
      },
    },
  });
  const summary = c.pickedDiseaseSummary('v1', {
    Disease1: 'Pituitary hormone deficiency, combined, 6 (613986)(AD)',
    Disease2: 'Microphthalmia, syndromic 5 (610125)(AD)',
    Disease3: 'Retinal dystrophy, early-onset (610125)(AD)',
  });
  assert.equal(summary.names, 'Pituitary hormone deficiency, combined, 6');
  assert.equal(summary.phenotypeMims, '613986');
});

test('a ticked disease missing from the workbook keeps its snapshot label', () => {
  const c = context({}, {
    report_disease_picks: {
      '600099|gone disease': { label: 'Gone disease (600099)(AR)', slot: 2 },
    },
  });
  const slots = c.pickedDiseaseSlots('v1', { Disease1: 'Other disease (600001)(AD)' });
  assert.equal(slots.length, 1);
  assert.equal(slots[0].idx, null);
  assert.equal(slots[0].text, 'Gone disease (600099)(AR)');
});

test('same phenotype MIM with different names keeps the ticked label', () => {
  const c = context({}, {
    report_disease_picks: {
      '610125|retinal dystrophy early onset': { label: 'Retinal dystrophy, early-onset (610125)(AD)', slot: 2 },
    },
  });
  const slots = c.pickedDiseaseSlots('v1', {
    Disease1: 'Microphthalmia, syndromic 5 (610125)(AD)',
    Disease2: 'Retinal dystrophy, early-onset (610125)(AD)',
  });
  assert.equal(JSON.stringify(slots.map(s => s.idx)), '[2]');
});
