const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');

const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');

function context(edits = {}) {
  const c = vm.createContext({
    state: { reports: { edits: { cnv1: edits } } },
    getEdit(id, field) { return id === 'cnv1' ? edits[field] : undefined; },
    escapeHtml(value) {
      return String(value ?? '').replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;');
    },
    escapeAttr(value) {
      return String(value ?? '').replaceAll('&', '&amp;').replaceAll('"', '&quot;').replaceAll('<', '&lt;');
    },
    diseaseAssociationSummary(a) {
      return `${a.display_name || ''}${a.phenotype_mim ? ` (${a.phenotype_mim})` : ''}${a.inheritance ? `(${a.inheritance})` : ''}`;
    },
    diseaseAssociationDetail(a) { return a.detail || a.display_name || ''; },
    diseaseSourceBadges() { return ''; },
  });
  const start = source.indexOf('function _cnvReportDiseaseItems');
  const end = source.indexOf('function _renderCnvSvBenign', start);
  vm.runInContext(source.slice(start, end), c);
  return c;
}

test('CNV phenotype uses selectable OMIM workbook disease rows', () => {
  const key = 'omim:COL1A1:120150:1';
  const c = context({ report_disease_items: { [key]: { label: 'Caffey disease' } } });
  const html = c._renderCnvGeneDiseases({
    gene: 'COL1A1', omim_id: '120150',
    disease_associations: [{
      id: 'omim-slot:1', omim_slot: 1, source_kind: 'omim',
      display_name: 'Caffey disease', phenotype_mim: '114000', inheritance: 'AD',
      detail: 'Caffey disease detail',
    }],
  }, 'cnv1');
  assert.match(html, /cnv-report-disease-pick/);
  assert.match(html, /Caffey disease/);
  assert.match(html, /checked/);
  assert.match(html, /data-phenotype-mim="114000"/);
});

test('CNV phenotype shows two OMIM diseases before the count-based expander', () => {
  const c = context();
  const disease_associations = Array.from({ length: 4 }, (_, idx) => ({
    id: `omim-slot:${idx + 1}`, omim_slot: idx + 1, source_kind: 'omim',
    display_name: `Disease ${idx + 1}`, phenotype_mim: `60000${idx + 1}`,
    inheritance: 'AD', detail: `Detail ${idx + 1}`,
  }));
  const html = c._renderCnvGeneDiseases({
    gene: 'GENE1', omim_id: '123456', disease_associations,
  }, 'cnv1');
  assert.equal((html.match(/cnv-gene-disease cnv-disease-extra/g) || []).length, 2);
  assert.match(html, /展開其餘 2 個疾病/);
  assert.match(html, /aria-expanded="false"/);
});

test('pathogenic overlap renders one report checkbox per disease', () => {
  const c = context();
  const html = c._renderCnvSvOverlap({
    sv_type: 'DEL',
    p_loss: { diseases: ['Disease A', 'Disease B'], sources: ['CLN:1'], coords: [] },
  }, 'cnv1');
  assert.equal((html.match(/cnv-report-disease-pick/g) || []).length, 2);
  assert.match(html, /Disease A/);
  assert.match(html, /Disease B/);
  assert.doesNotMatch(html, /p_gain/);
});

test('pathogenic overlap shows five diseases before the count-based expander', () => {
  const c = context();
  const diseases = Array.from({ length: 7 }, (_, idx) => `Disease ${idx + 1}`);
  const html = c._renderCnvSvOverlap({
    sv_type: 'DEL', p_loss: { diseases, sources: ['CLN:1'], coords: [] },
  }, 'cnv1');
  assert.equal((html.match(/cnv-overlap-disease-option cnv-disease-extra/g) || []).length, 2);
  assert.match(html, /展開其餘 2 個疾病/);
  assert.match(html, /cnv-sv-overlap-content-counted/);
  assert.doesNotMatch(html, /cnv-sv-overlap-toggle/);
});
