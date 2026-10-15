const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
function fn(name) {
  const start = source.indexOf(`function ${name}(`);
  const end = source.indexOf('\nfunction ', start + 1);
  return source.slice(start, end < 0 ? undefined : end);
}
function context() {
  const c = vm.createContext({
    escapeHtml: s => String(s).replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;'),
    _normalizeChrom: c => c,
  });
  // Limit evidence extraction before the following constant declarations.
  const begin = source.indexOf('function renderGbaEvidence(');
  const end = source.indexOf('// GIAB stratification', begin);
  vm.runInContext(source.slice(begin, end) + fn('_igvLocusFor'), c);
  return c;
}
test('GBA evidence keeps null breakpoint unknown and escapes caller strings', () => {
  const c = context();
  const rendered = c.renderGbaEvidence({gba_composition: '<script>x</script>', gba: {
    name: 'RecNciI', total_copy_number: 4, deletion_breakpoint_in_gene: null,
    allele_copy_number: 1, filter: 'PASS', evidence: [], review_required: [],
  }});
  assert.match(rendered, /deletion breakpoint in gene: 未知/);
  assert.match(rendered, /&lt;script&gt;/);
  assert.doesNotMatch(rendered, /<script>/);
  assert.match(rendered, /不等同 GBA1 的 DP／AD／VAF/);
  assert.equal(c.renderGbaEvidence({}), '');
});
test('GBA IGV opens the gene and the paired pseudogene windows', () => {
  const c = context();
  const locus = c._igvLocusFor({CHROM: 'chr1', POS: 155235203, END: 155235252,
    gba: {evidence: [
      {info: {JIDS: 'chr1_155235203_155235203,chr1_155214576_155214576'}},
      {info: {JIDS: 'chr1_155235252_155235252,chr1_155214625_155214625'}},
    ]}});
  assert.equal(locus, 'chr1:155235103-155235352 chr1:155214476-155214725');
});

test('TXT and PDF reuse allele name, second bullet and unknown composite frequency', () => {
  const c = context();
  Object.assign(c, {
    formatVariantTable: () => 'table', ZYG_CH: {},
    pickedDiseaseSummary: () => ({names: 'Gaucher disease', inheritance: '體染色體隱性遺傳'}),
    acmgClassCH: () => '致病性', variantClinSig: () => 'Pathogenic',
    consequenceEntry: () => ({}), clinvarText: () => '', acmgGuidelineText: () => '',
    wrapText: text => [text],
  });
  vm.runInContext(['parseHGVS', 'hgvsCellText', 'gnomadAFText', 'renderVariantBlock', 'pdfWriteVariant'].map(fn).join('\n'), c);
  const v = {HGVS: 'GBA1:NM_000157.4:RecNciI', gene_symbol: 'GBA1',
    gba: {kind: 'haplotype'}, gba_composition: '變異 RecNciI 包含位於同一等位基因上的組成。'};
  assert.equal(c.hgvsCellText(v), 'RecNciI');
  const txt = c.renderVariantBlock('id', v, 'causative').join('\n');
  assert.match(txt, /2\. 變異 RecNciI/);
  assert.match(txt, /3\. 此為致病性/);
  assert.doesNotMatch(txt, /顯示其為罕見/);
  const lines = [];
  c.pdfWriteVariant({subheading: x => lines.push(x), kv: () => {}, gap: () => {}, para: x => lines.push(x)}, 'id', v, 'causative');
  assert.match(lines.join('\n'), /2\. 變異 RecNciI/);
  assert.match(lines.join('\n'), /3\. 此為致病性/);
});
