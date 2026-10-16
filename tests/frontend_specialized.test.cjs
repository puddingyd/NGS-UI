const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const source = fs.readFileSync(require('node:path').join(__dirname, '../frontend/app.js'), 'utf8');
function fn(name) {
  const start = source.indexOf(`function ${name}(`);
  const end = source.indexOf('\nfunction ', start + 1);
  return source.slice(start, end < 0 ? undefined : end);
}
function context() {
  const c = vm.createContext({
    escapeHtml: s => String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;'),
    escapeAttr: s => String(s).replaceAll('&','&amp;').replaceAll('"','&quot;'),
    getEdit: () => '', _normalizeChrom: x => x,
  });
  vm.runInContext(['renderSpecializedSupport','renderSpecializedEvidence','hgvsCellText','parseHGVS','_igvLocusFor'].map(fn).join('\n'),c);
  return c;
}
test('SMN dosage support retains zero and never uses invented DP or heterozygosity', () => {
  const c = context();
  const v = {HGVS_C:'SMN1 拷貝數減少', specialized:{locus:'smn',kind:'copy_number',smn1_copy_number:0,smn2_copy_number:2,smn2_delta78_copy_number:1}};
  assert.match(c.renderSpecializedSupport(v), /SMN1 CN 0/);
  assert.doesNotMatch(c.renderSpecializedSupport(v), /DP|VAF/);
  assert.equal(c.hgvsCellText(v),'SMN1 拷貝數減少');
  assert.match(c.renderSpecializedEvidence(v,'id'),/SMN2Δ7–8 CN: 1/);
});
test('CYP joint support and ambiguous source are visible and escaped', () => {
  const c = context();
  const v = {specialized:{locus:'cyp21a2',kind:'variant',review_token:'abc',review_required:['<script>'], evidence:[{id:'x',format:{JDP:'80',JAD:'20,60',JAF:'0.75'}}]}};
  assert.match(c.renderSpecializedSupport(v),/JDP 80/);
  assert.match(c.renderSpecializedSupport(v),/不等同目標基因/);
  const html=c.renderSpecializedEvidence(v,'id');
  assert.match(html,/&lt;script&gt;/);
  assert.doesNotMatch(html,/<script>/);
  assert.match(html,/specialized-review/);
  v.specialized.hard_blocked=true;
  assert.doesNotMatch(c.renderSpecializedEvidence(v,'id'),/class="specialized-review"/);
});
test('HBA genotype is shown intact without assigning SEA from double deletion', () => {
  const c=context();
  const v={specialized:{locus:'hba',kind:'locus_genotype',genotype:'--/aa',total_copy_number:2,filter:'PASS',genotype_quality:54}};
  assert.match(c.renderSpecializedEvidence(v,'id'),/--\/aa/);
  assert.doesNotMatch(c.renderSpecializedEvidence(v,'id'),/SEA|Homozygous/);
});
test('IGV targets use only persisted target and homolog coordinates', () => {
  const c=context();
  assert.equal(c._igvLocusFor({specialized:{evidence:[{CHROM:'chr6',POS:32038514,info:{JIDS:'chr6_32005780_32005780'}}]}}),
    'chr6:32038414-32038614 chr6:32005680-32005880');
});
