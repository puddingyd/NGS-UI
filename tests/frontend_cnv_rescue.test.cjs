const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({
  '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
})[c]);
const context = vm.createContext({ escapeHtml: escape });
vm.runInContext(source.slice(source.indexOf('function _renderCnvRescueEvidence('),
  source.indexOf('function _renderCnvSvDetailBox(')), context);
const render = context._renderCnvRescueEvidence;
vm.runInContext(source.slice(source.indexOf('function _dragenJobStepLabel('),
  source.indexOf('function _toggleDragenLog(')), context);
const evidence = {
  rule: 'B', original: { chrom: 'chr1', pos: 1000, end: 8071, filter: 'cnvLength;cnvQual', qual: '5' },
  integrated: { chrom: 'chr1', pos: 1000, end: 7748 },
  sv_support: [{ original_sv_id: '<unsafe>', cnv_overlap: 3404 / 6748, sv_overlap: 1 }],
};

test('rescue is traceable without changing original FILTER or requiring review', () => {
  const html = render({ cnv_rescue: evidence });
  assert.match(html, /Rule B/);
  assert.match(html, /cnvLength;cnvQual/);
  assert.match(html, /50\.44% \/ SV 100\.00%/);
  assert.match(html, /&lt;unsafe&gt;/);
  assert.doesNotMatch(html, /<unsafe>|<button|<input|待檢視/);
});

test('unrescued events render nothing; merged parents show every rescued member', () => {
  assert.equal(render({}), '');
  const html = render({ is_merged_parent: true, cnv_rescue: evidence,
    cnv_rescue_events: [evidence, { ...evidence, rule: 'A', sv_support: [{ original_sv_id: 'full-match' }] }] });
  assert.match(html, /2 個片段/);
  assert.match(html, /Rule A/);
  assert.match(html, /Rule B/);
  assert.match(html, /full-match/);
  assert.equal((html.match(/cnv-rescue-event"/g) || []).length, 2);
});

test('rescue progress follows gene indexing for each sample in a batch', () => {
  assert.equal(context._dragenJobStepLabel({ step: 'post-processing:cnv-rescue' }), 'DRAGEN CNV rescue');
  for (const index of [0, 1]) {
    const state = { state: 'running', post_processing_sample_count: 2, post_processing_sample_index: index };
    const before = context._dragenProgressPercent({ ...state, step: 'sample-step:gene-index' });
    const rescue = context._dragenProgressPercent({ ...state, step: 'post-processing:cnv-rescue' });
    assert.ok(rescue >= before);
    assert.ok(rescue < 100);
  }
});
