const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const source = fs.readFileSync(require('node:path').join(__dirname, '../frontend/app.js'), 'utf8');

test('tertiary terminal state keeps polling PRS and refreshes index only once', async () => {
  const code = source.slice(source.indexOf('function _dragenStartPolling('), source.indexOf('async function _dragenRecoverActiveJob('));
  let tick, cleared = 0, refreshed = 0, seen = [];
  let state = {state:'done', prs_active:true};
  const context = {
    _DRAGEN_STATE:{pollTimer:null}, apiFetch:async () => state,
    _dragenSetJob:s => seen.push(s), loadIndex:async () => refreshed++,
    setInterval:f => {tick=f; return 1;}, clearInterval:() => cleared++,
  };
  vm.runInNewContext(code, context);
  context._dragenStartPolling('parent-1');
  await new Promise(resolve => setImmediate(resolve));
  assert.equal(cleared, 0);
  assert.equal(refreshed, 1);
  state = {state:'done', prs_active:false};
  await tick();
  assert.equal(cleared, 1);
  assert.equal(refreshed, 1);
  assert.equal(seen.length, 2);
});

test('selected PRS reserves postprocessing progress and reaches 100 only at parent done',()=>{
  const code=source.slice(source.indexOf('function _dragenJobStepLabel('),source.indexOf('function _toggleDragenLog('));
  const context={};vm.runInNewContext(code,context);
  const steps=['post-processing:clinvar-latest','sample-step:review-tsv','post-processing:cnv-rescue','promote-output'];
  const values=steps.map(step=>context._dragenProgressPercent({state:'running',step,with_prs:true}));
  for(const pct of [95,96,98,99])values.push(context._dragenProgressPercent({state:'running',step:'post-processing:prs',prs_execution:'postprocessing',prs_progress_pct:pct}));
  assert(values.every(p=>p<100));assert(values.every((p,i)=>i===0 || p>=values[i-1]));
  assert.equal(context._dragenProgressPercent({state:'done',step:'done',with_prs:true}),100);
  assert.equal(context._dragenProgressPercent({state:'running',step:'promote-output',with_prs:false}),99);
  assert.match(context._dragenJobStepLabel({step:'post-processing:prs',prs_sample:'S1',prs_phase:'祖源與品質檢查'}),/PRS 後處理 · S1 · 祖源與品質檢查/);
});

test('live Nextflow block mirrors the console layout and hides after the Nextflow phase', () => {
  const code = source.slice(source.indexOf('function _dragenNextflowStatusText('), source.indexOf('function _dragenFmtEta('));
  const el = {hidden: true, children: [], replaceChildren(...nodes) { this.children = nodes; }};
  const context = {
    document: {getElementById: () => el, createElement: () => ({className: '', textContent: ''})},
  };
  vm.runInNewContext(code, context);
  const rows = [
    {process: 'ADD_CALLERS_TAG', name: 'PREPARE_VCF:ADD_CALLERS_TAG', done: 36, total: 36, cached: 36, running: 0, status: 'done'},
    {process: 'VEP_ANNOTATE', name: 'ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE', done: 12, total: 36, cached: 0, running: 4, status: 'running'},
    {process: 'ACMG_CLASSIFY', name: 'ANNOTATE_SNV:ACMG_CLASSIFY', done: 0, total: 36, cached: 0, running: 0, status: 'pending'},
  ];
  const lines = context._dragenNextflowStatusText(rows);
  assert.equal(lines[0].text, '[  ✔ ] PREPARE_VCF:ADD_CALLERS_TAG            | 36 of 36, cached: 36 ✔');
  assert.equal(lines[1].text, '[  4 ] ANNOTATE_SNV:SNV_ANNOTATE:VEP_ANNOTATE | 12 of 36');
  assert.equal(lines[2].text, '[  - ] ANNOTATE_SNV:ACMG_CLASSIFY             | 0 of 36');
  context._dragenRenderNextflowStatus({state: 'running', step: 'nextflow', progress_phase: 'nextflow', nextflow_processes: rows});
  assert.equal(el.hidden, false);
  assert.deepEqual(el.children.map(c => c.className), ['nf-done', 'nf-running', 'nf-pending']);
  context._dragenRenderNextflowStatus({state: 'running', step: 'post-processing:genebe', progress_phase: 'postprocessing', nextflow_processes: rows});
  assert.equal(el.hidden, true);
});

test('the Log toggle opens and closes the Nextflow block and the log together', () => {
  const code = source.slice(source.indexOf('function _dragenNextflowStatusText('), source.indexOf('function _dragenFmtEta('))
    + source.slice(source.indexOf('function _toggleDragenLog('), source.indexOf('function _dragenStartPolling('));
  const els = {};
  for (const id of ['dragen-job-output', 'dragen-nf-section', 'dragen-nf-status', 'dragen-log-section', 'dragen-job-log'])
    els[id] = {hidden: true, scrollTop: 0, scrollHeight: 10, replaceChildren() {}};
  els['dragen-job-log-toggle'] = {textContent: '▶ Log', setAttribute() {}};
  const context = {document: {getElementById: id => els[id], createElement: () => ({})}};
  vm.runInNewContext(code, context);
  const running = {state: 'running', step: 'nextflow', progress_phase: 'nextflow',
    nextflow_processes: [{process: 'VEP_ANNOTATE', done: 1, total: 2, status: 'running', running: 1}]};
  context._dragenRenderNextflowStatus(running);
  assert.equal(els['dragen-job-output'].hidden, true);            // collapsed by default
  context._toggleDragenLog();
  assert.equal(els['dragen-job-output'].hidden, false);           // both shown
  assert.equal(els['dragen-nf-section'].hidden, false);
  assert.equal(els['dragen-log-section'].hidden, false);
  context._toggleDragenLog();
  assert.equal(els['dragen-job-output'].hidden, true);            // both hidden together
  context._toggleDragenLog();
  context._dragenRenderNextflowStatus({...running, step: 'post-processing:genebe', progress_phase: 'postprocessing'});
  assert.equal(els['dragen-job-output'].hidden, false);           // log stays open after Nextflow
  assert.equal(els['dragen-nf-section'].hidden, true);
});
