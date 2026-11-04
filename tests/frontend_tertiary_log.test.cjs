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
