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
