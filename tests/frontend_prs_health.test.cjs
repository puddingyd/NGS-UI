const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
const code = source.slice(source.indexOf('async function exportHealthDocx('), source.indexOf('function _reportPrintBlock('));

function setup(overrides={}) {
  const seen={urls:[],downloads:[],alerts:[],options:[]};
  const context={
    state:{currentLIS:'S1',index:[{LIS_ID:'S1',sample_id:'S1-dragen'}]}, API_BASE:'/api',
    apiFetch:async url=>{seen.readinessUrl=url;return {available:true};},
    _pickHealthReportSections:async options=>{seen.options.push(options);return ['prs','pgx'];},
    prsClinicalInformationPending:()=>false,
    _pgxChoiceSaving:false, _refreshPgxForExport:async()=>{},
    document:{querySelectorAll:()=>[]}, flushPendingSave:async()=>true,
    fetch:async url=>{seen.urls.push(url);return {ok:true,blob:async()=>({})};},
    downloadBlob:(blob,name)=>seen.downloads.push(name),alert:text=>seen.alerts.push(text),
    ...overrides,
  };
  vm.runInNewContext(code,context);
  return {context,seen};
}

test('health export checks fresh server readiness and passes selected PRS to the DOCX API',async()=>{
  const {context,seen}=setup();await context.exportHealthDocx();
  assert.equal(seen.readinessUrl,'/samples/S1-dragen/prs/report-readiness');
  assert.equal(seen.options[0].available,true);
  assert.match(seen.urls[0],/sections=prs%2Cpgx$/);
  assert.deepEqual(seen.downloads,['S1-dragen_health.docx']);
});

test('uncomputed clinical edits disable PRS even when a previous saved assessment exists',async()=>{
  let options;
  const {context,seen}=setup({prsClinicalInformationPending:()=>true,
    _pickHealthReportSections:async value=>{options=value;return ['pgx'];}});
  await context.exportHealthDocx();
  assert.equal(options.available,false);assert.match(options.reason,/重新計算/);
  assert.match(seen.urls[0],/sections=pgx$/);
});

test('switching cases while readiness is loading never opens or exports the old case',async()=>{
  let release;
  const {context,seen}=setup({apiFetch:()=>new Promise(resolve=>release=resolve)});
  const pending=context.exportHealthDocx();context.state.currentLIS='S2';release({available:true});await pending;
  assert.equal(seen.options.length,0);assert.equal(seen.urls.length,0);
});

test('clinical edits made while the picker is open block exporting an old PRS assessment',async()=>{
  let changed=false;
  const {context,seen}=setup({prsClinicalInformationPending:()=>changed,
    _pickHealthReportSections:async()=>{changed=true;return ['prs','pgx'];}});
  await context.exportHealthDocx();
  assert.equal(seen.urls.length,0);assert.match(seen.alerts[0],/重新計算/);
});
