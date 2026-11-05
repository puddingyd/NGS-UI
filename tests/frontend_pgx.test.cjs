const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
const code = source.slice(source.indexOf('function _pgxSourceWarning('), source.indexOf('function renderPharmcatReportBody('));
const row = {status:'different',badges:['結果不一致','LowGQ'],platform_genotypes:['*1/*33'],
  dragen_genotype:'*1/*1',candidates:['*1/*1'],fingerprint:'synthetic',selected_source:'platform'};
const pgx = {pharmcat_available:true,concordance:{genes:{CYP2C9:row}},report_view:{}};

function setup(overrides={}) {
  const seen={alerts:[],render:[],requests:[]};
  const context={
    state:{currentLIS:'SYNTHETIC',_auxLoadToken:1,index:[],data:{pgx:structuredClone(pgx)}},
    document:{querySelectorAll:()=>[],addEventListener:(name,callback)=>{seen.change=callback;}},
    escapeHtml:s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c])),
    alert:s=>seen.alerts.push(s),renderPharmcatBlock:id=>seen.render.push(id),
    apiPut:async(url,body)=>{seen.requests.push({url,body});return {pgx:{selected:true}};},
    apiFetch:async()=>({pgx:structuredClone(pgx)}),...overrides,
  };
  context.escapeAttr=context.escapeHtml;
  vm.createContext(context);vm.runInContext(code,context);
  return {context,seen,change:value=>seen.change({target:{matches:()=>true,dataset:{gene:'CYP2C9'},value}})};
}

test('same genes are silent; mismatch has explicit source selector and escaped evidence',()=>{
  const {context}=setup();
  assert.equal(context._pgxComparisonHtml({concordance:{genes:{NAT2:{status:'same',badges:[]}}}},'NAT2',true),'');
  const html=context._pgxComparisonHtml(pgx,'CYP2C9',true);
  assert.match(html,/人工確認後採用/);assert.match(html,/DRAGEN：\*1\/\*1/);
  assert.doesNotMatch(context._pgxComparisonHtml(pgx,'CYP2C9'),/<select/);
  assert.match(context._pgxComparisonHtml({concordance:{genes:{CYP2C9:{...row,dragen_genotype:'<script>'}}}},'CYP2C9',true),/&lt;script&gt;/);
});

test('source selection is saved immediately and both views receive the server projection',async()=>{
  const {context,seen,change}=setup();await change('dragen:*1/*1');
  assert.equal(seen.requests[0].body.confirmed,true);
  assert.equal(seen.requests[0].body.genotype,'*1/*1');
  assert.equal(context.state.data.pgx.selected,true);
  assert.deepEqual(seen.render,['sec-pharmcat','cat-pharmcat-c']);
});

test('failed save rolls back instead of exporting the unpersisted selection',async()=>{
  const {context,seen,change}=setup({apiPut:async()=>{throw new Error('failed');}});
  await change('dragen:*1/*1');
  assert.equal(context.state.data.pgx.concordance.genes.CYP2C9.selected_source,'platform');
  assert.match(seen.alerts[0],/未儲存/);
});

test('a late save response never overwrites a newly loaded case',async()=>{
  let release;const {context,change}=setup({apiPut:()=>new Promise(resolve=>{release=resolve;})});
  const pending=change('dragen:*1/*1');context.state.currentLIS='OTHER';context.state._auxLoadToken++;
  context.state.data={pgx:{other:true}};release({pgx:{selected:true}});await pending;
  assert.equal(context.state.data.pgx.other,true);
});

test('stale choices have a re-confirmation placeholder and can restore platform',()=>{
  const {context}=setup();
  const html=context._pgxComparisonHtml({concordance:{genes:{CYP2C9:{...row,stale:true,status:'same'}}}},'CYP2C9',true);
  assert.match(html,/value="unconfirmed" disabled selected/);assert.match(html,/value="platform"/);
});

test('exports refuse in-flight saves and refresh saved choices before output',async()=>{
  const {context}=setup();
  vm.runInContext('_pgxChoiceSaving = true',context);
  await assert.rejects(context._refreshPgxForExport(),/正在儲存/);
  vm.runInContext('_pgxChoiceSaving = false',context);
  context.apiFetch=async()=>({pgx:{review_blocked:true}});
  await assert.rejects(context._refreshPgxForExport(),/重新確認/);
});
