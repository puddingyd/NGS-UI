// Run with Playwright installed/resolvable (NODE_PATH if needed).
// Optional: CHROME_EXECUTABLE_PATH, PYTHON, PRS_SCREENSHOT_DIR.
// All endpoints are mocked; no real patient data or external server required.
const {chromium}=require('playwright');
const fs=require('node:fs'),assert=require('node:assert/strict'),path=require('node:path');
const http=require('node:http'),{execFileSync}=require('node:child_process');
const repo=path.resolve(__dirname,'..');
// Synthetic model and person only; scientific numerical tests live in Python.
const fixture=JSON.parse(execFileSync(process.env.PYTHON||'python3',['-c',`
import json
from pathlib import Path
from research.prs.absolute_risk.__main__ import evaluate
from backend.app.services.prs_presentation import explanation
p=Path('research/prs/examples')
b=json.loads((p/'synthetic_bundle.json').read_text())
s=json.loads((p/'synthetic_summary.json').read_text())
a=json.loads((p/'synthetic_person.json').read_text())
s['scores'][0].update(eas_percentile=97,intervene_group='>=95')
for row in b['bins']: row['hr']['>=95']=2
import sys
sys.path[:0]=['backend','tests']
from test_prs_multi import summary as multi_summary, payload as multi_payload
from app.services import prs_multi
ms=multi_summary();mp=multi_payload();pack=prs_multi.bundle()[0]
female=prs_multi.calculator.evaluate(pack,ms,prs_multi.calculator.person_for(ms,mp))
mp['sex']='male';mp['history']['AF']=None
male=prs_multi.calculator.evaluate(pack,ms,prs_multi.calculator.person_for(ms,mp))
print(json.dumps({'multi':{'female':{'result':female},'male':{'result':male},'explanation':prs_multi.explanation()}, 'results':{str(age):{'result':evaluate(b,s,{**a,'age':age,'sex':'male'},True)} for age in [50.5,75]}, 'explanation':explanation()}))
`],{cwd:repo,encoding:'utf8'}));
const results=fixture.results;
const html=fs.readFileSync(path.join(repo,'frontend/index.html'),'utf8');
const start=html.indexOf('<section id="prs-card"'),end=html.indexOf('</section>',start)+10;
assert(html.indexOf('id="prs-card"') < html.indexOf('id="cat-pharmcat-c"'));
assert(html.indexOf('data-target="prs-card"') < html.indexOf('data-target="cat-pharmcat-c"'));
const appSource=fs.readFileSync(path.join(repo,'frontend/app.js'),'utf8');
const pickerSource=appSource.slice(appSource.indexOf('function _pickHealthReportSections('),appSource.indexOf('async function exportHealthDocx('));
const pickerScript=`const escapeHtml=v=>String(v).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));const escapeAttr=escapeHtml;${pickerSource};window.testReportPicker=_pickHealthReportSections;`;
const preview=`<!doctype html><html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><link rel="stylesheet" href="/frontend/style.css"><link rel="stylesheet" href="/frontend/prs.css"><style>body{margin:20px;display:block}.container{margin:auto;max-width:1100px}.preview-note{text-align:center;color:#8a4b00;margin:14px}</style><p class="preview-note">UI 測試預覽：合成個案，非病人結果</p><main class="container">${html.slice(start,end)}<div id="cat-pharmcat-c" class="cat-block"><div class="block-header"><span class="title">PGx / PharmCAT</span></div><div class="block-body open">PGx 樣式比較</div></div></main><script src="/report-picker.js"></script><script type="module">import {resetPrsSample} from '/frontend/prs.js';window.resetPrsSample=resetPrsSample;resetPrsSample('SYNTHETIC_CASE');</script></html>`;
const server=http.createServer((req,res)=>{
 if(req.url==='/'){res.setHeader('Content-Type','text/html; charset=utf-8');res.end(preview);return;}
 if(req.url==='/report-picker.js'){res.setHeader('Content-Type','text/javascript');res.end(pickerScript);return;}
 const types={'/frontend/prs.js':'text/javascript','/frontend/prs.css':'text/css','/frontend/style.css':'text/css'};
 if(!types[req.url]){res.writeHead(404);res.end();return;}
 res.setHeader('Content-Type',types[req.url]);res.end(fs.readFileSync(path.join(repo,req.url)));
});
(async()=>{
 await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
 const browser=await chromium.launch({headless:true,...(process.env.CHROME_EXECUTABLE_PATH?{executablePath:process.env.CHROME_EXECUTABLE_PATH}:{})});
 try {
 const page=await browser.newPage({viewport:{width:1240,height:1000}}),errors=[];
 const screenshot=async name=>{if(process.env.PRS_SCREENSHOT_DIR){fs.mkdirSync(process.env.PRS_SCREENSHOT_DIR,{recursive:true});await page.screenshot({path:path.join(process.env.PRS_SCREENSHOT_DIR,name),fullPage:true});}};page.on('pageerror',e=>errors.push(String(e)));
 let delay=0,posts=0,started=0,cancelled=0,latest=null,ready=true,job=null,multiMode=false,qcFailed=null;
 const score={raw_score:1.23456789,percentile:97,group:'>=95',qc_ready:true,job_id:'published'};
 await page.route('**/api/samples/*/prs**',async route=>{
  const path=route.request().url(),post=route.request().method()==='POST';let response;
  if(path.endsWith('/assessments')) {posts++;const p=route.request().postDataJSON();assert.equal(Object.keys(p.history).length,6);response=multiMode?fixture.multi[p.sex]:results[String(p.age)];latest=response;}
  else if(path.endsWith('/jobs')) {started++;job={job_id:'abc',status:'running',step:'原始計分'};response=job;}
  else if(path.endsWith('/cancel')) {cancelled++;job={...job,status:'cancelled',step:'已取消'};response=job;}
  else if(path.endsWith('/log'))response={log:'SYNTHETIC JOB LOG',raw_log_available:true};
  else response={available:ready,issues:ready?[]:['祖源檢查未完成'],score,explanation:fixture.explanation,diseases:[{id:'CHD',name:'冠心病（CHD）',score}],latest:ready?latest:null,stale:false,pipeline:{job,can_start:!['running','queued'].includes(job?.status),candidates:[{id:'candidate',path:'synthetic.gvcf.gz'}]}};
  if(!post && multiMode && !path.endsWith('/log')) {
    response.diseases=Object.values(fixture.multi.female.result.diseases).map(r=>({id:r.id,name:r.name,score:{...score,qc_ready:score.qc_ready&&r.id!==qcFailed}}));response.explanation=fixture.multi.explanation;
    response.score=response.diseases[0].score;
  }
  if(delay)await new Promise(r=>setTimeout(r,delay));await route.fulfill({json:response}).catch(()=>{});
 });
 await page.goto(`http://127.0.0.1:${server.address().port}/`);await page.locator('#prs-details>summary').click();await page.waitForFunction(()=>!document.querySelector('#prs-submit').disabled);
 assert.match(await page.locator('#prs-details>summary').innerText(),/Polygenic risk rcore \(PRS\)/);
 assert.equal(await page.locator('#prs-card').evaluate(e=>getComputedStyle(e).backgroundColor),'rgb(255, 255, 255)');
 await page.mouse.move(0,0);
 assert.equal(await page.locator('#prs-details>summary').evaluate(e=>getComputedStyle(e).backgroundColor),await page.locator('#cat-pharmcat-c .block-header').evaluate(e=>getComputedStyle(e).backgroundColor));
 await page.evaluate(()=>{window.pickerResult=window.testReportPicker({available:false,reason:'尚未完成風險計算'});});
 const prsOption=page.locator('.health-report-options input[value="prs"]');
 assert.equal(await prsOption.isDisabled(),true);assert.equal(await prsOption.isChecked(),false);
 assert.match(await prsOption.locator('..').getAttribute('class'),/is-disabled/);
 assert.deepEqual(await page.locator('.health-report-options input').evaluateAll(es=>es.map(e=>e.value)),['acmg_sf','hereditary_cancer','stroke','carrier','prs','pgx']);
 await screenshot('report-unavailable.png');await page.locator('[data-act="ok"]').click();
 assert.deepEqual(await page.evaluate(()=>window.pickerResult),['acmg_sf','pgx']);
 await page.evaluate(()=>{window.pickerResult=window.testReportPicker({available:true});});
 assert.equal(await prsOption.isDisabled(),false);assert.equal(await prsOption.isChecked(),true);
 await screenshot('report-ready.png');await page.locator('[data-act="ok"]').click();
 assert.deepEqual(await page.evaluate(()=>window.pickerResult),['acmg_sf','prs','pgx']);
 assert.equal(await page.locator('#prs-disease-rows tr').count(),6);
 assert.match(await page.locator('#prs-disease-rows').innerText(),/1.234568/);
 assert((await page.locator('[data-risk]').allTextContents()).every(s=>s==='請輸入臨床資訊後計算'));
 const boxes=await page.locator('#prs-form > label,#prs-history-menu>summary,#prs-submit').evaluateAll(es=>es.map(e=>e.getBoundingClientRect().top));assert(Math.max(...boxes)-Math.min(...boxes)<12,'desktop clinical inputs must fit one row');
 await page.locator('#prs-explanation>summary').click();
 assert.match(await page.locator('#prs-method').innerText(),/PGS004879/);
 assert.match(await page.locator('#prs-method').innerText(),/610,305/);
 assert.match(await page.locator('#prs-limitations').innerText(),/不穩定心絞痛/);
 assert.match(await page.locator('#prs-limitations').innerText(),/穩定心絞痛/);
 assert.equal(await page.locator('#prs-sources a').count(),3);
 assert.doesNotMatch(await page.locator('#prs-card').innerText(),/Model A|intervene_estbb_tw2019_A_exact_v1/);
 await screenshot('table-pending.png');
 ready=false;score.qc_ready=false;job={job_id:'published',status:'running',step:'祖源與品質檢查'};
 await page.locator('#prs-refresh').click();await page.waitForFunction(()=>document.querySelector('#prs-score-status').textContent.includes('檢查進行中'));
 assert.equal(await page.locator('#prs-submit').isDisabled(),true);
 assert.match(await page.locator('#prs-disease-rows').innerText(),/1.234568/);
 ready=true;score.qc_ready=true;job={job_id:'published',status:'done',step:'完成'};
 await page.locator('#prs-refresh').click();await page.waitForFunction(()=>!document.querySelector('#prs-submit').disabled);
 await page.locator('#prs-age').fill('50.5');await page.locator('#prs-sex').selectOption('male');await page.locator('#prs-history-menu>summary').click();await page.locator('#prs-history').selectOption('no');await page.locator('#prs-history-menu>summary').click();await page.locator('#prs-submit').click();await page.locator('#prs-result').waitFor({state:'visible'});
 const r=latest.result;assert.equal(await page.locator('[data-disease="CHD"] [data-risk="ten"]').innerText(),`${(r.risk_10_year.risk*100).toFixed(2)}%`);assert.equal(await page.locator('#prs-chart path').count(),2);
 await screenshot('table-result.png');
 await page.setViewportSize({width:390,height:844});await page.waitForTimeout(100);assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),'mobile overflow outside scrollable table');await screenshot('table-mobile.png');
 await page.locator('#prs-age').fill('75');assert.equal(await page.locator('#prs-result').isVisible(),false);assert.equal(await page.locator('[data-disease="CHD"] [data-risk="ten"]').innerText(),'請輸入臨床資訊後計算');await page.locator('#prs-submit').click();await page.locator('#prs-result').waitFor({state:'visible'});assert.match(await page.locator('[data-disease="CHD"] [data-risk="ten"]').innerText(),/超出/);
 await page.locator('#prs-history-menu>summary').click();await page.locator('#prs-history').selectOption('yes');await page.locator('#prs-history-menu>summary').click();const n=posts;await page.locator('#prs-submit').click();assert.equal(posts,n);
 ready=false;await page.locator('#prs-refresh').click();await page.waitForFunction(()=>document.querySelector('#prs-submit').disabled);assert.match(await page.locator('#prs-disease-rows').innerText(),/1.234568/);assert.equal(await page.locator('#prs-result').isVisible(),false);
 await page.locator('#prs-run').click();await page.locator('#prs-cancel').waitFor({state:'visible'});assert.equal(started,1);
 assert.match(await page.locator('#prs-score-status').innerText(),/上一次已發布/);
 await page.locator('#prs-log-details>summary').click();await page.waitForFunction(()=>document.querySelector('#prs-log').textContent.includes('SYNTHETIC JOB LOG'));
 assert.match(await page.locator('#prs-raw-log').getAttribute('href'),/SYNTHETIC_CASE\/prs\/jobs\/abc\/raw-log/);
 await page.locator('#prs-cancel').click();await page.waitForTimeout(100);assert.equal(cancelled,1);
 multiMode=true;ready=true;latest=null;job=null;await page.setViewportSize({width:1240,height:1000});
 await page.locator('#prs-refresh').click();await page.waitForFunction(()=>!document.querySelector('#prs-submit').disabled);
 await page.locator('#prs-age').fill('50.5');await page.locator('#prs-sex').selectOption('female');
 await page.locator('#prs-history-menu>summary').click();
 for(const field of await page.locator('[data-history]').all())await field.selectOption('no');
 await page.locator('#prs-history-menu>summary').click();await page.locator('#prs-submit').click();await page.locator('#prs-result').waitFor({state:'visible'});
 assert.equal(await page.locator('#prs-curve-disease option').count(),6);
 assert((await page.locator('[data-risk="eighty"]').allTextContents()).every(s=>/^\d+\.\d+%$/.test(s)));
 await page.locator('#prs-curve-disease').selectOption('BREAST');assert.match(await page.locator('#prs-result-context').innerText(),/乳癌/);
 await screenshot('six-diseases-female.png');
 await page.locator('#prs-sex').selectOption('male');assert.equal(await page.locator('#prs-result').isVisible(),false);
 await page.locator('#prs-history-menu>summary').click();await page.locator('#prs-history-AF').selectOption('');await page.locator('#prs-history-menu>summary').click();
 await page.locator('#prs-submit').click();await page.locator('#prs-result').waitFor({state:'visible'});
 assert.equal(await page.locator('#prs-curve-disease option').count(),4);
 assert.match(await page.locator('[data-disease="BREAST"] [data-risk="ten"]').innerText(),/女性/);
 assert.match(await page.locator('[data-disease="AF"] [data-risk="ten"]').innerText(),/確認/);
 await page.setViewportSize({width:390,height:844});await page.waitForTimeout(100);
 assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth));await screenshot('six-diseases-mobile.png');
 qcFailed='CHD';latest=null;await page.locator('#prs-refresh').click();
 await page.waitForFunction(()=>document.querySelector('#prs-score-status').textContent.includes('其他疾病可依臨床資訊計算'));
 assert.match(await page.locator('#prs-score-status').innerText(),/冠心病/);assert.equal(await page.locator('#prs-submit').isDisabled(),false);
 delay=200;await page.locator('#prs-refresh').click();await page.evaluate(()=>window.resetPrsSample('OTHER_CASE'));await page.waitForTimeout(350);assert.equal(await page.locator('#prs-message').innerText(),'');assert.equal(await page.locator('#prs-result').isVisible(),false);
 assert.equal(errors.length,0,errors.join('\n'));console.log('PASS: six-disease rows/history/curves, per-disease QC, legacy CHD, report selection, mobile scroll, clinical invalidation, no ancestry keeps score, run/log/cancel, sample-switch race');
 } finally {await browser.close();await new Promise(resolve=>server.close(resolve));}
})().catch(e=>{console.error(e);process.exit(1)});
