// Run with Playwright installed/resolvable (NODE_PATH if needed).
// Optional: CHROME_EXECUTABLE_PATH, PYTHON, PRS_SCREENSHOT_DIR.
// All endpoints are mocked; no real patient data or external server required.
const {chromium}=require('playwright');
const fs=require('node:fs'),assert=require('node:assert/strict'),path=require('node:path');
const http=require('node:http'),{execFileSync}=require('node:child_process');
const repo=path.resolve(__dirname,'..');
// Synthetic model and person only; scientific numerical tests live in Python.
const results=JSON.parse(execFileSync(process.env.PYTHON||'python3',['-c',`
import json
from pathlib import Path
from research.prs.absolute_risk.__main__ import evaluate
p=Path('research/prs/examples')
b=json.loads((p/'synthetic_bundle.json').read_text())
s=json.loads((p/'synthetic_summary.json').read_text())
a=json.loads((p/'synthetic_person.json').read_text())
s['scores'][0].update(eas_percentile=97,intervene_group='>=95')
for row in b['bins']: row['hr']['>=95']=2
print(json.dumps({str(age):{'result':evaluate(b,s,{**a,'age':age,'sex':'male'},True)} for age in [50.5,75]}))
`],{cwd:repo,encoding:'utf8'}));
const html=fs.readFileSync(path.join(repo,'frontend/index.html'),'utf8');
const start=html.indexOf('<section id="prs-card"'),end=html.indexOf('</section>',start)+10;
const preview=`<!doctype html><html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><link rel="stylesheet" href="/frontend/style.css"><link rel="stylesheet" href="/frontend/prs.css"><style>body{margin:20px;display:block}.container{margin:auto;max-width:1100px}.preview-note{text-align:center;color:#8a4b00;margin:14px}</style><p class="preview-note">UI 測試預覽：合成個案，非病人結果</p><main class="container">${html.slice(start,end)}</main><script type="module">import {resetPrsSample} from '/frontend/prs.js';window.resetPrsSample=resetPrsSample;resetPrsSample('SYNTHETIC_CASE');</script></html>`;
const server=http.createServer((req,res)=>{
 if(req.url==='/'){res.setHeader('Content-Type','text/html; charset=utf-8');res.end(preview);return;}
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
 let delay=0,posts=0,started=0,cancelled=0,latest=null,ready=true,job=null;
 const score={raw_score:1.23456789,percentile:97,group:'>=95'};
 await page.route('**/api/samples/*/prs**',async route=>{
  const path=route.request().url(),post=route.request().method()==='POST';let response;
  if(path.endsWith('/assessments')) {posts++;const p=route.request().postDataJSON();response=results[String(p.age)];latest=response;}
  else if(path.endsWith('/jobs')) {started++;job={job_id:'abc',status:'running',step:'原始計分'};response=job;}
  else if(path.endsWith('/cancel')) {cancelled++;job={...job,status:'cancelled',step:'已取消'};response=job;}
  else if(path.endsWith('/log'))response={log:'SYNTHETIC JOB LOG'};
  else response={available:ready,issues:ready?[]:['祖源檢查未完成'],score,diseases:[{id:'CHD',name:'冠心病（CHD）',score}],latest:ready?latest:null,stale:false,pipeline:{job,can_start:!['running','queued'].includes(job?.status),candidates:[{id:'candidate',path:'synthetic.gvcf.gz'}]}};
  if(delay)await new Promise(r=>setTimeout(r,delay));await route.fulfill({json:response}).catch(()=>{});
 });
 await page.goto(`http://127.0.0.1:${server.address().port}/`);await page.locator('#prs-details>summary').click();await page.waitForFunction(()=>!document.querySelector('#prs-submit').disabled);
 assert.equal(await page.locator('#prs-disease-rows tr').count(),1);
 assert.match(await page.locator('#prs-disease-rows').innerText(),/1.234568/);
 assert.deepEqual(await page.locator('[data-risk]').allTextContents(),['請輸入臨床資訊後計算','請輸入臨床資訊後計算']);
 const boxes=await page.locator('#prs-form label,#prs-submit').evaluateAll(es=>es.map(e=>e.getBoundingClientRect().top));assert(Math.max(...boxes)-Math.min(...boxes)<12,'desktop clinical inputs must fit one row');
 await screenshot('table-pending.png');
 await page.locator('#prs-age').fill('50.5');await page.locator('#prs-sex').selectOption('male');await page.locator('#prs-history').selectOption('no');await page.locator('#prs-submit').click();await page.locator('#prs-result').waitFor({state:'visible'});
 const r=latest.result;assert.equal(await page.locator('[data-risk="ten"]').innerText(),`${(r.risk_10_year.risk*100).toFixed(2)}%`);assert.equal(await page.locator('#prs-chart path').count(),2);
 await screenshot('table-result.png');
 await page.setViewportSize({width:390,height:844});await page.waitForTimeout(100);assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),'mobile overflow outside scrollable table');await screenshot('table-mobile.png');
 await page.locator('#prs-age').fill('75');assert.equal(await page.locator('#prs-result').isVisible(),false);assert.equal(await page.locator('[data-risk="ten"]').innerText(),'請輸入臨床資訊後計算');await page.locator('#prs-submit').click();await page.locator('#prs-result').waitFor({state:'visible'});assert.match(await page.locator('[data-risk="ten"]').innerText(),/超出/);
 await page.locator('#prs-history').selectOption('yes');const n=posts;await page.locator('#prs-submit').click();assert.equal(posts,n);
 ready=false;await page.locator('#prs-refresh').click();await page.waitForFunction(()=>document.querySelector('#prs-submit').disabled);assert.match(await page.locator('#prs-disease-rows').innerText(),/1.234568/);assert.equal(await page.locator('#prs-result').isVisible(),false);
 await page.locator('#prs-run').click();await page.locator('#prs-cancel').waitFor({state:'visible'});assert.equal(started,1);
 await page.locator('#prs-log-details>summary').click();await page.waitForFunction(()=>document.querySelector('#prs-log').textContent.includes('SYNTHETIC JOB LOG'));
 await page.locator('#prs-cancel').click();await page.waitForTimeout(100);assert.equal(cancelled,1);
 delay=200;await page.locator('#prs-refresh').click();await page.evaluate(()=>window.resetPrsSample('OTHER_CASE'));await page.waitForTimeout(350);assert.equal(await page.locator('#prs-message').innerText(),'');assert.equal(await page.locator('#prs-result').isVisible(),false);
 assert.equal(errors.length,0,errors.join('\n'));console.log('PASS: table/one-line form/pending cells, 2 curves, mobile scroll, clinical invalidation, no ancestry keeps score, run/log/cancel, sample-switch race');
 } finally {await browser.close();await new Promise(resolve=>server.close(resolve));}
})().catch(e=>{console.error(e);process.exit(1)});
