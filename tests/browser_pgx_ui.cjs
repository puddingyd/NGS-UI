// Synthetic desktop/mobile interaction. No patient fixtures or remote services.
const {chromium}=require('playwright');
const fs=require('node:fs'),path=require('node:path'),http=require('node:http'),assert=require('node:assert/strict');
const {execFileSync}=require('node:child_process');
const repo=path.resolve(__dirname,'..');
const fixture=JSON.parse(execFileSync(process.env.PYTHON||'python3',['-c',`
import copy,json
from app.services import pgx_concordance as pc,docx_export
calls={'CYP2C9':('*1/*33','*1/*1'),'CYP2D6':('*36+*10/*36+*10','*10+*36/*10+*36'),
       'CYP2B6':('*1/*6','*1/*6;*4/*9'),'NAT2':('*4/*7','*7A/*13A;*7B/*4'),
       'CFTR':('Reference/Reference','ivacaftornon-responsiveCFTRsequence/ivacaftornon-responsiveCFTRsequence'),
       'G6PD':('B (reference)/B (reference)','B(wildtype)'),
       'NUDT15':('*1/*3','*1/*2;*3/*6'),
       'DPYD':('c.1627A>G (*5)','c.1627A>G/c.1896T>C'),
       'RYR1':('Reference/Reference','Reference/Reference')}
p={'genes':{},'pharmcat_available':True,'recommendations_pending':[],'concordance':{'genes':{}},'guideline_annotations':[]}
for g,(label,target) in calls.items():
 a=label.split('/')
 d={'label':label,'allele1_name':a[0],'allele2_name':a[1] if len(a)>1 else '', 'phenotypes':['Intermediate Metabolizer'] if g in ['CYP2C9','CYP2B6','CYP2D6','NAT2','NUDT15'] else [],'definition_version':'2026-02-09','source_diplotypes':[{'label':label}]}
 if g=='NUDT15':
  d.update(definition_source='CLINPGX',definition_version='2026-02-09-10-28',comparison_variants=[
   {'chr':'chr13','pos':48037782,'rsid':'rs746071566','reference':'AGGAGTC','call':'AGGAGTC/AGGAGTCGGAGTC'},
   {'chr':'chr13','pos':48045719,'rsid':'rs116855232','reference':'C','call':'C/T'}])
 if g=='DPYD':
  d.update(source_diplotypes=[{'label':label,'allele1_name':label},{'label':'c.1896T>C','allele1_name':'c.1896T>C'}],comparison_variants=[
   {'allele_names':[label],'call':'T/C','reference':'T'}, {'allele_names':['c.1896T>C'],'call':'A/G','reference':'A'}])
 payload={'gene':g,'diplotype':label,'details':d};p['genes'][g]=payload
 dr={'genotype':target,'starAlleleDatabaseSource':'PharmGKB Database Version: Snapshot-2022.01.01'}
 if g in ['CYP2C9','CFTR','RYR1']:dr['missingVariantSites']=[{'id':'synthetic:LowGQ'}]
 row=pc.compare_gene(g,payload,dr);row.update(fingerprint='synthetic',selected_source='platform',stale=False)
 p['concordance']['genes'][g]=row
p['report_view']=docx_export.build_pgx_report_view(p)
chosen=copy.deepcopy(p);pc._apply(chosen,'CYP2C9','*1/*1','Normal Metabolizer')
chosen['concordance']['genes']['CYP2C9'].update(selected_source='dragen',selected_genotype='*1/*1',reviewer='reviewer')
chosen['report_view']=docx_export.build_pgx_report_view(chosen)
print(json.dumps({'platform':p,'dragen':chosen}))
`],{cwd:repo,encoding:'utf8',env:{...process.env,PYTHONPATH:path.join(repo,'backend')}}));
const app=fs.readFileSync(path.join(repo,'frontend/app.js'),'utf8');
const code=app.slice(app.indexOf('function renderPharmcatBlock('),app.indexOf('const ROH_GRCH38_CHROM_LENGTHS'));
const harness=`const state={currentLIS:'SYNTHETIC',_auxLoadToken:1,index:[],data:{pgx:${JSON.stringify(fixture.platform)}}};
const toggledBlocks=new Set();const currentSampleTestType=()=> 'TITAN-WGS';
const escapeHtml=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));const escapeAttr=escapeHtml;
const apiPut=async(url,body)=>{const r=await fetch('/api'+url,{method:'PUT',body:JSON.stringify(body)});if(!r.ok)throw Error('failed');return r.json();};
const apiFetch=async url=>(await fetch('/api'+url)).json();
${code}
renderPharmcatBlock('sec-pharmcat');renderPharmcatBlock('cat-pharmcat-c');`;
let chosen='platform',rejectNext=false;
const server=http.createServer((req,res)=>{
 if(req.url==='/'){res.setHeader('Content-Type','text/html;charset=utf-8');res.end(`<!doctype html><meta name="viewport" content="width=device-width,initial-scale=1"><link rel="stylesheet" href="/style.css"><style>html,body{height:auto;overflow:visible}body{display:block;margin:16px}.case{max-width:1050px;margin:auto}</style><main class="case"><p>合成個案 UI 測試</p><section id="cat-pharmcat-c" class="cat-block"></section><section id="sec-pharmcat" class="cat-block"></section></main><script src="/harness.js"></script>`);return;}
 if(req.url==='/style.css'){res.setHeader('Content-Type','text/css');res.end(fs.readFileSync(path.join(repo,'frontend/style.css')));return;}
 if(req.url==='/harness.js'){res.setHeader('Content-Type','text/javascript');res.end(harness);return;}
 if(req.url.endsWith('/choice')){let body='';req.on('data',chunk=>body+=chunk);req.on('end',()=>{if(rejectNext){rejectNext=false;res.writeHead(409);res.end('{}');return;}const p=JSON.parse(body);assert(p.confirmed);chosen=p.source;res.setHeader('Content-Type','application/json');res.end(JSON.stringify({pgx:fixture[chosen],pharmcat:fixture[chosen]}));});return;}
 res.setHeader('Content-Type','application/json');res.end(JSON.stringify({pgx:fixture[chosen]}));
});
(async()=>{
 await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
 const browser=await chromium.launch({headless:true,executablePath:process.env.CHROME_EXECUTABLE_PATH||'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'});
 try{
  const page=await browser.newPage({viewport:{width:1280,height:950}}),errors=[];page.on('pageerror',e=>errors.push(String(e)));
  await page.goto(`http://127.0.0.1:${server.address().port}/`);
  const analysis=page.locator('#cat-pharmcat-c');
  const card=g=>analysis.locator('.pgx-analysis-gene').filter({has:page.locator('.pgx-gene-head strong',{hasText:new RegExp('^'+g+'$')})});
  assert.equal(await analysis.locator('.js-pgx-result-choice').count(),1);
  for(const g of ['CYP2D6','NAT2','G6PD'])assert.equal(await card(g).locator('.pgx-comparison-badge').count(),0);
  for(const g of ['NUDT15','DPYD'])assert.equal(await card(g).locator('.pgx-comparison').innerText(),'相位未定');
  assert.match(await card('DPYD').locator('.pgx-gene-genotype').innerText(),/c\.1627A>G.*c\.1896T>C/);
  assert.equal(await card('CYP2B6').locator('.pgx-comparison').innerText(),'相位未定');
  assert.equal(await card('CFTR').locator('.pgx-comparison').innerText(),'LowGQ');
  await card('CYP2C9').locator('select').selectOption('dragen:*1/*1');
  await page.waitForFunction(()=>document.querySelector('#cat-pharmcat-c .js-pgx-result-choice')?.value==='dragen:*1/*1'&&!document.querySelector('#cat-pharmcat-c .js-pgx-result-choice').disabled);
  assert.equal(await card('CYP2C9').locator('.pgx-gene-genotype').innerText(),'*1/*1');
  assert.match(await page.locator('#sec-pharmcat').innerText(),/已人工確認採用 DRAGEN/);
  assert.match(await analysis.innerText(),/用藥建議待重新註解/);
  if(process.env.PGX_SCREENSHOT_DIR){fs.mkdirSync(process.env.PGX_SCREENSHOT_DIR,{recursive:true});await page.screenshot({path:path.join(process.env.PGX_SCREENSHOT_DIR,'pgx-desktop.png'),fullPage:true});}
  await page.setViewportSize({width:390,height:844});
  assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
  if(process.env.PGX_SCREENSHOT_DIR)await card('CYP2C9').screenshot({path:path.join(process.env.PGX_SCREENSHOT_DIR,'pgx-mobile.png')});
  await card('CYP2C9').locator('select').selectOption('platform');
  await page.waitForFunction(()=>document.querySelector('#cat-pharmcat-c select')?.value==='platform'&&!document.querySelector('#cat-pharmcat-c select').disabled);
  assert.equal(await card('CYP2C9').locator('.pgx-gene-genotype').innerText(),'*1/*33');
  rejectNext=true;page.once('dialog',dialog=>dialog.accept());
  await card('CYP2C9').locator('select').selectOption('dragen:*1/*1');
  await page.waitForFunction(()=>document.querySelector('#cat-pharmcat-c select')?.value==='platform'&&!document.querySelector('#cat-pharmcat-c select').disabled);
  assert.equal(await card('CYP2C9').locator('.pgx-gene-genotype').innerText(),'*1/*33');
  assert.deepEqual(errors,[]);console.log('PGx desktop/mobile selection, restore, failed save and quiet badges PASS');
 }finally{await browser.close();server.close();}
})().catch(error=>{console.error(error);server.close();process.exitCode=1;});
