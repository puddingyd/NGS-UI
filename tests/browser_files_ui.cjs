// Synthetic endpoints; run with NODE_PATH pointing at a Playwright installation.
const {chromium}=require('playwright');
const fs=require('node:fs'),path=require('node:path'),http=require('node:http'),assert=require('node:assert/strict');
const repo=path.resolve(__dirname,'..');
const file=(name,folder='raw')=>({id:Buffer.from(name).toString('hex'),name,path:`/synthetic/${folder}/${name}`});
const fastqs=[];
for(let lane=1;lane<=8;lane++)for(let read=1;read<=2;read++)fastqs.push(file(`S1_S50_L00${lane}_R${read}_001.fastq.gz`));
const samples=[{sample_id:'S1',fastq:fastqs,bam:[file('S1.bam'),file('S1.bam.bai')],vcf:[file('S1.hard-filtered.vcf.gz'),{...file('S1.snv_indel.acmg.tsv','tertiary/S1-dragen/03_acmg'),tertiary_root:'root',relative_path:'S1-dragen/03_acmg/S1.snv_indel.acmg.tsv'}],other:['cnv.vcf.gz','sv.vcf.gz','cnv_sv.vcf.gz','hard-filtered.gvcf.gz','targeted.json'].map(s=>file('S1.'+s)),tertiary:[{name:'S1-dragen',root:'root',path:'S1-dragen',created_date:null,display_date:'2026-11-06',date_source:'tertiary_annotation',date_at:'2026-11-06T00:00:00+08:00',absolute_path:'/synthetic/tertiary/S1-dragen'},{name:'S1-nckuh',root:'root',path:'S1-nckuh',created_date:'2026-11-04',absolute_path:'/synthetic/tertiary/S1-nckuh'}]},
{sample_id:'RAW_ONLY',fastq:[file('RAW_ONLY_R1.fastq.gz'),file('RAW_ONLY_R2.fastq.gz')],bam:[],vcf:[],other:[],tertiary:[]},
...Array.from({length:53},(_,i)=>({sample_id:`WES${i}`,fastq:[],bam:[],vcf:[],other:[],tertiary:[]}))];
const server=http.createServer((req,res)=>{
  const pathname=new URL(req.url,'http://localhost').pathname;
  const relative=pathname==='/'?'index.html':pathname.slice(1);
  const target=path.resolve(repo,'frontend',relative);
  if(!target.startsWith(path.join(repo,'frontend')+path.sep)||!fs.existsSync(target)){res.writeHead(404);res.end();return;}
  res.setHeader('Content-Type',target.endsWith('.js')?'text/javascript':target.endsWith('.css')?'text/css':'text/html; charset=utf-8');
  res.end(fs.readFileSync(target));
});
(async()=>{
  await new Promise(r=>server.listen(0,'127.0.0.1',r));
  const browser=await chromium.launch({executablePath:process.env.CHROME_EXECUTABLE_PATH||'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',headless:true});
  const context=await browser.newContext({viewport:{width:1512,height:982},acceptDownloads:true,permissions:['clipboard-read','clipboard-write']});
  const page=await context.newPage(),errors=[];page.on('pageerror',e=>errors.push(String(e)));
  let refreshes=0, delayNext=false;
  const queries=[];
  await context.route('**/api/**',async route=>{
    const url=new URL(route.request().url()),p=url.pathname;
    if(p.startsWith('/api/files/download/'))return route.fulfill({status:200,headers:{'content-disposition':'attachment; filename="synthetic.fastq.gz"'},body:'SYNTHETIC DOWNLOAD'});
    if(p==='/api/files/tertiary-download')return route.fulfill({status:200,headers:{'content-disposition':'attachment; filename="S1.snv_indel.acmg.tsv"'},body:'CHROM\tPOS\nchr1\t123'});
    if(p==='/api/files/browse')return route.fulfill({contentType:'text/html',body:'<html><body><h1>三級資料</h1><a href="?root=root&path=S1-dragen/03_acmg">03_acmg/</a></body></html>'});
    if(p==='/api/files/refresh'){refreshes++;return route.fulfill({json:{}});}
    if(p==='/api/files'){
      const q=url.searchParams.get('q')||'',offset=Number(url.searchParams.get('offset')||0),limit=50;
      queries.push(q);
      if(delayNext){delayNext=false;await new Promise(r=>setTimeout(r,500));}
      const matched=samples.filter(s=>s.sample_id.includes(q));
      return route.fulfill({json:{total:matched.length,samples:matched.slice(offset,offset+limit),meta:{updated_at:'2026-11-06 02:00:00'}}});
    }
    let value={};
    if(p==='/api/auth/me')value={username:'SYNTHETIC'};
    else if(p==='/api/samples')value=[];
    else if(p==='/api/options')value={category_options:[],tag_suggestions:[]};
    else if(p==='/api/dragen/jobs')value={jobs:[]};
    return route.fulfill({json:value});
  });
  try {
    await page.goto(`http://127.0.0.1:${server.address().port}/`);
    await page.locator('#btn-file-inventory').waitFor({state:'visible'});
    const buttons=await page.locator('.topbar-actions button').allTextContents();
    assert(buttons.indexOf('三級分析')<buttons.indexOf('檔案清單')&&buttons.indexOf('檔案清單')<buttons.indexOf('上傳個案清單'));
    await page.locator('#btn-file-inventory').click();
    await page.waitForTimeout(350);
    assert.equal(await page.locator('#file-inventory-rows tr').count(),0);
    assert.deepEqual(queries,[],'opening must not request all samples');
    await page.locator('#file-inventory-refresh').click();
    await page.waitForFunction(()=>!document.getElementById('file-inventory-refresh').disabled);
    assert.deepEqual(queries,[],'refresh without a query must not load all samples');
    assert.equal(refreshes,1);
    await page.locator('#file-inventory-search').fill('S');
    await page.locator('#file-inventory-rows tr').first().waitFor();
    assert.equal(await page.locator('.file-inventory-table thead th').count(),6);
    assert.equal(await page.locator('#file-inventory-rows tr').count(),50);
    const row=page.locator('#file-inventory-rows tr').filter({has:page.locator('th',{hasText:/^S1$/})});
    assert.equal(await row.locator('td').nth(0).locator('.file-entry:visible').count(),2);
    await row.locator('summary').click();
    assert.equal(await row.locator('td').nth(0).locator('.file-entry:visible').count(),16);
    await row.locator('summary').click();
    assert.equal(await row.locator('td').nth(0).locator('.file-entry:visible').count(),2);
    await row.locator('.btn-copy').first().click();
    assert.equal(await page.evaluate(()=>navigator.clipboard.readText()),fastqs[0].path);
    const downloadPromise=page.waitForEvent('download');await row.locator('.file-entry a').first().click();
    const download=await downloadPromise;assert.equal(download.suggestedFilename(),'synthetic.fastq.gz');
    assert.deepEqual(await row.locator('.file-tertiary-link').allTextContents(),['2026-11-06','2026-11-04']);
    assert.equal(await row.locator('.file-tertiary-link svg').count(),0);
    await row.locator('.file-tertiary-entry .btn-copy').first().click();
    assert.equal(await page.evaluate(()=>navigator.clipboard.readText()),'/synthetic/tertiary/S1-dragen');
    const annotation=row.locator('a').filter({hasText:'S1.snv_indel.acmg.tsv'});
    assert.equal(await annotation.count(),1);
    const annotationDownload=page.waitForEvent('download');await annotation.click();
    assert.equal((await annotationDownload).suggestedFilename(),'S1.snv_indel.acmg.tsv');
    await row.locator('.btn-copy[data-copy$="S1.snv_indel.acmg.tsv"]').click();
    assert.equal(await page.evaluate(()=>navigator.clipboard.readText()),'/synthetic/tertiary/S1-dragen/03_acmg/S1.snv_indel.acmg.tsv');
    assert.match(await row.locator('.file-tertiary-link').first().getAttribute('title'),/三級註解紀錄日期/);
    const popupPromise=page.waitForEvent('popup');await row.locator('.file-tertiary-link').first().click();
    const popup=await popupPromise;await popup.waitForLoadState();assert.match(await popup.locator('body').innerText(),/03_acmg/);await popup.close();
    const out=process.env.FILE_SCREENSHOT_DIR||'/private/tmp/ngs-file-ui-preview';fs.mkdirSync(out,{recursive:true});
    await page.screenshot({path:path.join(out,'desktop.png')});
    await page.locator('#file-inventory-next').click();await page.waitForFunction(()=>document.getElementById('file-inventory-count').textContent==='51–54 / 54');
    await page.locator('#file-inventory-search').fill('RAW_ONLY');await page.waitForFunction(()=>document.querySelectorAll('#file-inventory-rows tr').length===1);
    const emptyCells=await page.locator('#file-inventory-rows tr td').allTextContents();assert(emptyCells.slice(1).every(s=>!s.trim()));
    await page.locator('#file-inventory-refresh').click();await page.waitForFunction(()=>!document.getElementById('file-inventory-refresh').disabled);assert.equal(refreshes,2);
    delayNext=true;
    const pendingSearch=page.waitForRequest(r=>new URL(r.url()).pathname==='/api/files'&&new URL(r.url()).searchParams.get('q')==='S1');
    await page.locator('#file-inventory-search').fill('S1');await pendingSearch;
    await page.locator('#file-inventory-search').fill('');await page.waitForTimeout(650);
    assert.equal(await page.locator('#file-inventory-rows tr').count(),0,'late search must not repopulate cleared rows');
    const queryCount=queries.length;
    await page.locator('#file-inventory-search').fill('   ');await page.waitForTimeout(350);
    assert.equal(queries.length,queryCount,'whitespace must not trigger a query');
    await page.locator('#file-inventory-search').fill('S1');await page.waitForFunction(()=>document.querySelector('#file-inventory-rows th')?.textContent==='S1');
    await page.setViewportSize({width:390,height:844});await page.screenshot({path:path.join(out,'mobile.png')});
    const searchBox=await page.locator('#file-inventory-search').boundingBox();
    assert(searchBox&&searchBox.y>=0&&searchBox.y+searchBox.height<844,'mobile search must remain visible');
    await page.locator('#file-inventory-search').click();
    assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),'table overflow must stay in its scroll container');
    assert(await page.locator('.file-table-scroll').evaluate(e=>e.scrollWidth>e.clientWidth));
    await row.locator('.file-tertiary-entry .btn-copy').first().click();
    assert.equal(await page.evaluate(()=>navigator.clipboard.readText()),'/synthetic/tertiary/S1-dragen');
    await page.keyboard.press('Escape');assert(await page.locator('#file-inventory-modal').evaluate(e=>e.classList.contains('hidden')));
    const countBeforeReopen=queries.length;
    await page.setViewportSize({width:1512,height:982});await page.locator('#btn-file-inventory').click();await page.waitForTimeout(350);
    assert.equal(await page.locator('#file-inventory-search').inputValue(),'');
    assert.equal(await page.locator('#file-inventory-rows tr').count(),0);
    assert.equal(queries.length,countBeforeReopen,'reopening starts blank');
    assert.deepEqual(errors,[]);
    console.log('PASS blank open/reopen/clear, no empty-query requests, late-response guard, dated links and folder path copy; existing file-list interactions');
    console.log(out);
  } finally {await context.close();await browser.close();server.close();}
})().catch(e=>{console.error(e);server.close();process.exitCode=1;});
