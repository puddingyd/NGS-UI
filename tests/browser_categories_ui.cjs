// Full-page synthetic interaction; run with NODE_PATH pointing at Playwright.
const {chromium}=require('playwright');
const fs=require('node:fs'),path=require('node:path'),http=require('node:http'),assert=require('node:assert/strict');
const repo=path.resolve(__dirname,'..');
let categories=['GI','Neurology','MCA'], counts={GI:4,Neurology:2,MCA:1}, revision='v1', saves=0, rejectSave=false, delaySave=false;
const server=http.createServer((req,res)=>{
  const p=new URL(req.url,'http://localhost').pathname;
  if(p.startsWith('/api/')) {
    if(p==='/api/options/categories') {
      let body='';req.on('data',b=>body+=b);req.on('end',async()=>{
        saves++;
        if(delaySave)await new Promise(r=>setTimeout(r,350));
        if(rejectSave){res.writeHead(409,{'Content-Type':'application/json'});res.end(JSON.stringify({detail:'分類清單已被更新，請關閉後重新開啟管理分類'}));return;}
        const payload=JSON.parse(body);assert.equal(payload.category_revision,revision);
        const mapping=Object.fromEntries(categories.map(n=>[n,''])),next={};
        for(const row of payload.categories){if(row.original!==null)mapping[row.original]=row.name;next[row.name]=row.original===null?0:counts[row.original];}
        categories=Object.keys(next).sort((a,b)=>next[b]-next[a]);counts=next;revision+='x';
        res.setHeader('Content-Type','application/json');res.end(JSON.stringify({category_options:categories,category_counts:counts,category_revision:revision,category_changes:mapping,updated_samples:5}));
      });return;
    }
    let value={};
    if(p==='/api/auth/me')value={username:'SYNTHETIC'};
    else if(p==='/api/samples')value=[{sample_id:'S1',lis_id:'S1',category:'GI',test_type:'WES'}];
    else if(p==='/api/options')value={category_options:categories,category_counts:counts,category_revision:revision};
    else if(p==='/api/samples/unregistered')value={items:[],updated_at:1};
    else if(p==='/api/dragen/jobs')value={jobs:[]};
    res.setHeader('Content-Type','application/json');res.end(JSON.stringify(value));return;
  }
  const target=path.resolve(repo,'frontend',p==='/'?'index.html':p.slice(1));
  if(!target.startsWith(path.join(repo,'frontend')+path.sep)||!fs.existsSync(target)){res.writeHead(404);res.end();return;}
  res.setHeader('Content-Type',target.endsWith('.js')?'text/javascript':target.endsWith('.css')?'text/css':'text/html; charset=utf-8');
  let body=fs.readFileSync(target,'utf8');
  if(p==='/app.js')body+='\nwindow.categoryTest={state,populateCategorySelect,showModal};';
  res.end(body);
});
(async()=>{
  await new Promise(r=>server.listen(0,'127.0.0.1',r));
  const browser=await chromium.launch({executablePath:process.env.CHROME_EXECUTABLE_PATH||'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',headless:true});
  try{
    const page=await browser.newPage({viewport:{width:1280,height:950}}),errors=[];
    page.on('pageerror',e=>errors.push(String(e)));
    await page.goto(`http://127.0.0.1:${server.address().port}/`);
    await page.waitForFunction(()=>window.categoryTest?.state.options);
    const checkLabelButton=async button=>{
      const layout=await button.evaluate(el=>{
        const label=el.previousElementSibling.getBoundingClientRect(),icon=el.getBoundingClientRect();
        return {gap:icon.left-label.right,center:Math.abs((icon.top+icon.bottom-label.top-label.bottom)/2),
          nextToLabel:el.parentElement.classList.contains('category-label-tools'),width:icon.width};
      });
      assert(layout.nextToLabel);assert(layout.gap>=0&&layout.gap<=4);assert(layout.center<=2);assert(layout.width<=18);
    };
    const modal=page.locator('#category-manager-modal');
    const open=async()=>{await page.locator('#new-case-modal [data-manage-categories]').click();await page.waitForFunction(()=>!document.getElementById('category-manager-save').disabled);};
    await page.locator('#btn-new-case').click();
    await checkLabelButton(page.locator('#new-case-modal [data-manage-categories]'));
    await page.locator('#new-case-modal').screenshot({path:'/private/tmp/ngs-category-new-case-label.png'});
    await page.locator('#new-case-name').fill('保留未送出的姓名');
    await page.locator('#new-case-category').selectOption('GI');
    await open();
    assert.deepEqual(await modal.locator('.category-manager-count').allTextContents(),['4','2','1']);
    await modal.locator('input').first().fill('Unsaved');
    await page.keyboard.press('Escape');
    assert(await modal.isHidden());assert(await page.locator('#new-case-modal').isVisible());
    assert.equal(await page.locator('#new-case-name').inputValue(),'保留未送出的姓名');
    assert.equal(saves,0);
    assert(await page.locator('#new-case-modal [data-manage-categories]').evaluate(e=>e===document.activeElement));
    await open();assert.equal(await modal.locator('input').first().inputValue(),'GI');
    // A save must translate an unsaved category selection without dropping other report edits.
    await page.evaluate(()=>{
      const {state,populateCategorySelect}=window.categoryTest;
      state.reports={category:'GI',comment:'未儲存的判讀文字',edits:{v1:{comment:'保留'}},category_revision:'v1'};
      state.data={meta:{Category:'GI'}};
      populateCategorySelect(document.getElementById('m-category'),'GI');
    });
    await modal.locator('input').first().fill('Gastro');
    await modal.locator('[data-remove-category="2"]').click();
    assert.match(await modal.locator('#category-manager-status').innerText(),/1 個個案將改為未分類/);
    await modal.locator('#category-manager-add').click();
    await modal.locator('input').last().fill('New');
    delaySave=true;
    await modal.locator('#category-manager-save').click();
    await page.keyboard.press('Escape');assert(await modal.isVisible());
    await modal.waitFor({state:'hidden'});
    assert.equal(await page.locator('#new-case-category').inputValue(),'Gastro');
    assert.equal(await page.locator('#m-category').inputValue(),'Gastro');
    assert.equal(await page.locator('#new-case-name').inputValue(),'保留未送出的姓名');
    assert.equal(await page.evaluate(()=>window.categoryTest.state.reports.comment),'未儲存的判讀文字');
    assert.equal(await page.evaluate(()=>window.categoryTest.state.reports.edits.v1.comment),'保留');
    assert.equal(await page.evaluate(()=>window.categoryTest.state.index[0].Category),'Gastro');
    assert.deepEqual(categories,['Gastro','Neurology','New']);
    // Duplicate names and server conflicts stay in the editor with their drafts.
    await open();await modal.locator('input').last().fill('gastro');
    await modal.locator('#category-manager-save').click();
    assert.match(await modal.locator('#category-manager-error').innerText(),/不可重複/);assert.equal(saves,1);
    await modal.locator('input').last().fill('New');rejectSave=true;
    await modal.locator('#category-manager-save').click();
    await modal.locator('#category-manager-error').waitFor({state:'visible'});
    assert.match(await modal.locator('#category-manager-error').innerText(),/已被更新/);
    assert.equal(await modal.locator('input').first().inputValue(),'Gastro');rejectSave=false;
    await modal.locator('[data-close]').click();
    // The basic-info button opens the same editor; deletion translates the active selection to blank.
    await page.locator('#new-case-modal [data-close="new-case-modal"]').first().click();
    await page.evaluate(()=>document.getElementById('sample-card').classList.remove('hidden'));
    await checkLabelButton(page.locator('#sample-card [data-manage-categories]'));
    await page.locator('#sample-card').screenshot({path:'/private/tmp/ngs-category-basic-label.png'});
    await page.locator('#sample-card [data-manage-categories]').click();
    await page.waitForFunction(()=>!document.getElementById('category-manager-save').disabled);
    await modal.locator('[data-remove-category="0"]').click();
    await modal.locator('#category-manager-save').click();await modal.waitFor({state:'hidden'});
    assert.equal(await page.locator('#m-category').inputValue(),'');
    assert.equal(await page.evaluate(()=>window.categoryTest.state.reports.category),'');
    // Reopen and inspect real desktop/mobile layout, focus containment, and persistence after reload.
    await page.locator('#sample-card [data-manage-categories]').click();
    await page.waitForFunction(()=>!document.getElementById('category-manager-save').disabled);
    await modal.locator('#category-manager-save').focus();await page.keyboard.press('Tab');
    assert(await modal.locator('input').first().evaluate(e=>e===document.activeElement));
    const out='/private/tmp/ngs-category-ui';fs.mkdirSync(out,{recursive:true});
    await page.screenshot({path:path.join(out,'desktop.png')});
    await page.setViewportSize({width:390,height:844});
    assert(await modal.locator('.category-manager-card').evaluate(el=>el.scrollWidth<=el.clientWidth));
    assert(await modal.locator('#category-manager-save').isVisible());
    await page.screenshot({path:path.join(out,'mobile.png')});
    await page.reload();await page.waitForFunction(()=>window.categoryTest?.state.options);
    assert.deepEqual(await page.evaluate(()=>window.categoryTest.state.options.category_options),['Neurology','New']);
    assert.deepEqual(errors,[]);
    console.log('Category editor desktop/mobile, nested modal, save/cancel, rename/delete, conflict and reload checks passed.');
  }finally{await browser.close();server.close();}
})().catch(e=>{console.error(e);process.exitCode=1;server.close();});
