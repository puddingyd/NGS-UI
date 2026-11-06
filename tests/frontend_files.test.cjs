const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname,'../frontend/files.js'),'utf8');
const modulePromise = import('data:text/javascript;base64,'+Buffer.from(source).toString('base64'));
const icons={copy:'<svg data-copy-icon></svg>',external:'<svg data-link-icon></svg>'};
const file=n=>({id:String(n),name:`S1_L${n}.fastq.gz`,path:`/raw/S1_L${n}.fastq.gz`});

test('six cells, two initial FASTQs, all remaining files in a native collapsible section',async()=>{
  const {renderFileRows}=await modulePromise;
  const html=renderFileRows([{sample_id:'S1',fastq:Array.from({length:16},(_,i)=>file(i+1)),bam:[{id:'b',name:'S1.bam.bai',path:'/bam/S1.bam.bai'}]}],icons);
  assert.equal((html.match(/<td[ >]/g)||[]).length,5);
  assert.equal((html.split('<details')[0].match(/file-entry/g)||[]).length,2);
  assert.equal((html.match(/file-entry/g)||[]).length,17);
  assert.match(html,/>more</);assert.match(html,/>less</);
  assert.match(html,/data-copy="\/bam\/S1.bam.bai"/);
  assert.doesNotMatch(html,/缺|不適用|repeats/);
});

test('empty cells stay blank and single/paired FASTQ has no more control',async()=>{
  const {renderFileRows}=await modulePromise;
  const html=renderFileRows([{sample_id:'RAW_ONLY',fastq:[file(1)]}],icons);
  assert.doesNotMatch(html,/<details|more|缺|找不到|—/);
  assert.match(html,/<td><\/td><td><\/td><td><\/td>/);
});

test('filenames and copy paths are escaped; dated directory links retain target and copy path',async()=>{
  const {renderFileRows}=await modulePromise;
  const html=renderFileRows([{sample_id:'<S>',fastq:[{id:'a/b',name:'<img>.fastq.gz',path:'" onclick="alert(1)'}],tertiary:[{root:'hash',path:'S-dragen',name:'S-dragen',created_date:'2026-11-06',created_at:'2026-11-06T00:00:00+08:00',absolute_path:'/tertiary/S-dragen'}]}],icons);
  assert.doesNotMatch(html,/<img>| onclick="/);
  assert.match(html,/&lt;img&gt;/);
  assert.match(html,/download\/a%2Fb/);
  assert.match(html,/browse\?root=hash&amp;path=S-dragen/);
  assert.match(html,/target="_blank" rel="noopener"/);
  assert.doesNotMatch(html,/data-link-icon|clinvar-external-link/);
  assert.match(html,/>2026-11-06<\/a>/);
  assert.match(html,/data-copy="\/tertiary\/S-dragen"/);
});

test('unavailable birth time keeps a usable folder name and path copy',async()=>{
  const {renderFileRows}=await modulePromise;
  const html=renderFileRows([{sample_id:'S',tertiary:[{root:'h',path:'S-dragen',name:'S-dragen',created_date:null,absolute_path:'/tertiary/S-dragen'}]}],icons);
  assert.match(html,/>S-dragen<\/a>/);
  assert.match(html,/data-copy="\/tertiary\/S-dragen"/);
  assert.doesNotMatch(html,/1970|undefined|缺/);
});
