// Separate genotype scoring from clinical risk assessment. Values are server-owned.
let sampleId = null, generation = 0, controller = null, ready = false;
let displayedResult = null, timer = null, currentJob = null, dirtyClinical = false, busy = false;
let assessmentResults = {}, scoredDiseases = [];
const diseases = [
  ['CHD','冠心病（CHD）'],['T2D','第二型糖尿病（T2D）'],['ALL_CANCERS','整體癌症（All cancers）'],
  ['BREAST','乳癌（Breast cancer）'],['AF','心房顫動（AF）'],['EPILEPSY','癲癇（Epilepsy）'],
];
const byId = id => document.getElementById(id);
const pending = '請輸入臨床資訊後計算';
const active = job => ['queued','running','cancelling'].includes(job?.status);
const today = () => {const d=new Date();return `${d.getFullYear()}-${String(d.getMonth()+1).padStart(2,'0')}-${String(d.getDate()).padStart(2,'0')}`;};
export function prsClinicalInformationPending(sid) {
  return sampleId===sid && dirtyClinical;
}
export function formatRisk(value) {
  return typeof value==='number' && Number.isFinite(value) && value>=0 && value<=1 ? `${(value*100).toFixed(2)}%` : '無法估算';
}
function resetRisks() {
  displayedResult=null;assessmentResults={};byId('prs-result').hidden=true;
  for (const cell of document.querySelectorAll('#prs-disease-rows [data-risk]'))cell.textContent=pending;
}
function renderRows(records=[]) {
  scoredDiseases=records.filter(r=>r.score).map(r=>r.id);
  const body=byId('prs-disease-rows');body.replaceChildren();
  for(const [id,name] of diseases) {
    const disease=records.find(r=>r.id===id)||{id,name,score:null};
    const row=document.createElement('tr');row.dataset.disease=disease.id;
    const score=disease.score;
    for(const value of [disease.name, score ? Number(score.raw_score).toFixed(6) : '尚未完成分析', score ? `第 ${Number(score.percentile).toFixed(2)} 百分位` : '尚未完成分析']) {
      const cell=document.createElement('td');cell.textContent=value;row.append(cell);
    }
    for(const key of ['ten','twenty','eighty']) {const cell=document.createElement('td');cell.dataset.risk=key;cell.textContent=pending;row.append(cell);}
    body.append(row);
  }
}
export function resetPrsSample(sid) {
  sampleId=sid||null;generation++;controller?.abort();clearTimeout(timer);
  ready=false;busy=false;displayedResult=null;currentJob=null;dirtyClinical=false;
  if(!byId('prs-card'))return;
  byId('prs-card').classList.toggle('hidden',!sampleId);
  byId('prs-details').open=false;byId('prs-form').reset();
  buildHistory();byId('prs-history-menu').open=false;
  byId('prs-result').hidden=true;byId('prs-message').textContent='';
  byId('prs-score-status').textContent='';renderExplanation();byId('prs-explanation').open=false;
  byId('prs-pipeline-message').textContent='';byId('prs-submit').disabled=true;
  byId('prs-source-label').hidden=true;byId('prs-run').disabled=true;byId('prs-cancel').hidden=true;byId('prs-new-run').hidden=true;
  byId('prs-log-details').open=false;byId('prs-log').textContent='尚無分析紀錄';byId('prs-raw-log').hidden=true;renderRows();
}
async function request(path='', body, signal=controller?.signal) {
  const response=await fetch(`/api/samples/${encodeURIComponent(sampleId)}/prs${path}`,{
    method:body?'POST':'GET',credentials:'same-origin',cache:'no-store',signal,
    headers:{Accept:'application/json',...(body?{'Content-Type':'application/json'}:{})},
    ...(body?{body:JSON.stringify(body)}:{}),
  });
  const data=await response.json().catch(()=>({}));
  if(!response.ok)throw new Error(response.status===401?'請先登入。':(typeof data.detail==='string'?data.detail:'無法取得 PRS 資料'));
  return data;
}
function svgNode(name, attrs, content) {
  const node = document.createElementNS('http://www.w3.org/2000/svg', name);
  for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, String(value));
  if (content !== undefined) node.textContent = content;
  return node;
}

function renderCurve(result) {
  const own = result.risk_to_80.curve;
  const reference = result.reference.risk_to_80.curve;
  const svg = byId('prs-chart');
  svg.replaceChildren();
  const width = Math.max(300, Math.min(720, svg.clientWidth || 720));
  const height = 280, left = 62, right = 22, top = 24, bottom = 44;
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  const first = result.person.age;
  const maximum = Math.max(...own.map(p => p.risk), ...reference.map(p => p.risk));
  const cap = Math.min(1, Math.max(.01, Math.ceil(maximum * 100 / 5) * .05));
  const x = age => left + (age-first)/(80-first)*(width-left-right);
  const y = risk => height-bottom-risk/cap*(height-top-bottom);
  const title=`目前至80歲的 ${result.name||'冠心病（CHD）'} 累積風險：此個案與 PRS 40–60% 參考組`;
  svg.setAttribute('aria-label',title);svg.append(svgNode('title', {},title));
  for (let i=0; i<=4; i++) {
    const value = cap*i/4;
    svg.append(svgNode('line', {x1:left,x2:width-right,y1:y(value),y2:y(value),stroke:'#dfe5ed'}));
    svg.append(svgNode('text', {x:left-9,y:y(value)+4,'text-anchor':'end',class:'prs-axis'}, `${(value*100).toFixed(1)}%`));
  }
  const intermediate = width < 420 ? [Math.round((first+80)/20)*10] : [30,40,50,60,70];
  const ticks = [...new Set([first, ...intermediate.filter(a => a>first+4 && a<76), 80])];
  for (const age of ticks) svg.append(svgNode('text', {x:x(age),y:height-17,'text-anchor':'middle',class:'prs-axis'}, `${age}歲`));
  for (const [curve,color,dash] of [[reference,'#66768d','7 5'],[own,'#126e92',null]]) {
    const attrs = {d:curve.map((p,i)=>`${i?'L':'M'}${x(p.age)},${y(p.risk)}`).join(' '),fill:'none',stroke:color,'stroke-width':3};
    if (dash) attrs['stroke-dasharray'] = dash;
    svg.append(svgNode('path',attrs));
    for (const point of curve) {
      const dot=svgNode('circle',{cx:x(point.age),cy:y(point.risk),r:3,fill:color});
      dot.append(svgNode('title',{},`${point.age}歲：${formatRisk(point.risk)}`));svg.append(dot);
    }
  }
  const body=byId('prs-curve-table');body.replaceChildren();
  own.forEach((point,i)=>{
    const row=document.createElement('tr');
    for (const value of [`${point.age}歲`,formatRisk(point.risk),formatRisk(reference[i]?.risk)]) {
      const cell=document.createElement('td');cell.textContent=value;row.append(cell);
    }
    body.append(row);
  });
}
function renderResult(snapshot) {
  const result=snapshot.result;
  assessmentResults=result.diseases||{CHD:{...result,available:true,name:'冠心病（CHD）'}};
  const select=byId('prs-curve-disease'), previous=select.value;select.replaceChildren();
  for(const [id,name] of diseases) {
    const r=assessmentResults[id];if(!r)continue;
    const row=document.querySelector(`#prs-disease-rows [data-disease="${id}"]`);if(!row)continue;
    row.querySelector('[data-risk="ten"]').textContent=!r.available?r.reason:r.risk_10_year?formatRisk(r.risk_10_year.risk):'超出模型年齡範圍';
    row.querySelector('[data-risk="twenty"]').textContent=!r.available?r.reason:!('risk_20_year' in r)?'請重新計算風險':r.risk_20_year?formatRisk(r.risk_20_year.risk):'超出模型年齡範圍';
    row.querySelector('[data-risk="eighty"]').textContent=r.available?formatRisk(r.risk_to_80.risk):r.reason;
    if(r.available)select.append(new Option(name,id));
  }
  if(assessmentResults[previous]?.available)select.value=previous;
  byId('prs-result').hidden=!select.options.length;
  selectCurve();
}
function selectCurve() {
  const r=assessmentResults[byId('prs-curve-disease').value];if(!r?.available)return;
  displayedResult=r;
  const twenty=!('risk_20_year' in r.reference)?'待重新計算':r.reference.risk_20_year?formatRisk(r.reference.risk_20_year.risk):'不適用';
  byId('prs-result-context').textContent=`${r.name} · ${r.person.assessment_date} 評估 · ${r.person.age}歲 · ${r.person.sex==='male'?'男性':'女性'} · 參考組：10年 ${r.reference.risk_10_year?formatRisk(r.reference.risk_10_year.risk):'不適用'}／20年 ${twenty}／至80歲 ${formatRisk(r.reference.risk_to_80.risk)}`;
  renderCurve(r);
}
function buildHistory() {
  const panel=byId('prs-history-fields');panel.replaceChildren();
  for(const [id,name] of diseases) {
    const label=document.createElement('label');label.append(document.createTextNode(name));
    const select=document.createElement('select');select.id=id==='CHD'?'prs-history':`prs-history-${id}`;select.dataset.history=id;
    for(const [text,value] of [['否','no'],['是','yes'],['尚未確認','']])select.append(new Option(text,value,value==='no',value==='no'));
    label.append(select);panel.append(label);
  }
  updateHistoryCount();
}
function updateHistoryCount() {
  const fields=[...document.querySelectorAll('[data-history]')];
  const count=fields.filter(f=>f.value).length;byId('prs-history-count').textContent=count?`已填 ${count}/6`:'尚未確認';
}
function clinicalHistory() {
  return Object.fromEntries([...document.querySelectorAll('[data-history]')].map(f=>[f.dataset.history,f.value==='no'?false:f.value==='yes'?true:null]));
}
function renderExplanation(explanation={}) {
  for(const [id,texts] of [['prs-method',explanation.method],['prs-limitations',explanation.limitations]]) {
    byId(id).replaceChildren();
    for(const text of texts||[]) {const li=document.createElement('li');li.textContent=text;byId(id).append(li);}
  }
  const sources=byId('prs-sources');sources.replaceChildren();
  for(const source of explanation.sources||[]) {
    if(sources.childNodes.length)sources.append(document.createTextNode(' · '));
    const link=document.createElement('a');link.href=source.url;link.textContent=source.label;link.target='_blank';link.rel='noopener';sources.append(link);
  }
}
function renderScoreStatus(score, records=[]) {
  let text='';
  if(score) {
    if(active(currentJob) && score.job_id!==currentJob.job_id)text='正在重新分析；表格目前保留上一次已發布的計分結果。';
    else if(records.length>1 && records.some(r=>r.score?.qc_ready) && records.some(r=>!r.score?.qc_ready))
      text='PRS 計分已完成；以下疾病尚未通過祖源或品質檢查，無法計算風險：'+records.filter(r=>!r.score?.qc_ready).map(r=>r.name).join('、')+'。其他疾病可依臨床資訊計算。';
    else if(!score.qc_ready)text=active(currentJob)
      ? 'PRS 計分已完成；祖源與品質檢查進行中，完成並符合要求後才能計算風險。'
      : 'PRS 計分已完成；祖源或品質檢查尚未符合要求，目前只能查看分數與百分位。';
    else text='PRS 計分與祖源／品質檢查已完成。';
  }
  byId('prs-score-status').textContent=text;
}
function renderPipeline(pipeline={}) {
  currentJob=pipeline.job||null;
  const labels={queued:'排隊中',running:'執行中',done:'完成',done_with_warnings:'計分完成，請查看提醒',failed:'失敗',cancelled:'已取消'};
  byId('prs-pipeline-message').textContent=[currentJob ? `${labels[currentJob.status]||currentJob.status}：${currentJob.step}` : '尚無平台分析紀錄',currentJob?.error,currentJob?.warning,pipeline.setup_issue].filter(Boolean).join('；');
  byId('prs-run').disabled=!pipeline.can_start;byId('prs-run').textContent=currentJob?'重新分析／重試':'開始 PRS 分析';
  byId('prs-new-run').hidden=!pipeline.can_start || !['failed','cancelled','done_with_warnings'].includes(currentJob?.status);
  byId('prs-cancel').hidden=!active(currentJob);byId('prs-cancel').disabled=false;
  const choices=pipeline.candidates||[], select=byId('prs-source'), previous=select.value;
  select.replaceChildren();
  if(choices.length>1)select.append(new Option('請確認原始 gVCF',''));
  for(const choice of choices)select.append(new Option(choice.path,choice.id));
  if(choices.some(c=>c.id===previous))select.value=previous;
  byId('prs-source-label').hidden=choices.length<=1;
}
function schedule() {
  clearTimeout(timer);
  if(sampleId && byId('prs-details').open && active(currentJob))timer=setTimeout(()=>{if(!busy)load(false);else schedule();},3000);
}
async function load(restore=true) {
  if(!sampleId || busy)return;
  const token=++generation;controller?.abort();controller=new AbortController();
  ready=false;byId('prs-submit').disabled=true;
  try {
    const data=await request();if(token!==generation)return;
    ready=data.available;renderRows(data.diseases);resetRisks();renderPipeline(data.pipeline);renderExplanation(data.explanation);renderScoreStatus(data.score,data.diseases);
    byId('prs-submit').disabled=!ready;
    byId('prs-message').textContent=[data.scope_note,data.issues.join('；') || (data.stale?'資料已更新，請重新計算風險。':'請輸入臨床資訊後計算風險。')].filter(Boolean).join(' ');
    if(data.latest && !dirtyClinical) {
      const p=data.latest.result.person;byId('prs-age').value=p.age;byId('prs-sex').value=p.sex;
      for(const field of document.querySelectorAll('[data-history]')) {const v=p.history?p.history[field.dataset.history]:field.dataset.history==='CHD'?p.prior_chd:false;field.value=v===false?'no':v===true?'yes':'';}updateHistoryCount();
      renderResult(data.latest);
      const missingTwenty=Object.values(assessmentResults).some(r=>r.available&&(!('risk_20_year' in r)||!('risk_20_year' in r.reference)));
      byId('prs-message').textContent=[data.scope_note,missingTwenty?'顯示上次評估；請重新按「風險計算」補上20年風險，不需重跑基因分析。':'顯示上次評估；修改資料後請重新計算。'].filter(Boolean).join(' ');
    }
    if(byId('prs-log-details').open)await loadLog(token);
  } catch(error) {
    if(token===generation && error.name!=='AbortError') {resetRisks();byId('prs-message').textContent=error.message;}
  } finally {if(token===generation)schedule();}
}
async function loadLog(token=generation) {
  if(!currentJob)return;
  try {const result=await request(`/jobs/${currentJob.job_id}/log`);if(token===generation) {
    byId('prs-log').textContent=result.log||'等待執行…';
    const link=byId('prs-raw-log');link.hidden=!result.raw_log_available;
    link.href=`/api/samples/${encodeURIComponent(sampleId)}/prs/jobs/${encodeURIComponent(currentJob.job_id)}/raw-log`;
  }}
  catch(error){if(token===generation && error.name!=='AbortError')byId('prs-log').textContent=error.message;}
}
byId('prs-details')?.addEventListener('toggle',()=>{if(byId('prs-details').open)load();else clearTimeout(timer);});
byId('prs-refresh')?.addEventListener('click',()=>load());
byId('prs-log-details')?.addEventListener('toggle',()=>{if(byId('prs-log-details').open)loadLog();});
document.addEventListener('click',event=>{if(event.target.closest?.('[data-target="prs-card"]'))byId('prs-details').open=true;});
document.addEventListener('click',event=>{const menu=byId('prs-history-menu');if(menu&&!menu.contains(event.target))menu.open=false;});
byId('prs-history-menu')?.addEventListener('keydown',event=>{if(event.key==='Escape'){event.currentTarget.open=false;event.currentTarget.querySelector('summary').focus();}});
byId('prs-form')?.addEventListener('input',()=>{
  dirtyClinical=true;resetRisks();updateHistoryCount();
  if(busy) {generation++;controller?.abort();busy=false;}
  byId('prs-submit').disabled=!ready;byId('prs-message').textContent='臨床資訊已變更，請重新計算風險。';schedule();
});
window.addEventListener('resize',()=>{if(displayedResult && !byId('prs-result').hidden)renderCurve(displayedResult);});
byId('prs-form')?.addEventListener('submit',async event=>{
  event.preventDefault();if(!sampleId || !ready || busy)return;resetRisks();
  if(!scoredDiseases.some(id=>clinicalHistory()[id]===false)) {byId('prs-message').textContent='請至少確認一個已完成分析的疾病目前尚未罹患，才能計算首次發病風險。';return;}
  const token=++generation;controller?.abort();controller=new AbortController();busy=true;
  byId('prs-submit').disabled=true;byId('prs-message').textContent='計算並保存風險…';
  try {
    const result=await request('/assessments',{age:Number(byId('prs-age').value),sex:byId('prs-sex').value,assessment_date:today(),history:clinicalHistory()});
    if(token!==generation)return;dirtyClinical=false;renderResult(result);byId('prs-message').textContent='已保存本次風險評估。';
  } catch(error){if(token===generation && error.name!=='AbortError')byId('prs-message').textContent=error.message;}
  finally {if(token===generation){busy=false;byId('prs-submit').disabled=!ready;schedule();}}
});
byId('prs-curve-disease')?.addEventListener('change',selectCurve);
async function startScore(newRun=false) {
  if(busy)return;
  const token=++generation;controller?.abort();controller=new AbortController();busy=true;
  byId('prs-run').disabled=true;
  try {await request('/jobs',{candidate_id:byId('prs-source').value||null,new_run:newRun});if(token===generation){busy=false;await load(false);}}
  catch(error){if(token===generation){byId('prs-pipeline-message').textContent=error.message;byId('prs-run').disabled=false;}}
  finally {if(token===generation)busy=false;}
}
byId('prs-run')?.addEventListener('click',()=>startScore());
byId('prs-new-run')?.addEventListener('click',()=>startScore(true));
byId('prs-cancel')?.addEventListener('click',async()=>{
  if(!currentJob)return;const token=generation;byId('prs-cancel').disabled=true;
  try {await request(`/jobs/${currentJob.job_id}/cancel`,{});if(token===generation)byId('prs-pipeline-message').textContent='正在取消，已完成的計分結果會保留。';}
  catch(error){if(token===generation)byId('prs-pipeline-message').textContent=error.message;}
});
