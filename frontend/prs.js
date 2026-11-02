// Separate genotype scoring from clinical risk assessment. Values are server-owned.
let sampleId = null, generation = 0, controller = null, ready = false;
let displayedResult = null, timer = null, currentJob = null, dirtyClinical = false, busy = false;
const byId = id => document.getElementById(id);
const pending = '請輸入臨床資訊後計算';
const active = job => ['queued','running','cancelling'].includes(job?.status);
const today = () => {const d=new Date();return `${d.getFullYear()}-${String(d.getMonth()+1).padStart(2,'0')}-${String(d.getDate()).padStart(2,'0')}`;};
export function formatRisk(value) {
  return typeof value==='number' && Number.isFinite(value) && value>=0 && value<=1 ? `${(value*100).toFixed(2)}%` : '無法估算';
}
function resetRisks() {
  displayedResult=null;byId('prs-result').hidden=true;
  for (const cell of document.querySelectorAll('#prs-disease-rows [data-risk]'))cell.textContent=pending;
}
function renderRows(diseases=[{id:'CHD',name:'冠心病（CHD）',score:null}]) {
  const body=byId('prs-disease-rows');body.replaceChildren();
  for(const disease of diseases) {
    const row=document.createElement('tr');row.dataset.disease=disease.id;
    const score=disease.score;
    for(const value of [disease.name, score ? Number(score.raw_score).toFixed(6) : '尚未完成分析', score ? `第 ${Number(score.percentile).toFixed(2)} 百分位` : '尚未完成分析']) {
      const cell=document.createElement('td');cell.textContent=value;row.append(cell);
    }
    for(const key of ['ten','eighty']) {const cell=document.createElement('td');cell.dataset.risk=key;cell.textContent=pending;row.append(cell);}
    body.append(row);
  }
}
export function resetPrsSample(sid) {
  sampleId=sid||null;generation++;controller?.abort();clearTimeout(timer);
  ready=false;busy=false;displayedResult=null;currentJob=null;dirtyClinical=false;
  if(!byId('prs-card'))return;
  byId('prs-card').classList.toggle('hidden',!sampleId);
  byId('prs-details').open=false;byId('prs-form').reset();
  byId('prs-result').hidden=true;byId('prs-message').textContent='';
  byId('prs-pipeline-message').textContent='';byId('prs-submit').disabled=true;
  byId('prs-source-label').hidden=true;byId('prs-run').disabled=true;byId('prs-cancel').hidden=true;byId('prs-new-run').hidden=true;
  byId('prs-log-details').open=false;byId('prs-log').textContent='尚無分析紀錄';renderRows();
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
  svg.append(svgNode('title', {}, '目前至80歲的 CHD 累積風險：此個案與 PRS 40–60% 參考組'));
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
  const result=snapshot.result;displayedResult=result;byId('prs-result').hidden=false;
  const row=document.querySelector('#prs-disease-rows [data-disease="CHD"]');
  row.querySelector('[data-risk="ten"]').textContent=result.risk_10_year ? formatRisk(result.risk_10_year.risk) : '超出模型年齡範圍';
  row.querySelector('[data-risk="eighty"]').textContent=formatRisk(result.risk_to_80.risk);
  byId('prs-result-context').textContent=`CHD · ${result.person.assessment_date} 評估 · ${result.person.age}歲 · ${result.person.sex==='male'?'男性':'女性'} · 參考組：10年 ${formatRisk(result.reference.risk_10_year?.risk)}／至80歲 ${formatRisk(result.reference.risk_to_80.risk)}`;
  byId('prs-model').textContent=`台灣 GBD2019 · INTERVENE EstBB · ${result.model_id}`;
  renderCurve(result);byId('prs-limitations').replaceChildren();
  for(const text of result.limitations||[]) {const li=document.createElement('li');li.textContent=text;byId('prs-limitations').append(li);}
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
    ready=data.available;renderRows(data.diseases);resetRisks();renderPipeline(data.pipeline);
    byId('prs-submit').disabled=!ready;
    byId('prs-message').textContent=data.issues.join('；') || (data.stale?'資料已更新，請重新計算風險。':'請輸入臨床資訊後計算風險。');
    if(data.latest && !dirtyClinical) {
      const p=data.latest.result.person;byId('prs-age').value=p.age;byId('prs-sex').value=p.sex;byId('prs-history').value=p.prior_chd===false?'no':'';
      renderResult(data.latest);byId('prs-message').textContent='顯示上次評估；修改資料後請重新計算。';
    }
    if(byId('prs-log-details').open)await loadLog(token);
  } catch(error) {
    if(token===generation && error.name!=='AbortError') {resetRisks();byId('prs-message').textContent=error.message;}
  } finally {if(token===generation)schedule();}
}
async function loadLog(token=generation) {
  if(!currentJob)return;
  try {const result=await request(`/jobs/${currentJob.job_id}/log`);if(token===generation)byId('prs-log').textContent=result.log||'等待執行…';}
  catch(error){if(token===generation && error.name!=='AbortError')byId('prs-log').textContent=error.message;}
}
byId('prs-details')?.addEventListener('toggle',()=>{if(byId('prs-details').open)load();else clearTimeout(timer);});
byId('prs-refresh')?.addEventListener('click',()=>load());
byId('prs-log-details')?.addEventListener('toggle',()=>{if(byId('prs-log-details').open)loadLog();});
document.addEventListener('click',event=>{if(event.target.closest?.('[data-target="prs-card"]'))byId('prs-details').open=true;});
byId('prs-form')?.addEventListener('input',()=>{
  dirtyClinical=true;resetRisks();
  if(busy) {generation++;controller?.abort();busy=false;}
  byId('prs-submit').disabled=!ready;byId('prs-message').textContent='臨床資訊已變更，請重新計算風險。';schedule();
});
window.addEventListener('resize',()=>{if(displayedResult && !byId('prs-result').hidden)renderCurve(displayedResult);});
byId('prs-form')?.addEventListener('submit',async event=>{
  event.preventDefault();if(!sampleId || !ready || busy)return;resetRisks();
  if(byId('prs-history').value!=='no') {byId('prs-message').textContent='僅適用確認尚未罹患 CHD 的個案。';return;}
  const token=++generation;controller?.abort();controller=new AbortController();busy=true;
  byId('prs-submit').disabled=true;byId('prs-message').textContent='計算並保存風險…';
  try {
    const result=await request('/assessments',{age:Number(byId('prs-age').value),sex:byId('prs-sex').value,assessment_date:today(),prior_chd:false});
    if(token!==generation)return;dirtyClinical=false;renderResult(result);byId('prs-message').textContent='已保存本次風險評估。';
  } catch(error){if(token===generation && error.name!=='AbortError')byId('prs-message').textContent=error.message;}
  finally {if(token===generation){busy=false;byId('prs-submit').disabled=!ready;schedule();}}
});
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
