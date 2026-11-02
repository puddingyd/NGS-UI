// CHD research panel. All numerical values come from the authenticated API.
let sampleId = null;
let generation = 0;
let controller = null;
let ready = false;
let displayedResult = null;
const byId = id => document.getElementById(id);

export function formatRisk(value) {
  return typeof value === 'number' && Number.isFinite(value) && value >= 0 && value <= 1
    ? `${(value * 100).toFixed(2)}%` : '無法估算';
}

export function resetPrsSample(sid) {
  sampleId = sid || null;
  generation += 1;
  controller?.abort();
  ready = false;
  displayedResult = null;
  const card = byId('prs-card');
  if (!card) return;
  card.classList.toggle('hidden', !sampleId);
  byId('prs-details').open = false;
  byId('prs-form').reset();
  byId('prs-form').hidden = true;
  byId('prs-result').hidden = true;
  byId('prs-message').textContent = '';
  byId('prs-score').textContent = '';
  const now = new Date();
  byId('prs-date').value = `${now.getFullYear()}-${String(now.getMonth()+1).padStart(2,'0')}-${String(now.getDate()).padStart(2,'0')}`;
  byId('prs-date').max = byId('prs-date').value;
  byId('prs-submit').disabled = true;
}

async function request(path = '', body) {
  const response = await fetch(`/api/samples/${encodeURIComponent(sampleId)}/prs${path}`, {
    method: body ? 'POST' : 'GET', credentials: 'same-origin', cache: 'no-store',
    headers: {'Accept': 'application/json', ...(body ? {'Content-Type': 'application/json'} : {})},
    ...(body ? {body: JSON.stringify(body)} : {}), signal: controller.signal,
  });
  const data = await response.json().catch(() => ({}));
  if (!response.ok) throw new Error(response.status === 401 ? '請先登入，再開啟 CHD 風險。' :
    (typeof data.detail === 'string' ? data.detail : '無法取得 CHD 風險，請稍後重試。'));
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
  const result = snapshot.result;
  displayedResult = result;
  byId('prs-result').hidden = false;
  byId('prs-risk-ten').textContent = formatRisk(result.risk_10_year?.risk);
  byId('prs-risk-eighty').textContent = formatRisk(result.risk_to_80?.risk);
  byId('prs-ten-reference').textContent = result.risk_10_year
    ? `參考組 ${formatRisk(result.reference.risk_10_year.risk)}` : '完整10年期間超過80歲';
  byId('prs-eighty-reference').textContent = `參考組 ${formatRisk(result.reference.risk_to_80.risk)}`;
  byId('prs-result-context').textContent = `${result.person.assessment_date} 評估 · ${result.person.age}歲 · ${result.person.sex === 'male' ? '男性' : '女性'} · PRS 第${result.percentile.toFixed(2)}百分位`;
  byId('prs-model').textContent = `台灣 GBD2019 · INTERVENE EstBB · ${result.model_id}`;
  byId('prs-score').textContent = `PRS 第${result.percentile.toFixed(2)}百分位 · 分組 ${result.risk_to_80.group}%`;
  renderCurve(result);
  byId('prs-limitations').replaceChildren();
  for (const text of result.limitations || []) {
    const item=document.createElement('li');item.textContent=text;byId('prs-limitations').append(item);
  }
}

async function load() {
  if (!sampleId) return;
  const token=++generation;
  controller?.abort();controller=new AbortController();ready=false;
  displayedResult=null;
  byId('prs-result').hidden=true;
  byId('prs-score').textContent='';
  byId('prs-submit').disabled=true;
  byId('prs-message').textContent='讀取 PRS 與風險模型…';
  try {
    const data=await request();
    if(token!==generation)return;
    ready=data.available;
    byId('prs-form').hidden=!ready;
    byId('prs-submit').disabled=!ready;
    byId('prs-message').textContent=data.issues.join('；') || (data.stale ? '計分或個案資料已更新，請重新估算。' : '請確認評估資料後估算。');
    if(data.score)byId('prs-score').textContent=`PRS 第${data.score.percentile.toFixed(2)}百分位 · 分組 ${data.score.group}%`;
    byId('prs-result').hidden=true;
    if(data.latest) {
      const p=data.latest.result.person;
      byId('prs-age').value=p.age;byId('prs-sex').value=p.sex;
      byId('prs-date').value=p.assessment_date;byId('prs-history').value=p.prior_chd===false?'no':'';
      renderResult(data.latest);
      byId('prs-message').textContent='顯示上次評估；修改資料後請重新估算。';
    }
  } catch(error) {
    if(token===generation && error.name!=='AbortError')byId('prs-message').textContent=error.message;
  }
}

byId('prs-details')?.addEventListener('toggle',()=>{if(byId('prs-details').open)load();});
byId('prs-refresh')?.addEventListener('click',load);
byId('prs-form')?.addEventListener('input',()=>{
  // Do not display a previously calculated risk next to changed inputs.
  generation+=1;controller?.abort();
  byId('prs-result').hidden=true;byId('prs-submit').disabled=!ready;
  displayedResult = null;
  byId('prs-message').textContent='評估資料已變更，請重新估算。';
});
window.addEventListener('resize',()=>{
  if(displayedResult && !byId('prs-result').hidden)renderCurve(displayedResult);
});
byId('prs-form')?.addEventListener('submit',async event=>{
  event.preventDefault();if(!sampleId || !ready)return;
  byId('prs-result').hidden=true;
  if(byId('prs-history').value!=='no') {
    byId('prs-message').textContent='僅適用確認尚未罹患 CHD 的個案。';return;
  }
  const token=++generation;
  controller?.abort();controller=new AbortController();
  byId('prs-submit').disabled=true;byId('prs-message').textContent='計算並保存本次評估…';
  try {
    const result=await request('/assessments',{
      age:Number(byId('prs-age').value),sex:byId('prs-sex').value,
      assessment_date:byId('prs-date').value,prior_chd:false,
    });
    if(token!==generation)return;
    renderResult(result);byId('prs-message').textContent='已保存本次研究風險評估。';
  } catch(error) {
    if(token===generation && error.name!=='AbortError')byId('prs-message').textContent=error.message;
  } finally {if(token===generation)byId('prs-submit').disabled=!ready;}
});
