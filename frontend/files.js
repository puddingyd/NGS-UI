// Shared file list. Native <details> keeps FASTQ expansion keyboard accessible.
const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
let invalidate = () => {};
export function resetFileInventory() { invalidate(); }

function fileEntry(file, icons) {
  const href = file.tertiary_root
    ? `/api/files/tertiary-download?${new URLSearchParams({root:file.tertiary_root, path:file.relative_path})}`
    : `/api/files/download/${encodeURIComponent(file.id)}`;
  return `<div class="file-entry"><a href="${esc(href)}" title="${esc(file.path)}">${esc(file.name)}</a><button class="btn-copy" type="button" data-copy="${esc(file.path)}" title="複製路徑" aria-label="複製 ${esc(file.name)} 路徑">${icons.copy}</button></div>`;
}

export function renderFileRows(samples, icons) {
  return samples.map(sample => {
    const fastqs = sample.fastq || [];
    let fastq = fastqs.slice(0, 2).map(f => fileEntry(f, icons)).join('');
    if (fastqs.length > 2) fastq += `<details class="file-fastq-more"><summary><span class="when-closed">more</span><span class="when-open">less</span></summary>${fastqs.slice(2).map(f => fileEntry(f, icons)).join('')}</details>`;
    const columns = ['bam','vcf','other'].map(key => `<td>${(sample[key] || []).map(f => fileEntry(f, icons)).join('')}</td>`).join('');
    const folders = (sample.tertiary || []).map(folder => {
      const query = new URLSearchParams({root:folder.root, path:folder.path});
      const label = folder.display_date || folder.created_date || folder.name;
      const date = folder.date_at || folder.created_at;
      const dateLabel = folder.date_source === 'tertiary_annotation' ? '三級註解紀錄日期' : '建立時間';
      const title = `${folder.absolute_path || folder.name}${date ? `\n${dateLabel}：${date}` : ''}`;
      const copy = folder.absolute_path ? `<button class="btn-copy" type="button" data-copy="${esc(folder.absolute_path)}" title="複製三級資料夾路徑" aria-label="複製 ${esc(folder.name)} 三級資料夾路徑">${icons.copy}</button>` : '';
      return `<div class="file-entry file-tertiary-entry"><a class="file-tertiary-link" href="/api/files/browse?${esc(query)}" target="_blank" rel="noopener" title="${esc(title)}" aria-label="瀏覽 ${esc(folder.name)} 三級資料">${esc(label)}</a>${copy}</div>`;
    }).join('');
    return `<tr><th scope="row">${esc(sample.sample_id)}</th><td>${fastq}</td>${columns}<td class="file-tertiary-cell">${folders}</td></tr>`;
  }).join('');
}

export function setupFileInventory({apiFetch, showModal, icons}) {
  const modal = document.getElementById('file-inventory-modal');
  if (!modal) return;
  const input = document.getElementById('file-inventory-search');
  const rows = document.getElementById('file-inventory-rows');
  const status = document.getElementById('file-inventory-status');
  const meta = document.getElementById('file-inventory-meta');
  const refresh = document.getElementById('file-inventory-refresh');
  const previous = document.getElementById('file-inventory-prev');
  const next = document.getElementById('file-inventory-next');
  const count = document.getElementById('file-inventory-count');
  let token = 0, offset = 0, total = 0, timer, controller;
  const limit = 50;
  const busy = value => { refresh.disabled = value; previous.disabled = value || offset === 0; next.disabled = value || offset+limit >= total; };
  const clearResults = () => {
    offset = 0; total = 0; rows.innerHTML = ''; count.textContent = '';
    status.textContent = input.value.trim() ? '' : '請輸入 sample ID'; busy(false);
  };
  invalidate = () => {
    token++; clearTimeout(timer); controller?.abort();
    input.value = ''; offset = 0; total = 0;
    rows.innerHTML = ''; meta.textContent = ''; status.textContent = ''; count.textContent = '';
    modal.classList.add('hidden');
  };
  async function load(force = false) {
    const version = ++token;
    controller?.abort(); controller = new AbortController();
    const init = {signal:controller.signal};
    const term = input.value.trim();
    if (!term && !force) { clearResults(); return; }
    busy(true); status.textContent = force ? '正在更新檔案清單…' : '正在讀取…';
    try {
      if (force) {
        const updated = await apiFetch('/files/refresh', {...init, method:'POST'});
        if (version === token && updated?.updated_at) meta.textContent = `更新時間：${updated.updated_at}`;
      }
      if (version !== token) return;
      if (!term) { clearResults(); return; }
      const query = new URLSearchParams({q:term,offset,limit});
      const data = await apiFetch(`/files?${query}`, init);
      if (version !== token || modal.classList.contains('hidden')) return;
      if (!data) throw new Error('無法讀取檔案清單');
      total = data.total;
      rows.innerHTML = renderFileRows(data.samples || [], icons);
      meta.textContent = data.meta?.updated_at ? `更新時間：${data.meta.updated_at}` : '';
      count.textContent = total ? `${offset+1}–${Math.min(offset+limit,total)} / ${total}` : '0';
      status.textContent = total ? '' : '沒有符合的 sample';
    } catch (error) {
      if (version !== token || error.name === 'AbortError') return;
      status.textContent = '讀取失敗，請稍後重試。';
    } finally {
      if (version === token) busy(false);
    }
  }
  document.getElementById('btn-file-inventory')?.addEventListener('click', () => {
    token++; controller?.abort(); clearTimeout(timer);
    input.value = ''; meta.textContent = ''; clearResults();
    showModal('file-inventory-modal');
  });
  input.addEventListener('input', () => {
    token++; controller?.abort(); clearTimeout(timer); clearResults();
    if (!input.value.trim()) return;
    timer = setTimeout(() => void load(), 250);
  });
  refresh.addEventListener('click', () => { offset = 0; void load(true); });
  previous.addEventListener('click', () => { offset = Math.max(0,offset-limit); void load(); });
  next.addEventListener('click', () => { offset += limit; void load(); });
  new MutationObserver(() => {
    if (modal.classList.contains('hidden')) { token++; controller?.abort(); clearTimeout(timer); }
  }).observe(modal, {attributes:true, attributeFilter:['class']});
}
