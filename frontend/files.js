// Shared file list. Native <details> keeps FASTQ expansion keyboard accessible.
const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
let invalidate = () => {};
export function resetFileInventory() { invalidate(); }

function fileEntry(file, icons) {
  return `<div class="file-entry"><a href="/api/files/download/${encodeURIComponent(file.id)}" title="${esc(file.path)}">${esc(file.name)}</a><button class="btn-copy" type="button" data-copy="${esc(file.path)}" title="複製路徑" aria-label="複製 ${esc(file.name)} 路徑">${icons.copy}</button></div>`;
}

export function renderFileRows(samples, icons) {
  return samples.map(sample => {
    const fastqs = sample.fastq || [];
    let fastq = fastqs.slice(0, 2).map(f => fileEntry(f, icons)).join('');
    if (fastqs.length > 2) fastq += `<details class="file-fastq-more"><summary><span class="when-closed">more</span><span class="when-open">less</span></summary>${fastqs.slice(2).map(f => fileEntry(f, icons)).join('')}</details>`;
    const columns = ['bam','vcf','other'].map(key => `<td>${(sample[key] || []).map(f => fileEntry(f, icons)).join('')}</td>`).join('');
    const folders = (sample.tertiary || []).map(folder => {
      const query = new URLSearchParams({root:folder.root, path:folder.path});
      return `<a class="clinvar-external-link file-tertiary-link" href="/api/files/browse?${esc(query)}" target="_blank" rel="noopener" title="瀏覽 ${esc(folder.name)}" aria-label="瀏覽 ${esc(folder.name)} 三級資料">${icons.external}</a>`;
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
  invalidate = () => {
    token++; clearTimeout(timer); controller?.abort();
    input.value = ''; offset = 0; total = 0;
    rows.innerHTML = ''; meta.textContent = ''; status.textContent = ''; count.textContent = '';
    modal.classList.add('hidden');
  };
  const busy = value => { refresh.disabled = value; previous.disabled = value || offset === 0; next.disabled = value || offset+limit >= total; };
  async function load(force = false) {
    const version = ++token;
    controller?.abort(); controller = new AbortController();
    const init = {signal:controller.signal};
    busy(true); status.textContent = force ? '正在更新檔案清單…' : '正在讀取…';
    try {
      if (force) await apiFetch('/files/refresh', {...init, method:'POST'});
      if (version !== token) return;
      const query = new URLSearchParams({q:input.value.trim(),offset,limit});
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
    showModal('file-inventory-modal'); offset = 0; void load();
  });
  input.addEventListener('input', () => {
    token++; controller?.abort(); clearTimeout(timer); offset = 0;
    timer = setTimeout(() => void load(), 250);
  });
  refresh.addEventListener('click', () => { offset = 0; void load(true); });
  previous.addEventListener('click', () => { offset = Math.max(0,offset-limit); void load(); });
  next.addEventListener('click', () => { offset += limit; void load(); });
  new MutationObserver(() => {
    if (modal.classList.contains('hidden')) { token++; controller?.abort(); clearTimeout(timer); }
  }).observe(modal, {attributes:true, attributeFilter:['class']});
}
