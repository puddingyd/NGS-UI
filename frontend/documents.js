(function () {
  "use strict";

  const API = "/api/documents";
  const ACCEPT = ".pdf,.jpg,.jpeg,.png,.tif,.tiff,application/pdf,image/jpeg,image/png,image/tiff";
  const MRN_RE = /^[A-Za-z0-9_-]{1,32}$/;
  const PREVIEW_ZOOM_MIN = 0.5;
  const PREVIEW_ZOOM_MAX = 4;
  const PREVIEW_ZOOM_STEP = 0.25;
  const PREVIEW_PAN_STEP = 80;
  let options = {};
  let context = null;
  let pending = [];
  let previewUrl = "";
  let previewDocument = null;
  let previewPage = 0;
  let previewDocuments = [];
  let previewRequest = 0;
  let previewController = null;
  let previewZoom = 1;
  let previewPanX = 0;
  let previewPanY = 0;
  let previewDrag = null;

  function esc(value) {
    return String(value == null ? "" : value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  function formatBytes(value) {
    let size = Number(value) || 0;
    if (size < 1024) return `${size} B`;
    const units = ["KB", "MB", "GB", "TB"];
    let unit = 0;
    size /= 1024;
    while (size >= 1024 && unit < units.length - 1) {
      size /= 1024;
      unit += 1;
    }
    return `${size.toFixed(size >= 10 ? 1 : 2)} ${units[unit]}`;
  }

  function formatTime(value) {
    if (!value) return "—";
    try { return new Date(value).toLocaleString(); }
    catch (_error) { return String(value); }
  }

  function splitName(name, fallbackExt) {
    const clean = String(name || "Document").replace(/[\\/]/g, "_").trim();
    const match = clean.match(/^(.*?)(\.[^.]+)$/);
    if (match) return { base: match[1] || "Document", ext: match[2].toLowerCase() };
    return { base: clean || "Document", ext: fallbackExt || "" };
  }

  function fileExtension(file) {
    const name = String(file?.name || "").toLowerCase();
    for (const ext of [".jpeg", ".tiff", ".jpg", ".png", ".tif", ".pdf"]) {
      if (name.endsWith(ext)) return ext;
    }
    const type = String(file?.type || "").toLowerCase();
    if (type === "image/png") return ".png";
    if (type === "image/jpeg") return ".jpg";
    if (type === "image/tiff") return ".tiff";
    if (type === "application/pdf") return ".pdf";
    return "";
  }

  function supported(file) {
    return !!fileExtension(file);
  }

  function nativePreviewable(file) {
    const ext = fileExtension(file);
    return ext === ".png" || ext === ".jpg" || ext === ".jpeg";
  }

  function screenshotName() {
    const date = new Date();
    const pad = value => String(value).padStart(2, "0");
    return `Screenshot_${date.getFullYear()}${pad(date.getMonth() + 1)}${pad(date.getDate())}_${pad(date.getHours())}${pad(date.getMinutes())}${pad(date.getSeconds())}.png`;
  }

  function injectShell() {
    if (document.getElementById("patient-documents-modal")) return;
    const shell = document.createElement("div");
    shell.innerHTML = `
      <div id="patient-documents-modal" class="pdoc-overlay" hidden>
        <div class="pdoc-card" role="dialog" aria-modal="true" aria-labelledby="pdoc-title">
          <div class="pdoc-head">
            <div>
              <h2 id="pdoc-title">Documents</h2>
              <div id="pdoc-patient" class="pdoc-muted"></div>
            </div>
            <button type="button" class="pdoc-icon-btn" data-pdoc-close aria-label="關閉">&times;</button>
          </div>

          <section id="pdoc-login" class="pdoc-login" hidden>
            <h3>登入後管理病歷文件</h3>
            <form id="pdoc-login-form">
              <label>帳號<input id="pdoc-login-user" type="text" autocomplete="username" required></label>
              <label>密碼<input id="pdoc-login-password" type="password" autocomplete="current-password" required></label>
              <div id="pdoc-login-error" class="pdoc-error" aria-live="polite"></div>
              <button type="submit" class="btn btn-primary">登入</button>
            </form>
          </section>

          <div id="pdoc-content" hidden>
            <section class="pdoc-section">
              <div class="pdoc-section-head">
                <h3>新增文件</h3>
                <button id="pdoc-pick" type="button" class="btn btn-secondary">選擇檔案</button>
                <input id="pdoc-file-input" type="file" accept="${ACCEPT}" multiple hidden>
              </div>
              <div id="pdoc-paste-zone" class="pdoc-paste-zone" tabindex="0">
                <strong>拖曳檔案或貼上截圖到這裡</strong>
                <span>可拖入 PDF / JPG / PNG / TIFF，或點一下後按 Ctrl+V / ⌘V 貼圖</span>
              </div>
              <div id="pdoc-pending" class="pdoc-pending"></div>
              <div class="pdoc-upload-actions">
                <span id="pdoc-upload-status" class="pdoc-muted" aria-live="polite"></span>
                <button id="pdoc-upload-all" type="button" class="btn btn-primary" hidden>儲存全部</button>
              </div>
            </section>

            <section class="pdoc-section pdoc-existing-section">
              <div class="pdoc-section-head">
                <h3>已儲存文件</h3>
                <a id="pdoc-download-all" class="btn btn-secondary" target="_blank" rel="noopener" hidden>下載全部 ZIP</a>
                <button id="pdoc-refresh" type="button" class="btn btn-secondary">重新整理</button>
              </div>
              <div id="pdoc-list-status" class="pdoc-muted" aria-live="polite"></div>
              <div id="pdoc-list" class="pdoc-list"></div>
            </section>
          </div>
        </div>
      </div>

      <div id="patient-document-preview" class="pdoc-overlay pdoc-preview-overlay" hidden>
        <div class="pdoc-preview-card" role="dialog" aria-modal="true" aria-labelledby="pdoc-preview-title">
          <div class="pdoc-head">
            <div>
              <h2 id="pdoc-preview-title">圖片預覽</h2>
              <div id="pdoc-preview-page-label" class="pdoc-muted"></div>
            </div>
            <button type="button" class="pdoc-icon-btn" data-pdoc-preview-close aria-label="關閉">&times;</button>
          </div>
          <div class="pdoc-preview-stage">
            <button type="button" class="pdoc-file-nav" data-pdoc-file-prev aria-label="上一個檔案" title="上一個檔案（←）">&#10094;</button>
            <div id="pdoc-preview-body" class="pdoc-preview-body" aria-live="polite"><span class="pdoc-muted">載入中…</span></div>
            <button type="button" class="pdoc-file-nav" data-pdoc-file-next aria-label="下一個檔案" title="下一個檔案（→）">&#10095;</button>
          </div>
          <div id="pdoc-view-controls" class="pdoc-view-controls" role="toolbar" aria-label="圖片檢視工具">
            <div class="pdoc-tool-group" role="group" aria-label="縮放">
              <button type="button" class="pdoc-tool-btn" data-pdoc-zoom-out aria-label="縮小" title="縮小（−）">&#8722;</button>
              <output id="pdoc-zoom-label" class="pdoc-zoom-label" aria-live="polite">100%</output>
              <button type="button" class="pdoc-tool-btn" data-pdoc-zoom-in aria-label="放大" title="放大（＋）">&#43;</button>
              <button type="button" class="pdoc-tool-btn pdoc-reset-btn" data-pdoc-zoom-reset title="回到符合視窗的大小（0）">重設</button>
            </div>
            <div class="pdoc-tool-group pdoc-pan-tools" role="group" aria-label="移動圖片">
              <span class="pdoc-tool-label">移動</span>
              <button type="button" class="pdoc-tool-btn" data-pdoc-pan-left aria-label="圖片向左移" title="圖片向左移（Shift＋←）">&#8592;</button>
              <button type="button" class="pdoc-tool-btn" data-pdoc-pan-up aria-label="圖片向上移" title="圖片向上移（Shift＋↑）">&#8593;</button>
              <button type="button" class="pdoc-tool-btn" data-pdoc-pan-down aria-label="圖片向下移" title="圖片向下移（Shift＋↓）">&#8595;</button>
              <button type="button" class="pdoc-tool-btn" data-pdoc-pan-right aria-label="圖片向右移" title="圖片向右移（Shift＋→）">&#8594;</button>
            </div>
            <span class="pdoc-view-hint">放大後可拖曳；Ctrl/⌘＋滾輪縮放</span>
          </div>
          <div id="pdoc-preview-controls" class="pdoc-preview-controls" hidden>
            <button type="button" class="btn btn-secondary" data-pdoc-preview-prev>上一頁</button>
            <button type="button" class="btn btn-secondary" data-pdoc-preview-next>下一頁</button>
          </div>
        </div>
      </div>`;
    while (shell.firstElementChild) document.body.appendChild(shell.firstElementChild);
    wireShell();
  }

  function setStatus(id, message, error) {
    const element = document.getElementById(id);
    if (!element) return;
    element.textContent = message || "";
    element.classList.toggle("pdoc-error", !!error);
  }

  async function jsonRequest(url, init) {
    const response = await fetch(url, { credentials: "same-origin", ...(init || {}) });
    const body = await response.json().catch(() => ({}));
    if (response.status === 401) {
      showLogin();
      throw new Error("請先登入");
    }
    if (!response.ok) throw new Error(body.detail || `${response.status} ${response.statusText}`);
    return body;
  }

  function showLogin(message) {
    document.getElementById("pdoc-login").hidden = false;
    document.getElementById("pdoc-content").hidden = true;
    document.getElementById("pdoc-login-error").textContent = message || "";
    setTimeout(() => document.getElementById("pdoc-login-user")?.focus(), 0);
  }

  function showContent() {
    document.getElementById("pdoc-login").hidden = true;
    document.getElementById("pdoc-content").hidden = false;
  }

  async function ensureAuthenticated() {
    const response = await fetch("/api/auth/me", { credentials: "same-origin" });
    if (!response.ok) {
      showLogin();
      return false;
    }
    showContent();
    return true;
  }

  function addFiles(files, forcedName) {
    const rejected = [];
    Array.from(files || []).forEach(file => {
      if (!supported(file)) {
        rejected.push(file.name || "未知檔案");
        return;
      }
      const ext = fileExtension(file);
      const parts = splitName(forcedName || file.name || `Document${ext}`, ext);
      const record = {
        key: `${Date.now()}-${Math.random()}`,
        file,
        base: parts.base,
        ext,
        previewUrl: nativePreviewable(file) ? URL.createObjectURL(file) : "",
      };
      pending.push(record);
    });
    if (rejected.length) alert(`不支援以下檔案：\n${rejected.join("\n")}\n\n只支援 PDF、JPG、PNG、TIF、TIFF。`);
    renderPending();
  }

  function renderPending() {
    const host = document.getElementById("pdoc-pending");
    const upload = document.getElementById("pdoc-upload-all");
    upload.hidden = !pending.length;
    if (!pending.length) {
      host.innerHTML = "";
      return;
    }
    host.innerHTML = pending.map(item => {
      const thumb = item.previewUrl
        ? `<img src="${esc(item.previewUrl)}" alt="待上傳圖片預覽">`
        : `<span class="pdoc-file-kind">${item.ext === ".pdf" ? "PDF" : "TIFF"}</span>`;
      return `<div class="pdoc-pending-row" data-pending-key="${esc(item.key)}">
        <div class="pdoc-pending-thumb">${thumb}</div>
        <label class="pdoc-name-field">
          <span>檔名</span>
          <span class="pdoc-name-control"><input type="text" value="${esc(item.base)}" maxlength="180"><b>${esc(item.ext)}</b></span>
        </label>
        <span class="pdoc-size">${esc(formatBytes(item.file.size))}</span>
        <button type="button" class="btn btn-secondary" data-pdoc-remove-pending>移除</button>
      </div>`;
    }).join("");
  }

  function removePending(key) {
    const index = pending.findIndex(item => item.key === key);
    if (index < 0) return;
    if (pending[index].previewUrl) URL.revokeObjectURL(pending[index].previewUrl);
    pending.splice(index, 1);
    renderPending();
  }

  function clearPending() {
    pending.forEach(item => { if (item.previewUrl) URL.revokeObjectURL(item.previewUrl); });
    pending = [];
    renderPending();
  }

  async function uploadAll() {
    if (!pending.length || !context) return;
    const button = document.getElementById("pdoc-upload-all");
    const rows = Array.from(document.querySelectorAll(".pdoc-pending-row"));
    const names = new Map(rows.map(row => [
      row.dataset.pendingKey,
      String(row.querySelector("input")?.value || "").trim(),
    ]));
    if (Array.from(names.values()).some(name => !name)) {
      setStatus("pdoc-upload-status", "檔名不可為空", true);
      return;
    }
    button.disabled = true;
    const completed = [];
    let failure = "";
    for (let index = 0; index < pending.length; index += 1) {
      const item = pending[index];
      setStatus("pdoc-upload-status", `上傳中 ${index + 1}/${pending.length}：${names.get(item.key)}${item.ext}`);
      const data = new FormData();
      data.append("mrn", context.mrn);
      data.append("source_sample_id", context.sourceSampleId || "");
      data.append("display_name", `${names.get(item.key)}${item.ext}`);
      data.append("file", item.file, item.file.name || `${names.get(item.key)}${item.ext}`);
      try {
        await jsonRequest(API, { method: "POST", body: data });
        completed.push(item.key);
      } catch (error) {
        failure = `${names.get(item.key)}${item.ext}：${error.message || error}`;
        break;
      }
    }
    completed.forEach(removePending);
    button.disabled = false;
    if (failure) setStatus("pdoc-upload-status", `上傳停止：${failure}`, true);
    else setStatus("pdoc-upload-status", "文件已儲存");
    await loadList();
  }

  function existingExtension(documentInfo) {
    const parts = splitName(documentInfo.display_name, "");
    return parts;
  }

  function renderList(rows) {
    previewDocuments = rows.filter(row => row.previewable);
    const host = document.getElementById("pdoc-list");
    const downloadAll = document.getElementById("pdoc-download-all");
    if (downloadAll) {
      downloadAll.hidden = !rows.length || !context;
      downloadAll.href = context
        ? `${API}/archive.zip?mrn=${encodeURIComponent(context.mrn)}`
        : "";
    }
    if (!rows.length) {
      host.innerHTML = `<div class="pdoc-empty">（尚無文件）</div>`;
      return;
    }
    host.innerHTML = rows.map(row => {
      const parts = existingExtension(row);
      const preview = row.previewable
        ? `<button type="button" class="btn btn-secondary" data-pdoc-preview>預覽</button>`
        : "";
      return `<article class="pdoc-row" data-document-id="${esc(row.id)}" data-document='${esc(JSON.stringify(row))}'>
        <div class="pdoc-row-main">
          <div class="pdoc-existing-name">
            <input type="text" value="${esc(parts.base)}" maxlength="180" disabled>
            <b>${esc(parts.ext)}</b>
          </div>
          <div class="pdoc-meta">
            ${esc(formatTime(row.created_at))} · ${esc(row.created_by_username || "—")} · ${esc(row.file_format)} · ${esc(formatBytes(row.size_bytes))}
            ${Number(row.image_pages || 1) > 1 ? ` · ${esc(row.image_pages)} 頁` : ""}
          </div>
        </div>
        <div class="pdoc-row-actions">
          ${preview}
          <a class="btn btn-secondary" href="${API}/${encodeURIComponent(row.id)}/download" target="_blank" rel="noopener">下載</a>
          <button type="button" class="btn btn-secondary" data-pdoc-rename>修改檔名</button>
          <button type="button" class="btn pdoc-delete-btn" data-pdoc-delete>刪除</button>
        </div>
      </article>`;
    }).join("");
  }

  async function loadList() {
    if (!context) return;
    const currentContext = context;
    setStatus("pdoc-list-status", "載入中…");
    try {
      const rows = await jsonRequest(`${API}?mrn=${encodeURIComponent(context.mrn)}`);
      if (context !== currentContext) return;
      renderList(Array.isArray(rows) ? rows : []);
      setStatus("pdoc-list-status", "");
    } catch (error) {
      if (context !== currentContext) return;
      setStatus("pdoc-list-status", `載入失敗：${error.message || error}`, true);
    }
  }

  async function toggleRename(rowElement) {
    const input = rowElement.querySelector(".pdoc-existing-name input");
    const button = rowElement.querySelector("[data-pdoc-rename]");
    if (input.disabled) {
      input.disabled = false;
      rowElement.classList.add("is-editing");
      button.textContent = "儲存檔名";
      input.focus();
      input.select();
      return;
    }
    const documentInfo = JSON.parse(rowElement.dataset.document || "{}");
    const ext = rowElement.querySelector(".pdoc-existing-name b")?.textContent || "";
    const base = input.value.trim();
    if (!base) {
      setStatus("pdoc-list-status", "檔名不可為空", true);
      return;
    }
    button.disabled = true;
    try {
      await jsonRequest(`${API}/${encodeURIComponent(documentInfo.id)}`, {
        method: "PATCH",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ display_name: `${base}${ext}` }),
      });
      await loadList();
    } catch (error) {
      setStatus("pdoc-list-status", `修改失敗：${error.message || error}`, true);
      button.disabled = false;
    }
  }

  async function deleteDocument(rowElement) {
    const info = JSON.parse(rowElement.dataset.document || "{}");
    if (!confirm(`確定刪除「${info.display_name || "這個文件"}」？`)) return;
    try {
      await jsonRequest(`${API}/${encodeURIComponent(info.id)}`, { method: "DELETE" });
      await loadList();
    } catch (error) {
      setStatus("pdoc-list-status", `刪除失敗：${error.message || error}`, true);
    }
  }

  function releasePreviewUrl() {
    if (previewUrl) URL.revokeObjectURL(previewUrl);
    previewUrl = "";
  }

  function clamp(value, minimum, maximum) {
    return Math.min(maximum, Math.max(minimum, value));
  }

  function previewPanBounds() {
    const body = document.getElementById("pdoc-preview-body");
    const image = body?.querySelector("img");
    if (!body || !image || previewZoom <= 1) return { x: 0, y: 0 };
    const availableWidth = Math.max(0, body.clientWidth - 28);
    const availableHeight = Math.max(0, body.clientHeight - 28);
    return {
      x: Math.max(0, (image.offsetWidth * previewZoom - availableWidth) / 2),
      y: Math.max(0, (image.offsetHeight * previewZoom - availableHeight) / 2),
    };
  }

  function renderPreviewViewportControls() {
    const bounds = previewPanBounds();
    const hasImage = !!document.querySelector("#pdoc-preview-body img");
    const label = document.getElementById("pdoc-zoom-label");
    if (label) label.value = `${Math.round(previewZoom * 100)}%`;
    const setDisabled = (selector, disabled) => {
      const button = document.querySelector(selector);
      if (button) button.disabled = disabled;
    };
    setDisabled("[data-pdoc-zoom-out]", !hasImage || previewZoom <= PREVIEW_ZOOM_MIN);
    setDisabled("[data-pdoc-zoom-in]", !hasImage || previewZoom >= PREVIEW_ZOOM_MAX);
    setDisabled("[data-pdoc-zoom-reset]", !hasImage || (
      previewZoom === 1 && previewPanX === 0 && previewPanY === 0
    ));
    setDisabled("[data-pdoc-pan-left]", !hasImage || bounds.x === 0 || previewPanX <= -bounds.x);
    setDisabled("[data-pdoc-pan-right]", !hasImage || bounds.x === 0 || previewPanX >= bounds.x);
    setDisabled("[data-pdoc-pan-up]", !hasImage || bounds.y === 0 || previewPanY <= -bounds.y);
    setDisabled("[data-pdoc-pan-down]", !hasImage || bounds.y === 0 || previewPanY >= bounds.y);
  }

  function applyPreviewViewport() {
    const body = document.getElementById("pdoc-preview-body");
    const image = body?.querySelector("img");
    if (!body || !image) {
      renderPreviewViewportControls();
      return;
    }
    const bounds = previewPanBounds();
    previewPanX = clamp(previewPanX, -bounds.x, bounds.x);
    previewPanY = clamp(previewPanY, -bounds.y, bounds.y);
    image.style.transform = `translate3d(${previewPanX}px, ${previewPanY}px, 0) scale(${previewZoom})`;
    body.classList.toggle("is-zoomed", previewZoom > 1);
    renderPreviewViewportControls();
  }

  function resetPreviewViewport() {
    previewZoom = 1;
    previewPanX = 0;
    previewPanY = 0;
    previewDrag = null;
    const body = document.getElementById("pdoc-preview-body");
    body?.classList.remove("is-zoomed", "is-dragging");
    applyPreviewViewport();
  }

  function setPreviewZoom(value) {
    previewZoom = Math.round(clamp(value, PREVIEW_ZOOM_MIN, PREVIEW_ZOOM_MAX) * 100) / 100;
    if (previewZoom <= 1) {
      previewPanX = 0;
      previewPanY = 0;
    }
    applyPreviewViewport();
  }

  function panPreview(deltaX, deltaY) {
    if (previewZoom <= 1) return;
    previewPanX += deltaX;
    previewPanY += deltaY;
    applyPreviewViewport();
  }

  async function loadPreviewPage(page) {
    if (!previewDocument) return;
    const pages = Number(previewDocument.image_pages || 1);
    if (!Number.isInteger(page) || page < 0 || page >= pages) return;
    const request = ++previewRequest;
    previewController?.abort();
    previewController = new AbortController();
    const currentDocument = previewDocument;
    previewPage = page;
    resetPreviewViewport();
    renderPreviewControls();
    const body = document.getElementById("pdoc-preview-body");
    body.innerHTML = `<span class="pdoc-muted">載入中…</span>`;
    renderPreviewViewportControls();
    releasePreviewUrl();
    try {
      const response = await fetch(
        `${API}/${encodeURIComponent(currentDocument.id)}/preview?page=${page}`,
        { credentials: "same-origin", signal: previewController.signal },
      );
      if (request !== previewRequest) return;
      if (response.status === 401) {
        closePreview();
        showLogin();
        return;
      }
      if (!response.ok) {
        const error = await response.json().catch(() => ({}));
        throw new Error(error.detail || `預覽失敗 (${response.status})`);
      }
      const blob = await response.blob();
      if (request !== previewRequest) return;
      previewUrl = URL.createObjectURL(blob);
      body.innerHTML = `<img src="${esc(previewUrl)}" alt="${esc(currentDocument.display_name)}" draggable="false">`;
      const image = body.querySelector("img");
      image?.addEventListener("load", applyPreviewViewport, { once: true });
      if (image?.complete) applyPreviewViewport();
    } catch (error) {
      if (request !== previewRequest || error.name === "AbortError") return;
      body.innerHTML = `<div class="pdoc-error">${esc(error.message || error)}</div>`;
      renderPreviewViewportControls();
    }
  }

  function renderPreviewControls() {
    const pages = Number(previewDocument?.image_pages || 1);
    const controls = document.getElementById("pdoc-preview-controls");
    const label = document.getElementById("pdoc-preview-page-label");
    const index = previewDocuments.findIndex(row => row.id === previewDocument?.id);
    controls.hidden = pages <= 1;
    label.textContent = [
      index >= 0 ? `檔案 ${index + 1} / ${previewDocuments.length}` : "",
      pages > 1 ? `第 ${previewPage + 1} / ${pages} 頁` : "",
    ].filter(Boolean).join(" · ");
    document.querySelector("[data-pdoc-file-prev]").disabled = index <= 0;
    document.querySelector("[data-pdoc-file-next]").disabled = index < 0 || index >= previewDocuments.length - 1;
    controls.querySelector("[data-pdoc-preview-prev]").disabled = previewPage <= 0;
    controls.querySelector("[data-pdoc-preview-next]").disabled = previewPage >= pages - 1;
  }

  function openPreview(rowElement) {
    showPreviewDocument(JSON.parse(rowElement.dataset.document || "{}"));
    document.querySelector("[data-pdoc-preview-close]")?.focus();
  }

  function showPreviewDocument(documentInfo) {
    previewDocument = documentInfo;
    previewPage = 0;
    document.getElementById("pdoc-preview-title").textContent = previewDocument.display_name || "圖片預覽";
    document.getElementById("patient-document-preview").hidden = false;
    loadPreviewPage(0);
  }

  function movePreviewDocument(offset) {
    if (!previewDocument) return;
    const index = previewDocuments.findIndex(row => row.id === previewDocument.id);
    if (index < 0 || !previewDocuments[index + offset]) return;
    showPreviewDocument(previewDocuments[index + offset]);
  }

  function closePreview() {
    previewRequest += 1;
    previewController?.abort();
    previewController = null;
    const previousId = previewDocument?.id;
    releasePreviewUrl();
    previewDocument = null;
    resetPreviewViewport();
    document.getElementById("patient-document-preview").hidden = true;
    document.getElementById("pdoc-preview-body").innerHTML = "";
    if (previousId) {
      const row = Array.from(document.querySelectorAll(".pdoc-row"))
        .find(item => item.dataset.documentId === previousId);
      row?.querySelector("[data-pdoc-preview]")?.focus();
    }
  }

  async function login(event) {
    event.preventDefault();
    const username = document.getElementById("pdoc-login-user").value.trim();
    const password = document.getElementById("pdoc-login-password").value;
    const errorElement = document.getElementById("pdoc-login-error");
    errorElement.textContent = "";
    try {
      const response = await fetch("/api/auth/login", {
        method: "POST",
        credentials: "same-origin",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ username, password }),
      });
      const body = await response.json().catch(() => ({}));
      if (!response.ok) throw new Error(response.status === 401 ? "帳號或密碼錯誤" : (body.detail || "登入失敗"));
      document.getElementById("pdoc-login-password").value = "";
      showContent();
      if (typeof options.onLogin === "function") options.onLogin(body);
      await loadList();
    } catch (error) {
      errorElement.textContent = error.message || String(error);
    }
  }

  function wireShell() {
    document.querySelector("[data-pdoc-close]")?.addEventListener("click", close);
    document.querySelector("[data-pdoc-preview-close]")?.addEventListener("click", closePreview);
    document.getElementById("pdoc-login-form")?.addEventListener("submit", login);
    document.getElementById("pdoc-pick")?.addEventListener("click", () => {
      const input = document.getElementById("pdoc-file-input");
      input.value = "";
      input.click();
    });
    document.getElementById("pdoc-file-input")?.addEventListener("change", event => addFiles(event.target.files));
    document.getElementById("pdoc-upload-all")?.addEventListener("click", uploadAll);
    document.getElementById("pdoc-refresh")?.addEventListener("click", loadList);

    const pasteZone = document.getElementById("pdoc-paste-zone");
    pasteZone?.addEventListener("click", () => pasteZone.focus());
    for (const eventName of ["dragenter", "dragover"]) {
      pasteZone?.addEventListener(eventName, event => {
        event.preventDefault();
        event.stopPropagation();
        if (event.dataTransfer) event.dataTransfer.dropEffect = "copy";
        pasteZone.classList.add("is-dragover");
      });
    }
    pasteZone?.addEventListener("dragleave", event => {
      event.preventDefault();
      if (!pasteZone.contains(event.relatedTarget)) {
        pasteZone.classList.remove("is-dragover");
      }
    });
    pasteZone?.addEventListener("drop", event => {
      event.preventDefault();
      event.stopPropagation();
      pasteZone.classList.remove("is-dragover");
      const files = event.dataTransfer?.files || [];
      if (!files.length) {
        setStatus("pdoc-upload-status", "沒有偵測到可上傳的檔案", true);
        return;
      }
      addFiles(files);
      setStatus("pdoc-upload-status", `已加入 ${files.length} 個拖曳檔案，請確認檔名後儲存`);
    });
    pasteZone?.addEventListener("paste", event => {
      let images = Array.from(event.clipboardData?.files || []).filter(file => String(file.type || "").startsWith("image/"));
      if (!images.length) {
        images = Array.from(event.clipboardData?.items || [])
          .filter(item => item.kind === "file" && String(item.type || "").startsWith("image/"))
          .map(item => item.getAsFile())
          .filter(Boolean);
      }
      if (!images.length) {
        setStatus("pdoc-upload-status", "剪貼簿裡沒有圖片", true);
        return;
      }
      event.preventDefault();
      images.forEach((file, index) => addFiles([file], index ? screenshotName().replace(".png", `_${index + 1}.png`) : screenshotName()));
      setStatus("pdoc-upload-status", `已貼上 ${images.length} 張圖片，請確認檔名後儲存`);
    });

    document.getElementById("pdoc-pending")?.addEventListener("click", event => {
      const button = event.target.closest("[data-pdoc-remove-pending]");
      if (!button) return;
      removePending(button.closest(".pdoc-pending-row")?.dataset.pendingKey || "");
    });
    document.getElementById("pdoc-list")?.addEventListener("click", event => {
      const row = event.target.closest(".pdoc-row");
      if (!row) return;
      if (event.target.closest("[data-pdoc-preview]")) openPreview(row);
      else if (event.target.closest("[data-pdoc-rename]")) toggleRename(row);
      else if (event.target.closest("[data-pdoc-delete]")) deleteDocument(row);
    });
    document.querySelector("[data-pdoc-preview-prev]")?.addEventListener("click", () => loadPreviewPage(previewPage - 1));
    document.querySelector("[data-pdoc-preview-next]")?.addEventListener("click", () => loadPreviewPage(previewPage + 1));
    document.querySelector("[data-pdoc-file-prev]")?.addEventListener("click", () => movePreviewDocument(-1));
    document.querySelector("[data-pdoc-file-next]")?.addEventListener("click", () => movePreviewDocument(1));
    document.querySelector("[data-pdoc-zoom-out]")?.addEventListener("click", () => setPreviewZoom(previewZoom - PREVIEW_ZOOM_STEP));
    document.querySelector("[data-pdoc-zoom-in]")?.addEventListener("click", () => setPreviewZoom(previewZoom + PREVIEW_ZOOM_STEP));
    document.querySelector("[data-pdoc-zoom-reset]")?.addEventListener("click", resetPreviewViewport);
    document.querySelector("[data-pdoc-pan-left]")?.addEventListener("click", () => panPreview(-PREVIEW_PAN_STEP, 0));
    document.querySelector("[data-pdoc-pan-right]")?.addEventListener("click", () => panPreview(PREVIEW_PAN_STEP, 0));
    document.querySelector("[data-pdoc-pan-up]")?.addEventListener("click", () => panPreview(0, -PREVIEW_PAN_STEP));
    document.querySelector("[data-pdoc-pan-down]")?.addEventListener("click", () => panPreview(0, PREVIEW_PAN_STEP));
    const previewBody = document.getElementById("pdoc-preview-body");
    previewBody?.addEventListener("pointerdown", event => {
      if (previewZoom <= 1 || !event.target.closest("img")) return;
      event.preventDefault();
      previewDrag = {
        pointerId: event.pointerId,
        startX: event.clientX,
        startY: event.clientY,
        panX: previewPanX,
        panY: previewPanY,
      };
      previewBody.setPointerCapture(event.pointerId);
      previewBody.classList.add("is-dragging");
    });
    previewBody?.addEventListener("pointermove", event => {
      if (!previewDrag || previewDrag.pointerId !== event.pointerId) return;
      previewPanX = previewDrag.panX + event.clientX - previewDrag.startX;
      previewPanY = previewDrag.panY + event.clientY - previewDrag.startY;
      applyPreviewViewport();
    });
    const stopPreviewDrag = event => {
      if (!previewDrag || previewDrag.pointerId !== event.pointerId) return;
      if (previewBody.hasPointerCapture(event.pointerId)) previewBody.releasePointerCapture(event.pointerId);
      previewDrag = null;
      previewBody.classList.remove("is-dragging");
    };
    previewBody?.addEventListener("pointerup", stopPreviewDrag);
    previewBody?.addEventListener("pointercancel", stopPreviewDrag);
    previewBody?.addEventListener("wheel", event => {
      if ((!event.ctrlKey && !event.metaKey) || !previewBody.querySelector("img")) return;
      event.preventDefault();
      setPreviewZoom(previewZoom + (event.deltaY < 0 ? PREVIEW_ZOOM_STEP : -PREVIEW_ZOOM_STEP));
    }, { passive: false });
    previewBody?.addEventListener("dblclick", event => {
      if (!event.target.closest("img")) return;
      setPreviewZoom(previewZoom === 1 ? 2 : 1);
    });
    document.getElementById("patient-documents-modal")?.addEventListener("click", event => {
      if (event.target.id === "patient-documents-modal") close();
    });
    document.getElementById("patient-document-preview")?.addEventListener("click", event => {
      if (event.target.id === "patient-document-preview") closePreview();
    });
    document.addEventListener("keydown", event => {
      const previewIsOpen = !document.getElementById("patient-document-preview")?.hidden;
      const isTyping = !!event.target.closest("input, textarea, select, [contenteditable]");
      if (previewIsOpen && !isTyping && !event.ctrlKey && !event.metaKey && !event.altKey) {
        if (event.shiftKey && ["ArrowLeft", "ArrowRight", "ArrowUp", "ArrowDown"].includes(event.key)) {
          event.preventDefault();
          const movement = {
            ArrowLeft: [-PREVIEW_PAN_STEP, 0], ArrowRight: [PREVIEW_PAN_STEP, 0],
            ArrowUp: [0, -PREVIEW_PAN_STEP], ArrowDown: [0, PREVIEW_PAN_STEP],
          }[event.key];
          panPreview(...movement);
          return;
        }
        if (["+", "="].includes(event.key)) {
          event.preventDefault();
          setPreviewZoom(previewZoom + PREVIEW_ZOOM_STEP);
          return;
        }
        if (["-", "_"].includes(event.key)) {
          event.preventDefault();
          setPreviewZoom(previewZoom - PREVIEW_ZOOM_STEP);
          return;
        }
        if (event.key === "0") {
          event.preventDefault();
          resetPreviewViewport();
          return;
        }
      }
      if (!document.getElementById("patient-document-preview")?.hidden
          && ["ArrowLeft", "ArrowRight"].includes(event.key)
          && !event.altKey && !event.ctrlKey && !event.metaKey && !event.shiftKey
          && !event.target.closest("input, textarea, select, [contenteditable]")) {
        event.preventDefault();
        movePreviewDocument(event.key === "ArrowLeft" ? -1 : 1);
        return;
      }
      if (event.key !== "Escape") return;
      if (!document.getElementById("patient-document-preview")?.hidden) closePreview();
      else if (!document.getElementById("patient-documents-modal")?.hidden) close();
    });
    window.addEventListener("resize", applyPreviewViewport);
  }

  async function open(event) {
    event?.preventDefault?.();
    event?.stopPropagation?.();
    injectShell();
    const raw = typeof options.getContext === "function" ? options.getContext() : {};
    const mrn = String(raw?.mrn || "").trim();
    if (!MRN_RE.test(mrn)) {
      alert("請先填寫有效的病歷號，再開啟 Documents。");
      return;
    }
    context = {
      mrn,
      sourceSampleId: String(raw?.sourceSampleId || "").trim(),
    };
    renderList([]);
    document.getElementById("pdoc-patient").textContent = `MRN：${mrn}`;
    document.getElementById("patient-documents-modal").hidden = false;
    setStatus("pdoc-upload-status", "");
    setStatus("pdoc-list-status", "");
    if (await ensureAuthenticated()) await loadList();
  }

  function close() {
    closePreview();
    clearPending();
    context = null;
    previewDocuments = [];
    const modal = document.getElementById("patient-documents-modal");
    if (modal) modal.hidden = true;
  }

  function init(initOptions) {
    options = { ...options, ...(initOptions || {}) };
    injectShell();
    const button = typeof options.button === "string"
      ? document.querySelector(options.button)
      : options.button;
    button?.addEventListener("click", open);
  }

  window.PatientDocuments = { init, open, close };
})();
