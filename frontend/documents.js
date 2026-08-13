(function () {
  "use strict";

  const API = "/api/documents";
  const ACCEPT = ".pdf,.jpg,.jpeg,.png,.tif,.tiff,application/pdf,image/jpeg,image/png,image/tiff";
  const MRN_RE = /^[A-Za-z0-9_-]{1,32}$/;
  let options = {};
  let context = null;
  let pending = [];
  let previewUrl = "";
  let previewDocument = null;
  let previewPage = 0;

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
          <div id="pdoc-preview-body" class="pdoc-preview-body"><span class="pdoc-muted">載入中…</span></div>
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
    setStatus("pdoc-list-status", "載入中…");
    try {
      const rows = await jsonRequest(`${API}?mrn=${encodeURIComponent(context.mrn)}`);
      renderList(Array.isArray(rows) ? rows : []);
      setStatus("pdoc-list-status", "");
    } catch (error) {
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

  async function loadPreviewPage(page) {
    if (!previewDocument) return;
    const body = document.getElementById("pdoc-preview-body");
    body.innerHTML = `<span class="pdoc-muted">載入中…</span>`;
    releasePreviewUrl();
    try {
      const response = await fetch(
        `${API}/${encodeURIComponent(previewDocument.id)}/preview?page=${page}`,
        { credentials: "same-origin" },
      );
      if (response.status === 401) {
        closePreview();
        showLogin();
        return;
      }
      if (!response.ok) {
        const error = await response.json().catch(() => ({}));
        throw new Error(error.detail || `預覽失敗 (${response.status})`);
      }
      previewUrl = URL.createObjectURL(await response.blob());
      previewPage = page;
      body.innerHTML = `<img src="${esc(previewUrl)}" alt="${esc(previewDocument.display_name)}">`;
      renderPreviewControls();
    } catch (error) {
      body.innerHTML = `<div class="pdoc-error">${esc(error.message || error)}</div>`;
    }
  }

  function renderPreviewControls() {
    const pages = Number(previewDocument?.image_pages || 1);
    const controls = document.getElementById("pdoc-preview-controls");
    const label = document.getElementById("pdoc-preview-page-label");
    controls.hidden = pages <= 1;
    label.textContent = pages > 1 ? `第 ${previewPage + 1} / ${pages} 頁` : "";
    controls.querySelector("[data-pdoc-preview-prev]").disabled = previewPage <= 0;
    controls.querySelector("[data-pdoc-preview-next]").disabled = previewPage >= pages - 1;
  }

  function openPreview(rowElement) {
    previewDocument = JSON.parse(rowElement.dataset.document || "{}");
    previewPage = 0;
    document.getElementById("pdoc-preview-title").textContent = previewDocument.display_name || "圖片預覽";
    document.getElementById("patient-document-preview").hidden = false;
    loadPreviewPage(0);
  }

  function closePreview() {
    releasePreviewUrl();
    previewDocument = null;
    document.getElementById("patient-document-preview").hidden = true;
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
    document.getElementById("patient-documents-modal")?.addEventListener("click", event => {
      if (event.target.id === "patient-documents-modal") close();
    });
    document.getElementById("patient-document-preview")?.addEventListener("click", event => {
      if (event.target.id === "patient-document-preview") closePreview();
    });
    document.addEventListener("keydown", event => {
      if (event.key !== "Escape") return;
      if (!document.getElementById("patient-document-preview")?.hidden) closePreview();
      else if (!document.getElementById("patient-documents-modal")?.hidden) close();
    });
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
