// NGS OMIM Capture — background service worker.
//
// All platform traffic goes through here: the content script on omim.org
// cannot reach the hospital intranet (mixed content / CORS), while an
// extension worker with host permission can. Auth is a personal capture
// token created on the platform (/omim-capture/), never the login cookie.

const DEFAULT_PLATFORM = "http://192.168.84.91:8765";
const OUTBOX_KEY = "outbox";
const OUTBOX_MAX = 50;

async function settings() {
  const s = await chrome.storage.local.get(["platformUrl", "token"]);
  return {
    platformUrl: String(s.platformUrl || DEFAULT_PLATFORM).replace(/\/+$/, ""),
    token: String(s.token || "").trim(),
  };
}

async function api(path, { method = "GET", body } = {}) {
  const { platformUrl, token } = await settings();
  if (!token) {
    return { ok: false, status: 0, error: "尚未設定 capture token：請到擴充功能設定頁貼上平台產生的 token。" };
  }
  let response;
  try {
    response = await fetch(platformUrl + path, {
      method,
      headers: {
        Authorization: `Bearer ${token}`,
        ...(body ? { "Content-Type": "application/json" } : {}),
      },
      body: body ? JSON.stringify(body) : undefined,
      credentials: "omit",
    });
  } catch (err) {
    return { ok: false, status: 0, network: true, error: `連不到平台 ${platformUrl}（${err.message || err}）` };
  }
  let data = null;
  try { data = await response.json(); } catch (_) { /* non-JSON error page */ }
  if (!response.ok) {
    const detail = data && data.detail ? data.detail : response.statusText;
    return { ok: false, status: response.status, error: `平台回應 ${response.status}：${detail}` };
  }
  return { ok: true, status: response.status, data };
}

async function readOutbox() {
  const s = await chrome.storage.local.get([OUTBOX_KEY]);
  return Array.isArray(s[OUTBOX_KEY]) ? s[OUTBOX_KEY] : [];
}

async function writeOutbox(items) {
  await chrome.storage.local.set({ [OUTBOX_KEY]: items.slice(-OUTBOX_MAX) });
}

// Retry captures that could not be delivered earlier. The server dedupes
// identical content, so a retry after an ambiguous failure is harmless.
async function flushOutbox() {
  const items = await readOutbox();
  if (!items.length) return { sent: 0, left: 0 };
  const left = [];
  let sent = 0;
  for (const item of items) {
    const result = await api("/api/omim/captures", { method: "POST", body: item });
    if (result.ok) sent += 1;
    else if (result.network) left.push(item);
    // Non-network errors (4xx) are dropped: retrying would fail the same way.
  }
  await writeOutbox(left);
  return { sent, left: left.length };
}

async function capture(payload) {
  const result = await api("/api/omim/captures", { method: "POST", body: payload });
  if (!result.ok && result.network) {
    const items = await readOutbox();
    items.push(payload);
    await writeOutbox(items);
    return { ...result, queued: true, outbox: Math.min(items.length, OUTBOX_MAX) };
  }
  if (result.ok) {
    const flushed = await flushOutbox();
    return { ...result, outboxSent: flushed.sent, outbox: flushed.left };
  }
  return result;
}

chrome.runtime.onMessage.addListener((message, _sender, sendResponse) => {
  (async () => {
    switch (message && message.type) {
      case "status":
        return api(`/api/omim/status?mim=${encodeURIComponent(message.mim)}`
          + `&edit_date=${encodeURIComponent(message.editDate || "")}`);
      case "capture":
        return capture(message.payload);
      case "adopt": {
        const [mim, name] = String(message.captureId || "").split("/");
        return api(`/api/omim/captures/${encodeURIComponent(mim)}/${encodeURIComponent(name)}/adopt`,
          { method: "POST" });
      }
      case "ping":
        return api("/api/omim/ping");
      case "outbox": {
        const flushed = await flushOutbox();
        return { ok: true, data: flushed };
      }
      case "open-options":
        await chrome.runtime.openOptionsPage();
        return { ok: true };
      default:
        return { ok: false, error: "unknown message" };
    }
  })().then(sendResponse, err => sendResponse({ ok: false, error: String(err) }));
  return true; // async response
});

chrome.action.onClicked.addListener(() => chrome.runtime.openOptionsPage());
