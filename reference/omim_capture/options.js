const DEFAULT_PLATFORM = "http://192.168.84.91:8765";
const $ = id => document.getElementById(id);

function show(ok, message) {
  const box = $("result");
  box.className = ok ? "ok" : "err";
  box.textContent = message;
}

function originPattern(url) {
  const u = new URL(url);
  return `${u.protocol}//${u.hostname}/*`;
}

chrome.storage.local.get(["platformUrl", "token", "queueScope", "queueSkip"], s => {
  $("platformUrl").value = s.platformUrl || DEFAULT_PLATFORM;
  $("token").value = s.token || "";
  $("queueScope").value = s.queueScope || "all";
  $("skipCount").textContent = Array.isArray(s.queueSkip) ? s.queueSkip.length : 0;
});

$("queueScope").addEventListener("change", e => chrome.storage.local.set({ queueScope: e.target.value }));
$("clearSkip").addEventListener("click", async e => {
  e.preventDefault();
  await chrome.storage.local.set({ queueSkip: [] });
  $("skipCount").textContent = "0";
});

$("save").addEventListener("click", async () => {
  let platformUrl = $("platformUrl").value.trim().replace(/\/+$/, "") || DEFAULT_PLATFORM;
  const token = $("token").value.trim();
  try {
    new URL(platformUrl);
  } catch (_) {
    show(false, "平台網址格式不正確。");
    return;
  }
  // Ask for host permission when the platform is not the built-in default.
  const pattern = originPattern(platformUrl);
  const granted = await chrome.permissions.contains({ origins: [pattern] })
    || await chrome.permissions.request({ origins: [pattern] });
  if (!granted) {
    show(false, `沒有取得連線到 ${pattern} 的權限。`);
    return;
  }
  await chrome.storage.local.set({ platformUrl, token });
  chrome.runtime.sendMessage({ type: "ping" }, result => {
    if (chrome.runtime.lastError) {
      show(false, chrome.runtime.lastError.message);
    } else if (result && result.ok) {
      show(true, `已連線，登入帳號：${result.data.username}`);
      chrome.runtime.sendMessage({ type: "outbox" });
    } else {
      show(false, (result && result.error) || "連線失敗");
    }
  });
});

$("open").addEventListener("click", () => {
  const base = ($("platformUrl").value.trim() || DEFAULT_PLATFORM).replace(/\/+$/, "");
  chrome.tabs.create({ url: `${base}/omim-capture/` });
});
