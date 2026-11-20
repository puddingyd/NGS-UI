// Synthetic WES/WGS export dialogs; no patient data or running backend.
const { chromium } = require('playwright');
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');

const root = path.resolve(__dirname, '..');
const source = fs.readFileSync(path.join(root, 'frontend/app.js'), 'utf8');
const start = source.indexOf('function _pickGeneListMode(');
const end = source.indexOf('async function handleLogout()', start);
const harness = `
  let testType = 'WES';
  function currentSampleTestType() { return testType; }
  function escapeHtml(value) {
    return String(value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  }
  ${source.slice(start, end)}
  function openDialog(type, opts = {}) {
    testType = type;
    window.pickPromise = _pickGeneListMode(opts);
  }
`;

(async () => {
  const browser = await chromium.launch({
    headless: true,
    executablePath: process.env.CHROME_EXECUTABLE_PATH || '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
  });
  try {
    const page = await browser.newPage({ viewport: { width: 1280, height: 900 } });
    await page.setContent('<!doctype html><meta name="viewport" content="width=device-width,initial-scale=1">');
    await page.addStyleTag({ content: fs.readFileSync(path.join(root, 'frontend/style.css'), 'utf8') });
    await page.addScriptTag({ content: harness });
    await page.evaluate(() => openDialog('WES'));
    const note = page.locator('.wes-gene-list-exclusions');
    assert.equal(await note.innerText(), 'WES 基因清單排除：mtDNA、SMN1/SMN2、STRC、CYP21A2、NCF1、IKBKG、OPN1LW/OPN1MW、MUC1');
    const layout = await note.evaluate(el => {
      const question = el.parentElement.querySelector('p');
      const buttons = el.nextElementSibling;
      const range = document.createRange();
      range.selectNodeContents(el);
      return {
        font: getComputedStyle(el).fontSize,
        questionFont: getComputedStyle(question).fontSize,
        color: getComputedStyle(el).color,
        questionColor: getComputedStyle(question).color,
        leftDelta: el.getBoundingClientRect().left - question.getBoundingClientRect().left,
        aboveButtons: el.getBoundingClientRect().bottom <= buttons.getBoundingClientRect().top,
        lines: range.getClientRects().length,
        modalWidth: el.parentElement.getBoundingClientRect().width,
      };
    });
    assert.equal(layout.font, layout.questionFont);
    assert.equal(layout.color, layout.questionColor);
    assert.equal(layout.leftDelta, 0);
    assert.equal(layout.aboveButtons, true);
    assert(layout.lines > 1);
    assert.equal(layout.modalWidth, 520);
    await page.screenshot({ path: '/private/tmp/wes-report-modal-desktop.png' });
    await page.locator('input[value="merged"]').check();
    await page.locator('[data-act="ok"]').click();
    assert.equal(await page.evaluate(() => window.pickPromise), 'merged');

    for (const type of ['WGS', 'TITAN-WGS']) {
      await page.evaluate(type => openDialog(type), type);
      assert.equal(await note.count(), 0);
      await page.locator('[data-act="cancel"]').click();
      assert.equal(await page.evaluate(() => window.pickPromise), null);
    }

    await page.setViewportSize({ width: 390, height: 844 });
    await page.evaluate(() => openDialog('WES', { title: '輸出 PDF', okLabel: '輸出 PDF' }));
    assert.equal(await note.isVisible(), true);
    assert.equal(await page.locator('[data-act="ok"]').isVisible(), true);
    assert(await page.locator('.modal-card').evaluate(el => el.scrollWidth <= el.clientWidth));
    await page.screenshot({ path: '/private/tmp/wes-report-modal-mobile.png' });
    await page.locator('[data-act="ok"]').click();
    assert.equal(await page.evaluate(() => window.pickPromise), 'grouped');
    console.log('WES/WGS dialogs, note layout, mobile layout, and export/cancel behavior passed.');
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
