const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const test = require('node:test');

const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');
const syncSource = source.slice(source.indexOf('let _newCaseEmrRequestId ='),
  source.indexOf('function renderNewCasePhenoEditor()'));

function setup({ saved = null, clinical = null, emr = { phenotype: {} } } = {}) {
  const calls = [];
  const classList = () => {
    const classes = new Set();
    return { add: key => classes.add(key), remove: key => classes.delete(key), contains: key => classes.has(key) };
  };
  const fields = Object.fromEntries([
    'new-case-mrn', 'new-case-lis-id', 'new-case-error', 'new-case-modal',
    'new-case-clinical-preview-section', 'new-case-clinical-preview', 'btn-new-case-emr',
  ].map(id => [id, { value: '', classList: classList(), addEventListener() {} }]));
  fields['new-case-mrn'].value = '00123456';
  fields['new-case-lis-id'].value = 'NEW-SAMPLE';
  const sex = { value: '' };
  const edit = {
    hpo: [{ phenotype: 'HP:OLD' }], panels: [{ name: 'OldPanel' }],
    source: 'Web phenotype input tool', edited: false,
  };
  const context = vm.createContext({
    document: { getElementById: id => fields[id], querySelector: () => sex },
    newCaseEdit: edit,
    renderNewCasePhenoEditor() {}, renderNewCaseEmrRef() {},
    apiFetch: async url => {
      calls.push(url);
      const response = url.startsWith('/emr/') ? emr : url.includes('/clinical-presentation/') ? clinical : saved;
      if (response instanceof Error) throw response;
      return response;
    },
  });
  vm.runInContext(syncSource, context);
  return { context, edit, fields, sex, calls };
}

test('sync reads both saved files by MRN only; saved phenotype wins over EMR and stale cache', async () => {
  const state = setup({
    saved: { hpo: [{ phenotype: 'HP:0001250', label: 'Seizure', weight: 3 }], panels: [{ name: 'Neuro', weight: 2 }] },
    clinical: { content: '<clinical text>\n' },
    emr: { consultation: { sex: 'F' }, phenotype: { found: true, hpo: [{ phenotype: 'HP:EMR' }] } },
  });
  await state.context.syncNewCaseEmr();
  assert.deepEqual(state.calls, ['/phenotype-tool/load?mrn=00123456',
    '/phenotype-tool/clinical-presentation/load?mrn=00123456', '/emr/00123456']);
  assert.equal(state.edit.hpo[0].phenotype, 'HP:0001250');
  assert.equal(state.edit.hpo[0].weight, 3);
  assert.equal(state.edit.panels[0].name, 'Neuro');
  assert.equal(state.edit.panels[0].weight, 2);
  assert.equal(state.edit.edited, true);
  assert.equal(state.edit.emrPhenotype.hpo[0].phenotype, 'HP:EMR');
  assert.equal(state.sex.value, 'F');
  assert.equal(state.fields['new-case-clinical-preview'].value, '<clinical text>\n');
  assert.equal(state.fields['new-case-clinical-preview-section'].hidden, false);
});

test('an explicitly empty saved snapshot remains authoritative over EMR', async () => {
  const state = setup({ saved: { hpo: [], panels: [] }, emr: { phenotype: { hpo: [{ phenotype: 'HP:EMR' }] } } });
  await state.context.syncNewCaseEmr();
  assert.equal(state.edit.hpo.length, 0);
  assert.equal(state.edit.panels.length, 0);
  assert.equal(state.edit.edited, true);
});

test('saved data and clinical text still load when EMR fails', async () => {
  const state = setup({ saved: { panels: [{ name: 'Neuro' }] }, clinical: { content: 'Saved text' }, emr: new Error('offline') });
  await state.context.syncNewCaseEmr();
  assert.equal(state.edit.panels[0].name, 'Neuro');
  assert.equal(state.fields['new-case-clinical-preview'].value, 'Saved text');
  assert.match(state.fields['new-case-error'].textContent, /EMR 同步失敗/);
  assert.equal(state.fields['btn-new-case-emr'].disabled, false);
});

test('no saved phenotype falls back to EMR and removes stale panels', async () => {
  const state = setup({ emr: { phenotype: { hpo: [{ phenotype: 'HP:EMR' }] } } });
  await state.context.syncNewCaseEmr();
  assert.equal(state.edit.hpo[0].phenotype, 'HP:EMR');
  assert.equal(state.edit.panels.length, 0);
  assert.equal(state.fields['new-case-clinical-preview-section'].hidden, true);
});

test('a file-read failure is not treated as absence and does not replace chips with EMR', async () => {
  const state = setup({ saved: new Error('500'), clinical: new Error('500'), emr: { phenotype: { hpo: [{ phenotype: 'HP:EMR' }] } } });
  await state.context.syncNewCaseEmr();
  assert.equal(state.edit.hpo[0].phenotype, 'HP:OLD');
  assert.match(state.fields['new-case-error'].textContent, /HPO／panel 讀取失敗/);
  assert.match(state.fields['new-case-error'].textContent, /Clinical presentation 讀取失敗/);
});

for (const change of ['mrn', 'sample', 'reopen', 'close']) {
  test(`a late response after ${change} cannot update the form`, async () => {
    let resolve;
    const saved = new Promise(done => { resolve = done; });
    const state = setup({ saved, clinical: { content: 'Old patient text' } });
    const pending = state.context.syncNewCaseEmr();
    if (change === 'mrn') state.fields['new-case-mrn'].value = 'DIFFERENT';
    if (change === 'sample') state.fields['new-case-lis-id'].value = 'OTHER';
    if (change === 'reopen') state.context.resetNewCaseSync();
    if (change === 'close') state.fields['new-case-modal'].classList.add('hidden');
    resolve({ hpo: [{ phenotype: 'HP:LATE' }] });
    await pending;
    assert.equal(state.edit.hpo[0].phenotype, 'HP:OLD');
    assert.equal(state.fields['new-case-clinical-preview'].value, '');
  });
}

test('empty MRN does not send requests', async () => {
  const state = setup();
  state.fields['new-case-mrn'].value = ' ';
  await state.context.syncNewCaseEmr();
  assert.equal(state.calls.length, 0);
  assert.match(state.fields['new-case-error'].textContent, /請先填 MRN/);
});
