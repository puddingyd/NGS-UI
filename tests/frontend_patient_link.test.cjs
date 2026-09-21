const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const test = require('node:test');

const source = fs.readFileSync(path.join(__dirname, '../frontend/phenotype/app.js'), 'utf8');
const main = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');

function form({ mrn = '00123456', code = '8BB126WE0092', clinical = '', hpo = [], conflict = false } = {}) {
  const calls = [];
  const statuses = [];
  const fields = {
    'patient-mrn': { value: mrn },
    'patient-code': { value: code },
    'clinical-presentation-text': { value: clinical },
  };
  const button = { disabled: false };
  const context = vm.createContext({
    document: { getElementById: id => fields[id], querySelectorAll: () => [button] },
    selectedFixedPanels: new Map(),
    loadedClinicalPresentationSidecar: Boolean(clinical),
    clinicalPresentationLastSaved: clinical,
    clinicalPresentationLastSavedPath: '/old/clinical.txt',
    loadedPhenotypeSidecar: false,
    clinicalAutosaveTimer: null,
    clearTimeout: () => {},
    _collectCustomPanels: () => [],
    _collectHpoAndPanelLines: () => hpo,
    updatePreview: () => {},
    showStatus: (message, type) => statuses.push({ message, type }),
    fetch: async (url, options) => {
      calls.push({ url, payload: JSON.parse(options.body) });
      const link = url.endsWith('/patient-link');
      return {
        ok: !(link && conflict),
        json: async () => link
          ? conflict ? { detail: '此檢體編號已連結其他病歷號' } : { lis_id: '26WE0092', mrn }
          : { path: `/saved/${mrn}.txt` },
      };
    },
  });
  vm.runInContext(source.slice(source.indexOf('function _clinicalPresentationFields()'),
    source.indexOf('function scheduleClinicalPresentationAutosave()')), context);
  vm.runInContext(source.slice(source.indexOf('async function generateFile()'),
    source.indexOf('function showInlineClinicalStatus(')), context);
  return { context, calls, statuses, button };
}

for (const [name, options, expected] of [
  ['HPO', { hpo: ['HP:0001250\tSeizure\t1'] }, ['/api/phenotype-tool/save', '/api/phenotype-tool/patient-link']],
  ['unchanged clinical text', { clinical: 'Clinical test text' }, ['/api/phenotype-tool/clinical-presentation/save', '/api/phenotype-tool/patient-link']],
  ['identifiers alone', {}, ['/api/phenotype-tool/patient-link']],
]) {
  test(`Save links both identifiers with ${name}`, async () => {
    const { context, calls, statuses, button } = form(options);
    await context.generateFile();
    assert.deepEqual(calls.map(call => call.url), expected);
    assert.deepEqual(calls.at(-1).payload, { mrn: '00123456', code: '8BB126WE0092' });
    assert.equal(statuses.at(-1).type, 'success');
    assert.match(statuses.at(-1).message, /已連結檢體 26WE0092/);
    assert.equal(button.disabled, false);
  });
}

test('Saving with only MRN keeps the existing phenotype flow without a link', async () => {
  const { context, calls } = form({ code: '', hpo: ['HP:0001250\tSeizure\t1'] });
  await context.generateFile();
  assert.deepEqual(calls.map(call => call.url), ['/api/phenotype-tool/save']);
});

test('Clinical autosave never creates an identity link', async () => {
  const { context, calls } = form({ clinical: 'Clinical test text' });
  context.clinicalPresentationLastSaved = '';
  await context.saveClinicalPresentationSidecar();
  assert.deepEqual(calls.map(call => call.url), ['/api/phenotype-tool/clinical-presentation/save']);
});

test('Conflicting linkage is shown as a save error', async () => {
  const { context, statuses, button } = form({ conflict: true });
  await context.generateFile();
  assert.equal(statuses.at(-1).type, 'error');
  assert.match(statuses.at(-1).message, /連結未儲存/);
  assert.equal(button.disabled, false);
});
