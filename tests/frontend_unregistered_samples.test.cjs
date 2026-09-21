const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const test = require('node:test');
const source = fs.readFileSync(path.join(__dirname, '../frontend/app.js'), 'utf8');

function deferred() {
  let resolve, reject;
  const promise = new Promise((yes, no) => { resolve = yes; reject = no; });
  return { promise, resolve, reject };
}

function setup(api) {
  const calls = [];
  const timers = [];
  const classList = () => {
    const values = new Set();
    return { add: k => values.add(k), remove: k => values.delete(k), contains: k => values.has(k) };
  };
  const fields = Object.fromEntries([
    'new-case-mrn', 'new-case-name', 'new-case-lis-id', 'new-case-lis-id-search', 'new-case-lis-id-dropdown',
    'new-case-modal', 'new-case-error', 'new-case-list-status', 'new-case-dept-hint', 'new-case-pheno-source', 'btn-new-case-emr',
  ].map(id => [id, { value: '', textContent: '', classList: classList() }]));
  fields['new-case-lis-id-dropdown'].classList.add('hidden');
  const testType = { value: 'WES' };
  const context = vm.createContext({
    document: { getElementById: id => fields[id], querySelector: () => testType },
    clearTimeout() {}, setTimeout(fn) { timers.push(fn); return timers.length; },
    apiFetch: async url => { calls.push(url); return api(url); },
    resetNewCaseSync() {}, _updateNewCaseEmrLink() {}, renderNewCasePhenoEditor() {},
    _renderNewCaseLisDropdown() {},
    inferNewCaseTestType: entry => entry.source_vcf_size > 100 * 1024 * 1024 ? 'WGS' : '',
  });
  vm.runInContext(source.slice(source.indexOf('let _unregisteredById ='), source.indexOf('const NEW_CASE_WGS_VCF_SIZE_BYTES')), context);
  vm.runInContext('globalThis.edit = newCaseEdit; globalThis.cache = _unregisteredCache;', context);
  vm.runInContext(source.slice(source.indexOf('function _setUnregisteredList('), source.indexOf('document.getElementById("btn-new-case")?.addEventListener')), context);
  vm.runInContext(source.slice(source.indexOf('function _fillNewCaseAutoField('), source.indexOf('// EMR sync button on the modal:')), context);
  const seed = rows => {
    context.cache.list = rows;
    context._setUnregisteredList(rows);
  };
  const select = id => { fields['new-case-lis-id'].value = id; return context._applyNewCaseLisSelection(id); };
  return { context, fields, testType, calls, timers, seed, select, edit: context.edit };
}

const row = (id, mrn = '') => ({ lis_id: id, roster: { mrn, name: id } });
const detail = (id, mrn = '') => ({ ...row(id, mrn), phenotype: { mrn, hpo: [{ phenotype: `HP:${id}`, weight: 2 }], panels: [{ name: 'Neuro' }] }, source_vcf_size: 200 * 1024 * 1024 });

test('cached list is returned immediately while a slow refresh runs', async () => {
  const wait = deferred();
  const state = setup(() => wait.promise);
  state.seed([row('S1')]);
  const result = await state.context._loadUnregisteredSamples();
  assert.equal(result[0].lis_id, 'S1');
  assert.deepEqual(state.calls, ['/samples/unregistered']);
  wait.resolve({ items: [row('S2')], updated_at: 1 });
  await state.context._fetchUnregisteredSamples();
  assert.equal(state.context.cache.list[0].lis_id, 'S2');
});

test('roster changes refresh MRN without checking a separate revision or rescanning endpoint', async () => {
  const state = setup(() => ({ items: [row('S1', 'NEW-MRN')], updated_at: 1 }));
  state.seed([row('S1')]);
  await state.context._loadUnregisteredSamples();
  await state.context._fetchUnregisteredSamples();
  assert.equal(state.context.cache.list[0].roster.mrn, 'NEW-MRN');
  assert.ok(state.calls.every(url => url === '/samples/unregistered'));
});

test('a cold response schedules one poll and then displays the completed index', async () => {
  let ready = false;
  const state = setup(() => ready ? { items: [row('S1')], updated_at: 1 } : { items: [], refreshing: true });
  await state.context._loadUnregisteredSamples();
  assert.equal(state.timers.length, 1);
  assert.match(state.fields['new-case-list-status'].textContent, /背景更新/);
  ready = true;
  state.timers[0]();
  await state.context._fetchUnregisteredSamples();
  assert.equal(state.context.cache.list[0].lis_id, 'S1');
});

test('background failure preserves the list and manual refresh really requests a rescan', async () => {
  let fail = true;
  const state = setup(() => { if (fail) throw new Error('offline'); return { items: [row('S2')], updated_at: 1 }; });
  state.seed([row('S1')]);
  await assert.rejects(state.context._fetchUnregisteredSamples());
  assert.equal(state.context.cache.list[0].lis_id, 'S1');
  assert.match(state.fields['new-case-list-status'].textContent, /上次清單/);
  fail = false;
  await state.context._loadUnregisteredSamples({ force: true });
  assert.equal(state.calls.at(-1), '/samples/unregistered?refresh=true');
  assert.equal(state.context.cache.list[0].lis_id, 'S2');
});

test('a late background reply cannot overwrite a newer forced refresh', async () => {
  const old = deferred();
  const state = setup(url => url.includes('refresh=true') ? { items: [row('NEW')], updated_at: 2 } : old.promise);
  const pending = state.context._fetchUnregisteredSamples();
  await state.context._fetchUnregisteredSamples({ force: true });
  old.resolve({ items: [row('OLD')], updated_at: 1 });
  await pending;
  assert.equal(state.context.cache.list[0].lis_id, 'NEW');
});

test('a list request begun before registration cannot re-add the registered row', async () => {
  const wait = deferred();
  const state = setup(() => wait.promise);
  state.seed([row('S1')]);
  const pending = state.context._fetchUnregisteredSamples();
  state.context._removeUnregisteredFromCache('S1');
  wait.resolve({ items: [row('S1')], updated_at: 1 });
  await pending;
  assert.equal(state.context.cache.list.length, 0);
});

test('selection fills roster immediately, then retrieves only that sample and restores phenotype/type', async () => {
  const wait = deferred();
  const state = setup(() => wait.promise);
  state.seed([row('S1', 'MRN1'), row('S2', 'MRN2')]);
  const pending = state.select('S1');
  assert.equal(state.fields['new-case-mrn'].value, 'MRN1');
  assert.equal(state.edit.detailLoading, true);
  assert.deepEqual(state.calls, ['/samples/unregistered/S1']);
  wait.resolve(detail('S1', 'MRN1'));
  await pending;
  assert.equal(state.edit.hpo[0].phenotype, 'HP:S1');
  assert.equal(state.edit.panels[0].name, 'Neuro');
  assert.equal(state.testType.value, 'WGS');
  assert.equal(state.edit.detailLoading, false);
});

test('quickly switching samples replaces auto-filled MRN and ignores the old detail reply', async () => {
  const old = deferred();
  const state = setup(url => url.endsWith('S1') ? old.promise : detail('S2', 'MRN2'));
  state.seed([row('S1', 'MRN1'), row('S2', 'MRN2')]);
  const first = state.select('S1');
  await state.select('S2');
  old.resolve(detail('S1', 'MRN1'));
  await first;
  assert.equal(state.fields['new-case-mrn'].value, 'MRN2');
  assert.equal(state.edit.hpo[0].phenotype, 'HP:S2');
});

test('explicit MRN is sent to the detail lookup, and manual chips/type changes are preserved', async () => {
  const wait = deferred();
  const state = setup(() => wait.promise);
  state.seed([row('S1', 'MRN1')]);
  state.fields['new-case-mrn'].value = 'MANUAL';
  const pending = state.select('S1');
  state.edit.edited = true;
  state.edit.hpo = [{ phenotype: 'HP:MANUAL' }];
  state.testType.value = 'TITAN-WGS';
  wait.resolve(detail('S1', 'MRN1'));
  await pending;
  assert.equal(state.calls[0], '/samples/unregistered/S1?mrn=MANUAL');
  assert.equal(state.fields['new-case-mrn'].value, 'MANUAL');
  assert.equal(state.edit.hpo[0].phenotype, 'HP:MANUAL');
  assert.equal(state.testType.value, 'TITAN-WGS');
});

test('an unavailable or already registered selection becomes an error and leaves no stale chips', async () => {
  const state = setup(() => null);
  state.seed([row('S1')]);
  await state.select('S1');
  assert.match(state.edit.detailError, /已登錄/);
  assert.equal(state.edit.hpo.length, 0);
  assert.equal(state.context.cache.list.length, 0);
});

test('clearing identity fields while details load is not undone by the response', async () => {
  const wait = deferred();
  const state = setup(() => wait.promise);
  state.seed([row('S1', 'MRN1')]);
  const pending = state.select('S1');
  state.fields['new-case-mrn'].value = '';
  state.fields['new-case-name'].value = '';
  wait.resolve(detail('S1', 'MRN1'));
  await pending;
  assert.equal(state.fields['new-case-mrn'].value, '');
  assert.equal(state.fields['new-case-name'].value, '');
  assert.equal(state.edit.hpo.length, 0);
});
