const assert = require('node:assert/strict');
const { mkdtempSync, writeFileSync, readFileSync, rmSync } = require('node:fs');
const { tmpdir } = require('node:os');
const { join } = require('node:path');
const { spawnSync } = require('node:child_process');
const { test } = require('node:test');

function fixture(t) {
  const directory = mkdtempSync(join(tmpdir(), 'snapshot-action-'));
  t.after(() => rmSync(directory, { recursive: true, force: true }));
  const state = join(directory, 'state');
  const calls = join(directory, 'calls');
  writeFileSync(state, '');
  writeFileSync(calls, '');
  writeFileSync(join(directory, 'uv'), `#!${process.execPath}
require('node:fs').appendFileSync(process.env.TEST_CALLS,
  JSON.stringify({args: process.argv.slice(2), env: process.env}) + '\\n');
process.exit(Number(process.env.TEST_EXIT || 0));
`, { mode: 0o700 });
  const env = {
    PATH: directory, GITHUB_STATE: state, TEST_CALLS: calls,
    GITHUB_EVENT_NAME: 'push', GITHUB_REF: 'refs/heads/main',
    CCACHE_DIR: join(directory, 'cache'),
    'INPUT_KEY-PREFIX': 'pilot/',
    'INPUT_ACCESS-KEY-ID': 'test-key', 'INPUT_SECRET-ACCESS-KEY': 'test-secret',
  };
  function run(overrides = {}) {
    const result = spawnSync(process.execPath, [join(__dirname, 'action.cjs')], {
      env: { ...env, ...overrides }, encoding: 'utf8',
    });
    assert.equal(result.status, 0, result.stderr);
    assert.ok(!result.stdout.includes('test-secret'));
    return result.stdout;
  }
  return {
    env, run,
    state: () => readFileSync(state, 'utf8').trim().slice('snapshot='.length),
    calls: () => readFileSync(calls, 'utf8').trim().split('\n').filter(Boolean).map(JSON.parse),
  };
}

test('restore is anonymous; post reuses saved cache and credentials', (t) => {
  const f = fixture(t);
  f.run();
  const restore = f.calls()[0];
  assert.equal(restore.args.at(-1), 'restore');
  assert.ok(restore.args.includes('--locked'));
  assert.equal(restore.env.AWS_ACCESS_KEY_ID, '');
  assert.equal(restore.env.AWS_SECRET_ACCESS_KEY, '');
  assert.equal(restore.env['INPUT_SECRET-ACCESS-KEY'], undefined);
  f.run({ STATE_snapshot: f.state(), CCACHE_DIR: '/changed', 'INPUT_KEY-PREFIX': 'changed/' });
  const publish = f.calls()[1];
  assert.equal(publish.args.at(-1), 'publish');
  assert.equal(publish.env.CCACHE_DIR, f.env.CCACHE_DIR);
  assert.equal(publish.env.CACHE_PREFIX, 'pilot/');
  assert.equal(publish.env.AWS_SECRET_ACCESS_KEY, 'test-secret');
});

for (const overrides of [
  { GITHUB_EVENT_NAME: 'pull_request', GITHUB_REF: 'refs/pull/1/merge' },
  { GITHUB_EVENT_NAME: 'push', GITHUB_REF: 'refs/heads/feature' },
  { CCACHE_SNAPSHOT_SKIP_SAVE: 'true' },
]) {
  test(`post skips publication: ${JSON.stringify(overrides)}`, (t) => {
    const f = fixture(t);
    f.run(overrides);
    f.run({ ...overrides, STATE_snapshot: f.state() });
    assert.equal(f.calls().length, 1);
  });
}

test('restore failure remains nonfatal and still allows a later upload', (t) => {
  const f = fixture(t);
  assert.match(f.run({ TEST_EXIT: '1' }), /::warning::/);
  f.run({ STATE_snapshot: f.state() });
  assert.equal(f.calls()[1].args.at(-1), 'publish');
});

test('post failure remains nonfatal', (t) => {
  const f = fixture(t);
  f.run();
  assert.match(f.run({ STATE_snapshot: f.state(), TEST_EXIT: '1' }), /::warning::/);
});
