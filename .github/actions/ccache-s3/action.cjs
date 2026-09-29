const { appendFileSync } = require('node:fs');
const { resolve } = require('node:path');
const { spawnSync } = require('node:child_process');

// GitHub restores GITHUB_STATE entries as STATE_* only for this action's post hook.
const post = process.env.STATE_snapshot !== undefined;
const operation = post ? 'publish' : 'restore';
try {
  const config = post ? JSON.parse(process.env.STATE_snapshot) : {
    CCACHE_DIR: process.env.CCACHE_DIR,
    CACHE_PREFIX: process.env['INPUT_KEY-PREFIX'],
    CACHE_ENDPOINT: 'https://s3.cern.ch',
    CACHE_REGION: 'cern',
    CACHE_BUCKET: 'cache',
    AWS_ACCESS_KEY_ID: process.env['INPUT_ACCESS-KEY-ID'] || '',
    AWS_SECRET_ACCESS_KEY: process.env['INPUT_SECRET-ACCESS-KEY'] || '',
    AWS_SESSION_TOKEN: '',
    AWS_EC2_METADATA_DISABLED: 'true',
  };
  if (!post) {
    appendFileSync(process.env.GITHUB_STATE, `snapshot=${JSON.stringify(config)}\n`);
  }
  if (post && (process.env.GITHUB_EVENT_NAME !== 'push'
      || process.env.GITHUB_REF !== 'refs/heads/main'
      || process.env.CCACHE_SNAPSHOT_SKIP_SAVE === 'true')) {
    console.log('ccache publication skipped: read-only run or cleanup failed');
  } else {
    if (!config.CCACHE_DIR || !config.CACHE_PREFIX) {
      throw new Error('Missing cache configuration');
    }
    // Anonymous restore never passes write credentials to the transfer process.
    const env = { ...process.env, ...config };
    delete env['INPUT_ACCESS-KEY-ID'];
    delete env['INPUT_SECRET-ACCESS-KEY'];
    delete env.STATE_snapshot;
    if (!post) {
      env.AWS_ACCESS_KEY_ID = '';
      env.AWS_SECRET_ACCESS_KEY = '';
    }
    const result = spawnSync('uv', [
      'run', '--locked', '--no-project', '--no-build',
      resolve(__dirname, '../../../CI/ccache_snapshot.py'), operation,
    ], {
      env,
      stdio: 'inherit', timeout: 600000, killSignal: 'SIGKILL',
    });
    if (result.error || result.status !== 0) {
      throw new Error('Transfer process failed');
    }
  }
} catch {
  // Cache failures are nonfatal, and diagnostics must not expose saved credentials.
  console.log(`::warning::ccache snapshot ${operation} failed; continuing without remote cache updates`);
}
