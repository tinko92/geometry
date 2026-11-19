// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// Prepare without compilation: node test/strategies/run_exact_redesign_issues.js --prepare OUT
// After the parent freeze: node test/strategies/run_exact_redesign_issues.js --run OUT
// Optional EXACT_NEAREST_STRATEGY overrides the named nearest_output strategy type.
'use strict';
const fs = require('fs'), path = require('path'), crypto = require('crypto');
const cp = require('child_process'), assert = require('assert');
const repo = path.resolve(__dirname, '../..');
const action = process.argv[2] || '--help';
if (action === '--help') {
    console.log('Use --prepare OUT, then --run OUT only after the parent header freeze.');
    console.log('The final build runs both approximate_output and nearest_output (floating only).');
    process.exit(0);
}
assert(['--prepare', '--run'].includes(action), 'Expected --prepare or --run');
assert(process.argv[3], 'Specify a fresh output directory');
const out = path.resolve(process.argv[3]);
const boost = path.resolve(process.env.BOOST_ROOT || '/Users/tinko/dev/boost');
const cxx = process.env.CXX || 'clang++';
const sha = data => crypto.createHash('sha256').update(data).digest('hex');
const git = args => cp.execFileSync('git', args, {cwd: repo, maxBuffer: 256 * 1024 * 1024});
const text = file => fs.readFileSync(file, 'utf8');
const writeJSON = (file, value) => fs.writeFileSync(file, JSON.stringify(value, null, 2) + '\n');
const manifestFile = path.join(out, 'manifest.json');
const revisions = Object.fromEntries(Object.entries({before: '2886f2acd', pr1491: '72b2d832', develop: 'aaf0560'})
    .map(([name, revision]) => [name, git(['rev-parse', revision]).toString().trim()]));
const driver = 'test/strategies/exact_issue_comparison.cpp';
const evaluatorSource = 'test/strategies/exact_issue_comparison_evaluate.cpp';
const immutable = [evaluatorSource, 'test/strategies/exact_issue_cases.hpp',
    'test/strategies/exact_additional_issue_cases.hpp', 'test/algorithms/overlay/overlay_cases.hpp',
    'test/algorithms/overlay/multi_overlay_cases.hpp'];
const units = [[1, 'double'], [1, 'float'], [1, 'int32_t'], [1, 'int64_t'], [2, ''], [3, ''], [4, '']];
const unitName = ([group, coordinate]) => `g${group}${coordinate ? '-' + coordinate : ''}`;
const flags = ['-std=c++17', '-O0', '-g0', '-DNDEBUG', '-ffp-contract=off', '-pthread',
    '-Wno-constant-conversion'];
const limits = {compile: 600, before: 90, after: 90, after_nearest: 90, develop: 30, pr1491: 30, evaluate: 45};

function stableFile(file, data) {
    fs.mkdirSync(path.dirname(file), {recursive: true});
    if (fs.existsSync(file)) assert.equal(sha(fs.readFileSync(file)), sha(data), `Frozen file changed: ${file}`);
    else fs.writeFileSync(file, data);
}
function headerManifest(directory) {
    const entries = [];
    function visit(relative) {
        for (const name of fs.readdirSync(path.join(directory, relative)).sort()) {
            const file = path.join(relative, name), stat = fs.lstatSync(path.join(directory, file));
            assert(!stat.isSymbolicLink(), `Unexpected header symlink: ${file}`);
            if (stat.isDirectory()) visit(file);
            else if (stat.isFile()) entries.push(`${file}\t${sha(fs.readFileSync(path.join(directory, file)))}`);
        }
    }
    visit('');
    const contents = entries.sort().join('\n') + '\n';
    return {files: entries.length, sha256: sha(contents), contents};
}
function exportRevision(name, revision) {
    const destination = path.join(out, 'headers', name);
    fs.mkdirSync(destination, {recursive: true});
    if (!fs.existsSync(path.join(destination, 'include'))) {
        const archive = git(['archive', '--format=tar', revision, 'include']);
        const result = cp.spawnSync('tar', ['-xf', '-', '-C', destination], {input: archive});
        assert.equal(result.status, 0, `Header extraction failed: ${name}`);
    }
    const rows = git(['ls-tree', '-rz', revision, 'include']).toString().split('\0').filter(Boolean);
    for (const row of rows) {
        const tab = row.indexOf('\t'), [mode, kind, expected] = row.slice(0, tab).split(' ');
        const file = row.slice(tab + 1);
        assert.equal(kind, 'blob'); assert(['100644', '100755'].includes(mode), `Unexpected file mode ${mode}`);
        const bytes = fs.readFileSync(path.join(destination, file));
        const actual = crypto.createHash(expected.length === 40 ? 'sha1' : 'sha256')
            .update(`blob ${bytes.length}\0`).update(bytes).digest('hex');
        assert.equal(actual, expected, `Git blob mismatch: ${name}/${file}`);
    }
    const headers = headerManifest(path.join(destination, 'include'));
    assert.equal(headers.files, rows.length, `Untracked files in frozen ${name} headers`);
    stableFile(path.join(out, `${name}-headers.tsv`), headers.contents);
    return {files: headers.files, sha256: headers.sha256, revision};
}
function snapshotSources() {
    const hashes = {};
    for (const file of [driver, ...immutable]) {
        const bytes = fs.readFileSync(path.join(repo, file));
        if (immutable.includes(file))
            assert.equal(sha(bytes), sha(git(['show', `${revisions.before}:${file}`])), `Original fixture/evaluator changed: ${file}`);
        stableFile(path.join(out, 'sources', file), bytes);
        hashes[file] = sha(bytes);
    }
    stableFile(path.join(out, 'original-driver.cpp'), git(['show', `${revisions.before}:${driver}`]));
    const patch = cp.spawnSync('git', ['diff', '--no-index', '--', path.join(out, 'original-driver.cpp'),
        path.join(out, 'sources', driver)], {encoding: 'utf8'});
    assert(patch.status === 0 || patch.status === 1, 'Cannot record driver adaptation');
    stableFile(path.join(out, 'driver-adaptation.diff'), patch.stdout);
    stableFile(path.join(out, 'runner.js'), fs.readFileSync(__filename));
    return hashes;
}
function verifyInputs(manifest) {
    for (const [file, hash] of Object.entries(manifest.source_sha256)) {
        assert.equal(sha(fs.readFileSync(path.join(repo, file))), hash, `Live fixture/driver changed: ${file}`);
        assert.equal(sha(fs.readFileSync(path.join(out, 'sources', file))), hash, `Frozen source changed: ${file}`);
    }
    assert.equal(sha(fs.readFileSync(__filename)), manifest.runner_sha256, 'Runner changed after preparation');
}

fs.mkdirSync(out, {recursive: true});
if (action === '--prepare') {
    assert(!fs.existsSync(manifestFile), 'Output already prepared; retain it and choose a fresh directory');
    const manifest = {state: 'preparing', repo, boost, revisions, flags, limits, compile_jobs: 1,
        runner_sha256: sha(fs.readFileSync(__filename)), prepared: new Date().toISOString(), headers: {},
        contract: 'Original fixtures and evaluator; typed operation scopes surround each complete fixture action. ' +
            'Approximate output is separate from explicit nearest output. Integral conversion stays BG-compatible.'};
    manifest.source_sha256 = snapshotSources();
    for (const [name, revision] of Object.entries(revisions)) manifest.headers[name] = exportRevision(name, revision);
    manifest.state = 'prepared';
    writeJSON(manifestFile, manifest);
    console.log(`Prepared immutable inputs and three baseline trees in ${out}; no compilation or issue execution.`);
    process.exit(0);
}

const manifest = JSON.parse(text(manifestFile));
assert.equal(manifest.state, 'prepared', 'Do not overwrite an earlier run; prepare a fresh output directory');
assert.equal(boost, manifest.boost, 'BOOST_ROOT changed after preparation');
assert.deepEqual(revisions, manifest.revisions);
verifyInputs(manifest);
const nearestOverride = process.env.EXACT_NEAREST_STRATEGY;
const nearest = nearestOverride || 'bg::strategies::relate::cartesian_exact<C, bg::strategies::relate::nearest_output>';
assert(!/[\r\n]/.test(nearest), 'Invalid strategy type expression');
manifest.nearest_strategy_adapter = nearest;
manifest.compiler = cp.execFileSync(cxx, ['--version'], {encoding: 'utf8'});
manifest.started = new Date().toISOString();
manifest.commands = []; manifest.runs = []; manifest.state = 'running';
fs.mkdirSync(path.join(out, 'runs'), {recursive: true});
fs.mkdirSync(path.join(out, 'bin'), {recursive: true});

// A copy, not a live -I path. Detect edits racing this snapshot; later live
// geometry edits cannot contaminate any baseline or final executable.
const live = path.join(repo, 'include'), start = headerManifest(live);
fs.cpSync(live, path.join(out, 'headers', 'after', 'include'), {recursive: true});
const frozen = headerManifest(path.join(out, 'headers', 'after', 'include'));
assert.equal(frozen.sha256, start.sha256, 'Geometry headers changed during freezing');
assert.equal(headerManifest(live).sha256, start.sha256, 'Geometry headers changed during freezing');
stableFile(path.join(out, 'after-headers.tsv'), frozen.contents);
manifest.headers.after = {files: frozen.files, sha256: frozen.sha256};
manifest.headers.after_nearest = manifest.headers.after;
manifest.geometry_diff_sha256 = sha(git(['diff', '--binary', 'HEAD', '--', 'include']));
fs.writeFileSync(path.join(out, 'geometry-adaptation.diff'), git(['diff', '--binary', 'HEAD', '--', 'include']));
manifest.boost_head = cp.execFileSync('git', ['rev-parse', 'HEAD'], {cwd: boost, encoding: 'utf8'}).trim();
for (const [name, args] of [['boost-submodules.txt', ['submodule', 'status']],
    ['boost-working-state.txt', ['status', '--short', '--untracked-files=no']]])
    fs.writeFileSync(path.join(out, name), cp.execFileSync('git', args, {cwd: boost, encoding: 'utf8'}));
manifest.external_dependencies = 'Other Boost headers come from recorded BOOST_ROOT; compiler dependency files are retained. ' +
    'Geometry headers, fixtures, evaluator and driver are fully snapshotted.';
writeJSON(manifestFile, manifest);

function bounded(command, args, seconds, input) {
    return cp.spawnSync('perl', ['-e', 'alarm shift; exec @ARGV or die $!', String(seconds), command, ...args],
        {cwd: out, encoding: 'utf8', input, timeout: (seconds + 5) * 1000, maxBuffer: 64 * 1024 * 1024});
}
function compile(build, unit, evaluator = false) {
    const bin = path.join(out, 'bin', evaluator ? 'evaluate' : `${build}-${unitName(unit)}`);
    const headerBuild = build === 'after_nearest' ? 'after' : build;
    const args = [...flags, `-I${out}/headers/${headerBuild}/include`, `-I${boost}`, '-MMD', '-MF', bin + '.d'];
    if (!evaluator) {
        args.push(`-DCOMPARISON_GROUP=${unit[0]}`);
        if (['before', 'after', 'after_nearest'].includes(build)) args.push('-DCOMPARISON_EXACT');
        if (build === 'after_nearest') {
            args.push('-DCOMPARISON_NEAREST_OUTPUT');
            if (nearestOverride) args.push(`-DCOMPARISON_EXACT_STRATEGY(C)=${nearestOverride}`);
        }
        if (unit[1]) args.push(`-DCOMPARISON_COORDINATE=${unit[1].startsWith('int') ? 'std::' : ''}${unit[1]}`,
            `-DCOMPARISON_COORDINATE_NAME="${unit[1]}"`);
    }
    args.push(path.join(out, 'sources', evaluator ? evaluatorSource : driver), '-o', bin);
    console.log(`Compile ${path.basename(bin)}`);
    const result = bounded(cxx, args, limits.compile);
    fs.writeFileSync(bin + '.compile.log', (result.stdout || '') + (result.stderr || ''));
    manifest.commands.push({command: cxx, args, status: result.status, signal: result.signal});
    writeJSON(manifestFile, manifest);
    return result.status === 0 ? bin : null;
}
function list(bin) {
    assert(bin, 'The frozen-before build must compile to establish the unchanged fixture catalog');
    assert.equal(bounded(bin, ['--binary-check'], 10).status, 0, 'Binary input lifting invariant failed');
    const result = bounded(bin, ['--list'], 10);
    assert.equal(result.status, 0, 'Fixture listing failed');
    return result.stdout.trim().split('\n').filter(Boolean).map(row => {
        const [name, issue, coordinate] = row.split('\t');
        return {case: name, issue: Number(issue), coordinate};
    }).filter(item => item.issue !== 494);
}
const records = [], catalog = new Map();
const evaluator = compile('before', [0, ''], true);
assert(evaluator, 'Frozen evaluator compilation failed');
const evaluatorHash = sha(fs.readFileSync(evaluator));
manifest.evaluator_sha256 = evaluatorHash;
function runCase(build, bin, item, mode) {
    const stem = path.join(out, 'runs', `${build}-${item.case}-${mode}`);
    const result = bin ? bounded(bin, [item.case, mode], limits[build]) : {status: null, signal: null};
    fs.writeFileSync(stem + '.stdout', result.stdout || '');
    fs.writeFileSync(stem + '.stderr', result.stderr || '');
    let value;
    try { value = JSON.parse(result.stdout); }
    catch { value = {...item, status: !bin ? 'compile_failure' : result.signal ? 'timeout_or_signal' : 'process_failure'}; }
    if (bin && result.status === 0) {
        assert.equal(value.case, item.case, 'Wrong fixture executed');
        assert.equal(Number(value.issue), item.issue, 'Wrong issue executed');
        assert.equal(value.coordinate, item.coordinate, 'Wrong coordinate type executed');
    }
    value.build = build; value.mode = mode;
    value.output_policy = build === 'after' ? 'approximate' : build === 'after_nearest' ? 'nearest_floating_bg_integer'
        : build === 'before' ? 'legacy_nearest_including_integer' : 'default';
    value.process = {status: result.status, signal: result.signal, error: result.error?.message};
    const evaluated = bounded(evaluator, [], limits.evaluate, JSON.stringify(value));
    fs.writeFileSync(stem + '.evaluation.stderr', evaluated.stderr || '');
    try { value = JSON.parse(evaluated.stdout); }
    catch { value.evaluation = 'process_failure'; }
    value.evaluation_process = {status: evaluated.status, signal: evaluated.signal, error: evaluated.error?.message};
    fs.writeFileSync(stem + '.json', JSON.stringify(value, null, 2) + '\n');
    records.push(value);
    fs.appendFileSync(path.join(out, 'results.jsonl'), JSON.stringify(value) + '\n');
    manifest.runs.push({build, case: item.case, mode, process: value.process, evaluation: value.evaluation,
        evaluation_process: value.evaluation_process});
    writeJSON(manifestFile, manifest);
    console.log(`${build} ${item.case} ${mode}: ${value.status}/${value.evaluation}`);
}
try {
    for (const build of ['before', 'develop', 'pr1491', 'after', 'after_nearest']) {
        for (const unit of units) {
            const bin = compile(build, unit);
            const cases = build === 'before' ? list(bin) : catalog.get(unitName(unit));
            if (build === 'before') catalog.set(unitName(unit), cases);
            else if (bin) assert.deepEqual(list(bin), cases, `Fixture selection changed: ${build}/${unitName(unit)}`);
            for (const item of cases) for (const mode of ['before', 'after', 'after_nearest'].includes(build)
                ? ['native', 'retained'] : ['native']) runCase(build, bin, item, mode);
        }
        const items = [...catalog.values()].flat(), subset = records.filter(r => r.build === build);
        assert.equal(new Set(items.map(item => item.case)).size, 82);
        assert.equal(new Set(items.map(item => item.issue)).size, 47);
        assert.equal(new Set(subset.map(r => `${r.case}/${r.mode}`)).size, subset.length);
        assert.equal(subset.length, ['before', 'after', 'after_nearest'].includes(build) ? 164 : 82);
        manifest[build + '_counts'] = {records: subset.length, fixtures: 82, issues: 47};
    }
    verifyInputs(manifest);
    for (const [build, expected] of Object.entries(manifest.headers)) {
        const location = build === 'after_nearest' ? 'after' : build;
        assert.equal(headerManifest(path.join(out, 'headers', location, 'include')).sha256, expected.sha256,
            `Frozen headers changed: ${build}`);
    }
    assert.equal(sha(fs.readFileSync(evaluator)), evaluatorHash, 'Frozen evaluator binary changed');
    function differences(a, b, prefix = '') {
        if (a === b) return [];
        if (a && b && typeof a === 'object' && typeof b === 'object')
            return [...new Set([...Object.keys(a), ...Object.keys(b)])].sort()
                .flatMap(key => differences(a[key], b[key], prefix ? prefix + '.' + key : key));
        return [{path: prefix, before: a, after: b}];
    }
    const semantic = r => Object.fromEntries(Object.entries(r).filter(([key]) =>
        !['build', 'mode', 'process', 'evaluation_process', 'output_policy'].includes(key)));
    const previous = new Map(records.filter(r => r.build === 'before').map(r => [`${r.case}/${r.mode}`, r]));
    const changes = records.filter(r => ['after', 'after_nearest'].includes(r.build)).map(r => ({
        build: r.build, case: r.case, mode: r.mode,
        differences: differences(semantic(previous.get(`${r.case}/${r.mode}`)), semantic(r))
    })).filter(r => r.differences.length);
    const failures = records.filter(r => r.status !== 'ok' || r.evaluation !== 'ok'
        || String(r.process.status) !== '0' || r.evaluation_process.status !== 0);
    const isFalse = value => value === false || value === 'false';
    const inputField = key => /^(a|b|input[0-9]+)$/.test(key);
    const independent = new Set(['boundary_oracle', 'linear_oracle', 'rational_point', 'top_chain_area']);
    const assessments = records.map(r => {
        const inputs = Object.entries(r.geometries || {}).filter(([key]) => inputField(key));
        const invalidInputs = inputs.filter(([, g]) => isFalse(g.rational_valid)).map(([key]) => key);
        // The 669 control deliberately feeds the default library's crop back
        // into the tested algorithm. It is not a native-input exactness promise.
        if (r.case.endsWith('source_default_crop') && isFalse(r.geometries?.source?.rational_valid))
            invalidInputs.push('source_default_crop');
        const outputs = Object.entries(r.geometries || {}).filter(([key]) => !inputField(key)
            && key !== 'source' && !(Number(r.issue) === 1174 && key === 'result'));
        const invalidOutputs = outputs.filter(([, g]) => isFalse(g.rational_valid)).map(([key]) => key);
        const publicDiagnostics = outputs.filter(([, g]) => isFalse(g.public_valid)).map(([key]) => key);
        const failedChecks = Object.entries(r.checks || {}).filter(([, value]) => isFalse(value)).map(([key]) => key);
        const hasExactObservation = invalidOutputs.length > 0 || failedChecks.length > 0;
        const transportFailed = failures.includes(r);
        return {build: r.build, case: r.case, issue: r.issue, mode: r.mode, output_policy: r.output_policy,
            classification: r.mode === 'native' && r.exception_kind === 'numeric_conversion'
                ? 'native_conversion_exception' : transportFailed ? 'execution_or_evaluator_failure' : invalidInputs.length
                ? 'out_of_contract_diagnostic_input' : r.mode === 'retained'
                    ? hasExactObservation ? 'internal_exact_observation' : 'internal_exact_checks_pass'
                    : hasExactObservation ? 'native_representation_or_feedback_observation' : 'native_checks_pass',
            invalid_inputs: invalidInputs, rational_invalid_outputs: invalidOutputs,
            failed_independent_oracles: failedChecks.filter(key => independent.has(key)),
            failed_shared_library_checks: failedChecks.filter(key => !independent.has(key)),
            public_validator_diagnostics: publicDiagnostics};
    });
    const internal = assessments.filter(r => r.classification === 'internal_exact_observation');
    const representations = assessments.filter(r => ['native_representation_or_feedback_observation',
        'native_conversion_exception'].includes(r.classification));
    const currentInternal = internal.filter(r => ['after', 'after_nearest'].includes(r.build));
    writeJSON(path.join(out, 'assessments.json'), {records: assessments,
        internal_exact_observations: internal.length, native_representation_observations: representations.length,
        interpretation: 'The unchanged independent rational oracles are authoritative for retained exact results. ' +
            'Frozen rational validity and other shared-library checks are separately identified, not complete independent set goldens. ' +
            'Native-output validity/set-equality failures may follow permitted approximation, collapse, or feedback. ' +
            'No native representation observation alone proves an internal topology defect or an invalid conversion. ' +
            'Stored-approximation emission/numeric_cast compliance is checked by separate output-policy primitive tests, ' +
            'not reconstructed from post-validation retained records. No legacy integer ties-to-even equivalence is required.'});
    const inputDigests = new Map(), inputMismatches = [];
    for (const r of records.filter(r => r.status === 'ok' && r.geometries)) {
        const inputs = Object.keys(r.geometries).filter(inputField).sort()
            .map(key => [key, r.geometries[key].wkt]);
        const digest = sha(JSON.stringify(inputs)), old = inputDigests.get(r.case);
        if (old && old.sha256 !== digest) inputMismatches.push({build: r.build, case: r.case, mode: r.mode});
        else if (!old) inputDigests.set(r.case, {sha256: digest, fields: inputs.map(([key]) => key)});
    }
    writeJSON(path.join(out, 'input-digests.json'), Object.fromEntries(inputDigests));
    writeJSON(path.join(out, 'comparison.json'), {records: records.length, failures, changes, input_mismatches: inputMismatches,
        policy: 'Every raw field difference is retained. Approximate/native rounding, integer conversion, ' +
            'public-validator diagnostics and ring ordering may differ intentionally; a raw difference is not a topology bug.'});
    manifest.state = failures.length || inputMismatches.length ? 'completed_with_failures'
        : currentInternal.length ? 'complete_with_internal_exact_observations'
        : representations.length ? 'complete_with_representation_observations' : 'complete';
    manifest.finished = new Date().toISOString(); manifest.records = records.length;
    manifest.failures = failures.length; manifest.changed_records = changes.length;
    manifest.input_mismatches = inputMismatches.length;
    manifest.internal_exact_observations = internal.length;
    manifest.current_internal_exact_observations = currentInternal.length;
    manifest.native_representation_observations = representations.length;
    writeJSON(manifestFile, manifest);
    console.log(JSON.stringify({state: manifest.state, records: records.length, failures: failures.length,
        changed_records: changes.length, current_internal_exact_observations: currentInternal.length,
        native_representation_observations: representations.length}, null, 2));
    if (failures.length || inputMismatches.length || currentInternal.length) process.exitCode = 1;
} catch (error) {
    manifest.state = 'failed'; manifest.error = error.stack; writeJSON(manifestFile, manifest);
    throw error;
}
