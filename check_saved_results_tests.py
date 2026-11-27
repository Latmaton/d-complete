"""Regression checks for saved-record status, coverage, and validation.

These tests check the record checker, not the mathematical constructions.
All damaged records are temporary copies. Supplied evidence is left unchanged.
"""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path
import shutil
import tempfile
from zipfile import ZipFile

import check_saved_results as checker


ROOT = Path(__file__).resolve().parent


def run_check(arguments, expected, description, text=None):
    output = io.StringIO()
    with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
        try:
            actual = checker.main([str(argument) for argument in arguments])
        except SystemExit as error:
            actual = error.code
    result = output.getvalue()
    if actual != expected or (text is not None and text not in result):
        raise AssertionError(description + ': expected exit ' + str(expected)
                             + ', got ' + str(actual) + '\n' + result)
    print('Passed:', description)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def current_fixture(results):
    hashes = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
              for name in checker.CORE_FILES}
    for path in sorted(results.glob('*/summary.json'), reverse=True):
        try:
            summary = checker.read_json(path)
        except (OSError, ValueError):
            continue
        if (isinstance(summary, dict) and summary.get('status') == 'passed'
                and summary.get('pair') == [-2, -3]
                and summary.get('source_sha256') == hashes):
            return path.parent
    raise ValueError('No passed (-2,-3) record matching the current six core source files '
                     'was found in ' + str(results) + '. Supply its parent directory with --results.')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT / 'results',
                        help='results parent containing a passed (-2,-3) run of the current source')
    args = parser.parse_args(argv)
    try:
        complete = current_fixture(args.results)
    except ValueError as error:
        parser.error(str(error))

    with tempfile.TemporaryDirectory(prefix='d-complete-record-tests-') as temporary:
        root = Path(temporary)

        def new_directory(name):
            directory = root / name
            directory.mkdir()
            return directory

        def copied_record(name):
            directory = new_directory(name)
            run = directory / 'run'
            shutil.copytree(complete, run)
            return directory, run

        # Isolate one actual completed run. The supplied directory may contain
        # any number of other cases or attempts and is never modified.
        current, run = copied_record('one-case')
        run_check(['--results', current], 0, 'one copied current record',
                  'Covered cases: (-2, -3)')
        run_check(['--results', current, '--require-all'], 3,
                  'one current case does not cover all eight cases', 'Missing complete cases:')
        run_check(['--results', current, '--pair', -3, -2], 0,
                  'reversed bases select the fixed case', 'Required cases: (-2, -3)')
        run_check(['--results', current, '--pair', 5, 2], 3,
                  'a missing selected pair is incomplete', 'Missing complete cases: (5, 2)')
        run_check(['--pair', 13, 17], 2, 'unsupported bases are a usage error')
        run_check(['--pair', -2, -3, '--require-all'], 2,
                  'pair selection and all-case coverage are mutually exclusive')

        empty = new_directory('empty')
        run_check(['--results', empty], 3, 'an empty results directory is incomplete')
        run_check(['--results', root / 'missing'], 1, 'a missing results directory is an error')

        missing_summary = new_directory('missing-summary')
        (missing_summary / 'run').mkdir()
        run_check(['--results', missing_summary], 1, 'a run directory without a summary is an error')

        header_cases = [
            ('missing-status', {}),
            ('unknown-status', {'status': 'passsed', 'pair': [-2, -3]}),
            ('non-object-summary', []),
            ('missing-pair', {'status': 'running'}),
            ('non-integer-pair', {'status': 'running', 'pair': [-2.0, -3]}),
        ]
        for name, summary in header_cases:
            directory = new_directory(name)
            run = directory / 'run'
            run.mkdir()
            write_json(run / 'summary.json', summary)
            run_check(['--results', directory], 1, name)

        directory, run = copied_record('invalid-json')
        (run / 'summary.json').write_text('{invalid', encoding='utf-8')
        run_check(['--results', directory], 1, 'unparseable JSON is an error')

        for status in ('running', 'failed', 'interrupted'):
            directory, run = copied_record('status-' + status)
            summary = checker.read_json(run / 'summary.json')
            summary['status'] = status
            write_json(run / 'summary.json', summary)
            run_check(['--results', directory], 3, status + ' is incomplete')

        directory, run = copied_record('complete-and-incomplete')
        unfinished = directory / 'unfinished'
        unfinished.mkdir()
        write_json(unfinished / 'summary.json', {'status': 'running', 'pair': [5, 2]})
        run_check(['--results', directory], 3, 'a complete case does not hide an incomplete record')
        run_check(['--results', directory, '--pair', -2, -3], 0,
                  'an incomplete record outside the selected pair does not block that pair')

        directory, run = copied_record('missing-stages')
        summary = checker.read_json(run / 'summary.json')
        del summary['stages_seconds']
        write_json(run / 'summary.json', summary)
        run_check(['--results', directory], 1, 'a passed record needs all stage metadata')

        directory, run = copied_record('partial-stages')
        summary = checker.read_json(run / 'summary.json')
        del summary['stages_seconds']['terminals']
        write_json(run / 'summary.json', summary)
        run_check(['--results', directory], 1, 'a passed record cannot omit the terminal stage')

        for i, value in enumerate((True, -1, float('inf'), float('nan'), '1.0')):
            directory, run = copied_record('invalid-time-' + str(i))
            summary = checker.read_json(run / 'summary.json')
            summary['stages_seconds']['density'] = value
            write_json(run / 'summary.json', summary)
            run_check(['--results', directory], 1, 'reject invalid stage time ' + repr(value))

        directory, run = copied_record('zero-stage-time')
        summary = checker.read_json(run / 'summary.json')
        summary['stages_seconds']['density'] = 0
        write_json(run / 'summary.json', summary)
        run_check(['--results', directory], 0, 'zero is a valid nonnegative stage time')

        directory, run = copied_record('float-parameter')
        parameters = checker.read_json(run / 'parameters.json')
        parameters['parameters']['p'] = float(parameters['parameters']['p'])
        write_json(run / 'parameters.json', parameters)
        run_check(['--results', directory], 1, 'same-valued float cannot replace an integer parameter')

        directory, run = copied_record('nested-float-parameter')
        parameters = checker.read_json(run / 'parameters.json')
        parameters['scales'][0]['H'] = float(parameters['scales'][0]['H'])
        write_json(run / 'parameters.json', parameters)
        run_check(['--results', directory], 1, 'nested integer parameter types are checked')

        directory, run = copied_record('boolean-parameter')
        parameters = checker.read_json(run / 'parameters.json')
        parameters['scales'][0]['n'] = bool(parameters['scales'][0]['n'])
        write_json(run / 'parameters.json', parameters)
        run_check(['--results', directory], 1, 'Boolean cannot replace an integer parameter')

        directory, run = copied_record('bad-source-hash')
        summary = checker.read_json(run / 'summary.json')
        summary['source_sha256']['main.py'] = '0' * 64
        write_json(run / 'summary.json', summary)
        run_check(['--results', directory], 1, 'source hash mismatch remains an error')

        directory, run = copied_record('missing-output')
        (run / 'tails.csv').unlink()
        run_check(['--results', directory], 1, 'a passed record cannot omit a required output')

        directory, run = copied_record('duplicate-pair')
        for i in range(7):
            shutil.copytree(complete, directory / ('duplicate-' + str(i)))
        run_check(['--results', directory, '--require-all'], 3,
                  'eight runs of one pair do not cover eight different cases')

        invalid_archive = root / 'invalid.zip'
        invalid_archive.write_text('not a ZIP archive', encoding='utf-8')
        run_check(['--results', current, '--source-archive', invalid_archive], 1,
                  'an invalid source archive is an error')
        missing_members = root / 'missing-members.zip'
        with ZipFile(missing_members, 'w') as handle:
            handle.writestr('README.txt', 'No core sources.')
        run_check(['--results', current, '--source-archive', missing_members], 1,
                  'the selected source archive must contain all six core files')

    print('All saved-record regression checks passed. No supplied records were changed.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
