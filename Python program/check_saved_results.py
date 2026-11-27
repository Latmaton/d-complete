"""Check only the consistency of saved records, without replacing the construction of residue families and terminal lower bounds from scratch."""
import argparse
import csv
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import re
from zipfile import ZipFile

from parameters import (PARAMETERS, SCALES, INTERVALS, PAPER_TARGETS,
                        HEIGHT_LIMITS, POSITIVE_SETTINGS)

ROOT = Path(__file__).resolve().parent
CORE_FILES = ('main.py', 'parameters.py', 'residues.py',
              'terminals.py', 'arithmetic.py', 'reporting.py')
STATUSES = {'passed', 'running', 'failed', 'interrupted'}
STAGES = {'residues', 'tails', 'terminals', 'density'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value, name, minimum=None):
    require(not isinstance(value, bool), name + ' must not be Boolean')
    if isinstance(value, int):
        result = value
    elif isinstance(value, str) and re.fullmatch(r'-?\d+', value):
        result = int(value)
    else:
        raise ValueError(name + ' must be an exact integer')
    if minimum is not None:
        require(result >= minimum, name + ' is below the permitted range')
    return result


def rational(value, name):
    require(not isinstance(value, bool), name + ' must not be Boolean')
    if isinstance(value, int):
        return F(value)
    if isinstance(value, str) and re.fullmatch(r'-?\d+(?:/[1-9]\d*)?', value):
        return F(value)
    raise ValueError(name + ' must be an integer or an exact fraction string, not a floating-point number')


def json_form(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: json_form(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_form(item) for item in value]
    return value


def check_json_types_and_values(actual, expected, name):
    """Compare JSON data without allowing Boolean/integer/float coercion."""
    require(type(actual) is type(expected), name + ' has an incorrect JSON type')
    if isinstance(expected, dict):
        require(set(actual) == set(expected), name + ' has missing or extra fields')
        for key, value in expected.items():
            check_json_types_and_values(actual[key], value, name + '.' + key)
    elif isinstance(expected, list):
        require(len(actual) == len(expected), name + ' has an incorrect length')
        for i, (item, value) in enumerate(zip(actual, expected)):
            check_json_types_and_values(item, value, name + '[' + str(i) + ']')
    else:
        require(actual == expected, name + ' differs from the current fixed parameters')


def check_summary_header(summary):
    require(isinstance(summary, dict), 'summary.json must contain a JSON object')
    status = summary.get('status')
    require(isinstance(status, str) and status in STATUSES,
            'summary.status is missing or is not one of passed, running, failed, interrupted')
    pair = summary.get('pair')
    require(isinstance(pair, list) and len(pair) == 2
            and all(type(value) is int for value in pair),
            'summary.pair must contain two JSON integers')
    pair = tuple(pair)
    require(pair in PARAMETERS, 'The record does not belong to the eight fixed cases')
    return pair, status


def check_completed_stages(summary):
    timings = summary.get('stages_seconds')
    require(isinstance(timings, dict) and set(timings) == STAGES,
            'A passed record must contain timings for all four stages')
    for stage, seconds in timings.items():
        require(type(seconds) in (int, float) and seconds >= 0
                and (isinstance(seconds, int) or math.isfinite(seconds)),
                'The timing for ' + stage + ' must be a finite nonnegative number')


def unique_object(items):
    result = {}
    for key, value in items:
        require(key not in result, 'Duplicate JSON key: ' + key)
        result[key] = value
    return result


def read_json(path):
    with path.open(encoding='utf-8') as handle:
        return json.load(handle, object_pairs_hook=unique_object)


def read_csv(path, fields):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == fields, path.name + ' has incorrect column names or column order')
        rows = list(reader)
    require(all(None not in row and all(v is not None for v in row.values())
                for row in rows), path.name + ' has an incorrect row length')
    return rows


def check_parameters(folder, summary):
    require(isinstance(summary.get('pair'), list) and len(summary['pair']) == 2,
            'summary.pair must contain two integers')
    pair = tuple(integer(n, 'bases') for n in summary['pair'])
    require(pair in PARAMETERS, 'The record does not belong to the eight fixed cases')
    d, B, ell, z, C, m0, tmax, J, target = PARAMETERS[pair]
    require(rational(summary['paper_target'], 'paper target') == PAPER_TARGETS[pair],
            'The paper target in summary differs from parameters.py')
    require(rational(summary['verification_target'], 'internal target') == target,
            'The internal target in summary differs from parameters.py')
    saved = read_json(folder / 'parameters.json')
    chunk = integer(saved['parameters']['chunk'], 'chunk', 1)
    p, q = pair
    expected_params = dict(p=p, q=q, d=d, B=B, ell=ell, z=F(z), C=F(C),
                           m0=m0, t_max=tmax, J=J, target=F(target),
                           paper_target=F(PAPER_TARGETS[pair]),
                           D0=POSITIVE_SETTINGS['denominator'], P=abs(p)**d,
                           s=len(B), H_MAX=HEIGHT_LIMITS[pair], chunk=chunk)
    expected_scales = [dict(H=H, X=F(X), n=n, count_exponent=e,
                            tail_terms=min(H, tmax))
                       for H, X, n, e in SCALES[pair]]
    expected_intervals = [(F(left), F(right), choices)
                          for left, right, choices in INTERVALS[pair]]
    expected = dict(parameters=expected_params, scales=expected_scales,
                    intervals=expected_intervals,
                    positive_settings=POSITIVE_SETTINGS if pair == (5, 2) else None)
    check_json_types_and_values(saved, json_form(expected), 'parameters.json')
    for key, count in [('scale_count', len(expected_scales)),
                       ('interval_count', len(expected_intervals)),
                       ('endpoint_count', 2 * len(expected_intervals))]:
        require(integer(summary[key], key, 0) == count, key + ' differs')
    return pair, expected_params, expected_scales


def check_constants(folder, params):
    p, q, P = params['p'], params['q'], params['P']
    B, ell, s = params['B'], params['ell'], params['s']
    c = F(abs(q)**sum(ell), P**s)
    require(c < 1, 'The contraction constant is not less than 1')
    A = sum((F(B[i] * abs(q)**(sum(ell[i:]) - 1), P**(s-i))
             for i in range(s)), F(0)) / (1-c)
    kappa = 1 if p > 0 and q > 0 else 2
    order = next(e for e in range(1, abs(q)+1) if pow(p, e, abs(q)) == 1)
    theta = F(kappa * order, abs(q))
    expected = dict(c=c, A=A, kappa=kappa, theta=theta)
    saved = read_json(folder / 'constants.json')
    require(set(saved) == set(expected), 'constants.json has incomplete fields')
    for key, value in expected.items():
        require(rational(saved[key], key) == value, 'Constant ' + key + ' differs')
    return A, kappa, theta


def check_histograms(folder, params):
    saved = read_json(folder / 'residue_heights.json')
    hist = {}
    for key, row in saved.items():
        B = integer(key, 'size bound', 1)
        require(B not in hist, 'Duplicate size bound')
        hist[B] = {}
        for height, number in row.items():
            h = integer(height, 'family height', 0)
            require(h not in hist[B], 'Duplicate family height')
            require(h <= params['H_MAX'], 'The family height exceeds the search range')
            hist[B][h] = integer(number, 'height count', 0)
        require(hist[B].get(0) == 1 and sum(hist[B].values()) == params['P'],
                'The height-zero count or total coverage count differs')
    require(set(hist) == set(params['B']), 'The height distribution has incomplete size bounds')
    rows = read_csv(folder / 'residue_heights.csv', ['B', 'height', 'count'])
    csv_values = [(integer(row['B'], 'B'), integer(row['height'], 'height'),
                   integer(row['count'], 'count')) for row in rows]
    expected = [(B, h, number) for B in sorted(hist)
                for h, number in sorted(hist[B].items())]
    require(csv_values == expected, 'The residue-height CSV and JSON disagree')
    return hist


def check_tails(folder, params, hist):
    P, ell, z = params['P'], params['ell'], params['z']
    C, D0, s = params['C'], params['D0'], params['s']
    moments = [sum((number * z**(h-ell[i]) for h, number in hist[B].items()), F(0))/P
               for i, B in enumerate(params['B'])]
    g = F(1)
    for value in moments:
        g *= value
    partial, Lambda = F(1), F(1)
    for value in reversed(moments[1:]):
        partial *= value
        Lambda = max(Lambda, partial)
    require(g < 1 and Lambda <= C, 'The exponential-moment conditions do not hold')
    saved = read_json(folder / 'moments.json')
    require(set(saved) == {'g', 'Lambda'}, 'moments.json has incomplete fields')
    require(rational(saved['g'], 'g') == g and rational(saved['Lambda'], 'Lambda') == Lambda,
            'The saved exponential-moment product or Lambda differs')

    # Recompute finite tail bounds from the recorded height distributions without reconstructing residue representations.
    length = params['t_max'] + params['m0'] * sum(ell)
    w = [D0]
    for t in range(1, length+1):
        bound = D0 * min(F(1), C/z**t)
        w.append(-(-bound.numerator // bound.denominator))
    for step in range(params['m0'] * s):
        i = step % s
        length -= ell[i]
        new = [D0]
        for t in range(1, length+1):
            total = sum(number * (D0 if t+ell[i]-h <= 0 else w[t+ell[i]-h])
                        for h, number in hist[params['B'][i]].items())
            new.append(-(-total // P))
        w = new
    expected = [F(number, D0) for number in w]
    rows = read_csv(folder / 'tails.csv', ['t', 'upper_bound'])
    actual = [(integer(row['t'], 'tail-bound index', 0), rational(row['upper_bound'], 'tail bound'))
              for row in rows]
    require(actual == list(enumerate(expected)), 'The tail-bound CSV differs from recomputation using the height distributions')
    require(expected[0] == 1 and all(0 <= b <= a <= 1 for a, b in zip(expected, expected[1:])),
            'The tail bounds are not monotone or lie outside the allowed range')
    return expected


def check_terminals(folder, summary, scales):
    requests = {(row['H']-v, row['X']) for row in scales
                for v in range(row['tail_terms'])}
    rows = read_csv(folder / 'terminal_lower_bounds.csv', ['h', 'X', 'lower_bound'])
    terminal = {}
    for row in rows:
        key = integer(row['h'], 'terminal height', 1), rational(row['X'], 'terminal radius')
        require(key not in terminal, 'Duplicate terminal request')
        terminal[key] = rational(row['lower_bound'], 'terminal lower bound')
        require(terminal[key] >= 0, 'A terminal lower bound is negative')
    require(set(terminal) == requests, 'The terminal records have missing or extra requests')
    require(integer(summary['terminal_request_count'], 'terminal request count', 0) == len(requests),
            'The terminal request count in summary differs')
    return terminal


def check_scales(folder, params, sources, terminal, w, A):
    fields = ['index', 'H', 'X', 'n', 'count_exponent', 'tail_terms', 'a', 'b']
    rows = read_csv(folder / 'scales.csv', fields)
    require(len(rows) == len(sources), 'The scale CSV has an incorrect number of rows')
    result = []
    D0, p, q = params['D0'], abs(params['p']), abs(params['q'])
    for i, (row, source) in enumerate(zip(rows, sources), 1):
        require(integer(row['index'], 'scale index') == i, 'The scale index differs')
        for key, value in source.items():
            require(rational(row[key], 'Scale field ' + key) == value, 'Scale field ' + key + ' differs')
        H, X, n = source['H'], source['X'], source['n']
        boundary = (A*q**H+X)/p**n
        scaled = D0 * boundary
        a = F(scaled.numerator // scaled.denominator + 1, D0)
        beta = sum(((w[v]-w[v+1]) * terminal[H-v, X]
                    for v in range(source['tail_terms'])), F(0))
        scaled = D0**2 * beta / p**n
        b = F(scaled.numerator // scaled.denominator, D0**2)
        require(a > boundary and b > 0, 'The strict scale inequalities do not hold')
        require(rational(row['a'], 'a') == a and rational(row['b'], 'b') == b,
                'The scale coefficients disagree with the terminal lower bounds and tail bounds')
        result.append((a, b))
    return result


def check_endpoints(folder, summary, pair, params, scales, kappa, theta):
    records = read_json(folder / 'endpoints.json')
    intervals = INTERVALS[pair]
    require(len(records) == 2 * len(intervals), 'The endpoint JSON has an incorrect number of rows')
    csv_rows = read_csv(folder / 'endpoints.csv', ['interval_index', 'left', 'right', 'y', 'value'])
    require(len(csv_rows) == len(records), 'The endpoint CSV has an incorrect number of rows')
    minimum = None
    for v, (left, right, choices) in enumerate(intervals, 1):
        for side, y in enumerate((F(left), F(right))):
            index = 2*(v-1)+side
            record, csv_row = records[index], csv_rows[index]
            require(integer(record['interval_index'], 'interval index') == v, 'The interval index differs')
            for key, value in [('left', F(left)), ('right', F(right)), ('y', y)]:
                require(rational(record[key], key) == value, 'An endpoint or partition interval differs from the parameters')
            contributions = record['contributions']
            require(len(contributions) == params['J'], 'The number of individual endpoint contributions differs')
            total = F(0)
            for j, ((i, k), part) in enumerate(zip(choices, contributions)):
                for key, value in [('j', j), ('i', i), ('k', k)]:
                    require(integer(part[key], key) == value, 'The endpoint choices differ from the parameters')
                a, b = scales[i-1]
                factor = F(abs(pair[0]))**k
                penalty = theta * max(F(0), factor*a-y/abs(pair[1])**j)
                contribution = factor*b-penalty
                expected = dict(a=a, b=b, factor=factor, penalty=penalty, contribution=contribution)
                for key, value in expected.items():
                    require(rational(part[key], key) == value, 'Individual endpoint field ' + key + ' differs')
                total += contribution
            value = total / (kappa * (1-F(1, abs(pair[1])**params['J'])) * y)
            require(rational(record['numerator'], 'endpoint numerator') == total, 'The endpoint numerator does not equal the sum of individual contributions')
            require(rational(record['value'], 'endpoint value') == value, 'The endpoint denominator normalization differs')
            for key in ('interval_index', 'left', 'right', 'y', 'value'):
                require(rational(csv_row[key], key) == rational(record[key], key),
                        'The endpoint CSV and JSON disagree')
            minimum = value if minimum is None else min(minimum, value)
    require(rational(summary['lower_bound'], 'final lower bound') == minimum, 'The lower bound in summary is not the minimum endpoint value')
    require(minimum > params['target'] >= params['paper_target'], 'The final lower bound does not strictly exceed the target')
    return minimum


def check_case(folder, summary, hashes):
    check_completed_stages(summary)
    require(summary['source_sha256'] == hashes, 'The hashes of the six core source files differ from the checked version')
    pair, params, sources = check_parameters(folder, summary)
    A, kappa, theta = check_constants(folder, params)
    hist = check_histograms(folder, params)
    w = check_tails(folder, params, hist)
    terminal = check_terminals(folder, summary, sources)
    scales = check_scales(folder, params, sources, terminal, w, A)
    lower = check_endpoints(folder, summary, pair, params, scales, kappa, theta)
    return pair, lower


def audit(results, source_archive=None, require_all=False, pair=None):
    require(not (require_all and pair is not None), '--require-all and --pair cannot be combined')
    required = set(PARAMETERS) if require_all else ({pair} if pair is not None else set())
    if source_archive is None:
        hashes = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in CORE_FILES}
        print('Checking source hashes against the current package:', ROOT)
    else:
        source_archive = Path(source_archive)
        with ZipFile(source_archive) as archive:
            members = archive.namelist()
            for name in CORE_FILES:
                require(members.count(name) == 1,
                        'The source archive must contain exactly one member named ' + name)
            hashes = {name: hashlib.sha256(archive.read(name)).hexdigest()
                      for name in CORE_FILES}
        print('Checking source hashes against the explicitly selected archive:', source_archive)
    checked = skipped = errors = outside = 0
    covered = set()
    folders = sorted(path for path in results.iterdir() if path.is_dir())
    for folder in folders:
        try:
            path = folder / 'summary.json'
            require(path.is_file(), 'The run directory has no summary.json')
            summary = read_json(path)
            record_pair, status = check_summary_header(summary)
            if pair is not None and record_pair != pair:
                outside += 1
                print('Outside the selected pair: ', folder.name, 'bases', record_pair)
                continue
            if status != 'passed':
                skipped += 1
                print('Incomplete record: ', folder.name, 'bases', record_pair, 'status=' + status)
                continue
            record_pair, lower = check_case(folder, summary, hashes)
            checked += 1
            covered.add(record_pair)
            print('Record consistent: ', folder.name, 'bases', record_pair, 'lower bound', lower)
        except Exception as error:
            errors += 1
            print('Record inconsistent: ', folder.name, type(error).__name__ + ':', error)
    print('Consistent passed records:', checked, 'Incomplete or non-passed records:', skipped, 'Inconsistent records:', errors)
    print('Covered cases:', ', '.join(str(item) for item in PARAMETERS if item in covered) or 'none')
    missing = required - covered
    if required:
        print('Required cases:', ', '.join(str(item) for item in PARAMETERS if item in required))
        print('Missing complete cases:', ', '.join(str(item) for item in PARAMETERS if item in missing) or 'none')
    else:
        print('Full eight-case coverage was not required; use --require-all to require it.')
    if outside:
        print('Records outside the selected pair:', outside)
    if not folders:
        print('No run directories were found; no run results were checked.')
    elif checked == 0:
        print('No completed record passed the consistency checks.')
    print('Only saved records were checked; residue families and terminal lower bounds were not reconstructed.')
    if errors:
        return 1
    if skipped or missing or checked == 0:
        return 3
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT/'results',
                        help='parent directory of run results; defaults to the results directory in this package')
    parser.add_argument('--source-archive', type=Path,
                        help='explicit ZIP archive of the six original core source files; otherwise check the current package')
    coverage = parser.add_mutually_exclusive_group()
    coverage.add_argument('--require-all', action='store_true',
                          help='require complete, consistent records for all eight fixed cases')
    coverage.add_argument('--pair', nargs=2, type=int, metavar=('P', 'Q'),
                          help='check records for this pair and require complete coverage; reversed bases are accepted')
    args = parser.parse_args(argv)
    pair = tuple(args.pair) if args.pair is not None else None
    if pair is not None:
        if pair not in PARAMETERS and pair[::-1] in PARAMETERS:
            pair = pair[::-1]
        if pair not in PARAMETERS:
            parser.error('no fixed parameters for this pair')
    try:
        require(args.results.is_dir(), 'The results directory does not exist: ' + str(args.results))
        if args.source_archive is not None:
            require(args.source_archive.is_file(),
                    'The source archive does not exist: ' + str(args.source_archive))
        return audit(args.results, args.source_archive, args.require_all, pair)
    except Exception as error:
        print('Cannot check saved records:', type(error).__name__ + ':', error)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
