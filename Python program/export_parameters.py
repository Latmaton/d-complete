"""Export all parameters from parameters.py using only the standard library and no previous computation results."""
import argparse
import csv
import hashlib
import io
from fractions import Fraction
from pathlib import Path

from parameters import (
    CHUNK, HEIGHT_LIMITS, INTERVALS, PAPER_TARGETS, PARAMETERS,
    POSITIVE_SETTINGS, SCALES, TERMINAL_LIMITS,
)

ROOT = Path(__file__).resolve().parent


def exact(value):
    """Convert integers and Fractions directly to exact text without using floating-point arithmetic."""
    if isinstance(value, (int, Fraction)):
        return str(value)
    if isinstance(value, (tuple, list)):
        return ', '.join(exact(item) for item in value)
    if isinstance(value, str):
        return value
    raise TypeError('Unsupported parameter type: ' + type(value).__name__)


def csv_text(headers, rows):
    output = io.StringIO(newline='')
    writer = csv.writer(output, lineterminator='\n')
    writer.writerow(headers)
    for row in rows:
        writer.writerow([exact(value) for value in row])
    return output.getvalue()


def markdown_table(headers, rows):
    lines = ['| ' + ' | '.join(headers) + ' |',
             '| ' + ' | '.join('---' for _ in headers) + ' |']
    for row in rows:
        lines.append('| ' + ' | '.join(exact(value) for value in row) + ' |')
    return '\n'.join(lines)


def validate_inputs():
    pairs = set(PARAMETERS)
    if not pairs or set(SCALES) != pairs or set(INTERVALS) != pairs:
        raise ValueError('Parameters, scales, and intervals must have the same base pairs')
    if set(HEIGHT_LIMITS) != pairs or set(PAPER_TARGETS) != pairs:
        raise ValueError('A base pair is missing from the height limits or paper targets')
    for pair, values in PARAMETERS.items():
        p, q = pair
        d, B, ell, z, C, m0, tmax, J, target = values
        if len(B) != len(ell) or not B or min(B) <= 0 or min(ell) <= 0:
            raise ValueError('Invalid periodic size bounds or height budgets: ' + str(pair))
        if target < PAPER_TARGETS[pair]:
            raise ValueError('The internal target is weaker than the paper target: ' + str(pair))
        if not SCALES[pair] or not INTERVALS[pair]:
            raise ValueError('The scale or interval table is empty: ' + str(pair))
        for H, X, n, e in SCALES[pair]:
            if H < 1 or X <= 0 or n < 0 or e < 0:
                raise ValueError('Invalid scale range: ' + str(pair))
        end = INTERVALS[pair][0][0]
        for left, right, choices in INTERVALS[pair]:
            if left != end or not 0 < left < right or len(choices) != J:
                raise ValueError('Intervals are not contiguous or have the wrong number of choices: ' + str(pair))
            for i, k in choices:
                if not 1 <= i <= len(SCALES[pair]) or not isinstance(k, int):
                    raise ValueError('Invalid scale index or exponent in an interval choice: ' + str(pair))
            end = right
        if end != abs(p) * INTERVALS[pair][0][0]:
            raise ValueError('The intervals do not cover a complete multiplicative period: ' + str(pair))


def terminal_mode(pair):
    if pair == (5, 2):
        return 'positive_intervals'
    if pair in TERMINAL_LIMITS:
        return 'signed_window'
    if pair in ((7, -2), (-7, -2)):
        return 'two_singletons'
    raise ValueError('No terminal algorithm is specified: ' + str(pair))


def build_files():
    """Return {relative filename: complete text} without performing any large-scale mathematical construction."""
    validate_inputs()
    files = {}
    summary_rows = []
    global_rows = [['CHUNK', CHUNK],
                   ['scale_denominator_D0', POSITIVE_SETTINGS['denominator']]]
    for name, value in POSITIVE_SETTINGS.items():
        global_rows.append(['POSITIVE_SETTINGS.' + name, value])
    files['common_settings.csv'] = csv_text(['parameter', 'value'], global_rows)

    source_hash = hashlib.sha256((ROOT / 'parameters.py').read_bytes()).hexdigest()
    total_scales = sum(len(rows) for rows in SCALES.values())
    total_intervals = sum(len(rows) for rows in INTERVALS.values())
    total_periods = sum(len(values[1]) for values in PARAMETERS.values())
    total_choices = sum(len(row[2]) for rows in INTERVALS.values() for row in rows)

    for pair, values in PARAMETERS.items():
        p, q = pair
        d, B, ell, z, C, m0, tmax, J, target = values
        folder = str(p) + '_' + str(q)
        mode = terminal_mode(pair)
        limit = TERMINAL_LIMITS.get(pair, '')
        parameter_rows = [
            ['p', p], ['q', q], ['d', d], ['P', abs(p)**d], ['s', len(B)],
            ['z', z], ['C', C], ['m0', m0], ['t_max', tmax], ['J', J],
            ['internal_target', target], ['paper_target', PAPER_TARGETS[pair]],
            ['height_limit', HEIGHT_LIMITS[pair]], ['terminal_mode', mode],
            ['terminal_limit', limit], ['chunk', CHUNK],
            ['scale_count', len(SCALES[pair])],
            ['interval_count', len(INTERVALS[pair])],
        ]
        period_rows = [[i, bound, budget]
                       for i, (bound, budget) in enumerate(zip(B, ell), 1)]
        scale_rows = [[i, H, X, n, e]
                      for i, (H, X, n, e) in enumerate(SCALES[pair], 1)]
        interval_rows = []
        interval_md_rows = []
        choice_rows = []
        interval_headers = ['v', 'left', 'right']
        for j in range(J):
            interval_headers += ['i_' + str(j), 'k_' + str(j)]
        for v, (left, right, choices) in enumerate(INTERVALS[pair], 1):
            flat = [v, left, right]
            descriptions = []
            for j, (i, k) in enumerate(choices):
                flat += [i, k]
                descriptions.append('j=' + str(j) + ': (' + str(i) + ', ' + str(k) + ')')
                choice_rows.append([v, j, i, k])
            interval_rows.append(flat)
            interval_md_rows.append([v, left, right, '; '.join(descriptions)])

        files[folder + '/parameters.csv'] = csv_text(['parameter', 'value'], parameter_rows)
        files[folder + '/period.csv'] = csv_text(['i', 'B_i', 'ell_i'], period_rows)
        files[folder + '/scales.csv'] = csv_text(['i', 'H', 'X', 'n', 'e'], scale_rows)
        files[folder + '/intervals.csv'] = csv_text(interval_headers, interval_rows)
        files[folder + '/choices.csv'] = csv_text(['v', 'j', 'i', 'k'], choice_rows)
        special = ('This case internally requires a bound strictly greater than `9001/10000`; the paper states the weaker bound strictly greater than `9/10`.'
                   if pair == (5, 2) else
                   'The internal and paper targets agree for this case. The main verification requires the computed lower bound to exceed this value strictly.')
        if mode == 'signed_window':
            window_note = 'The terminal construction uses the integer window `[-' + str(limit) + ', ' + str(limit) + ']`. Paths leaving the window may be discarded, so the resulting count is a lower bound.'
        elif mode == 'positive_intervals':
            window_note = 'The terminal construction uses the five-step interval recurrence for positive bases. Its `POSITIVE_SETTINGS` are listed in the main supplement; it does not use the `terminal_limit` of the general signed algorithm.'
        else:
            window_note = 'The terminal request is `(h,X)=(1,7)`. The two explicit one-term representations `1` and `p` give the lower bound `2`. No integer window is needed.'
        files[folder + '/README.md'] = (
            '# Complete parameters for bases (' + str(p) + ', ' + str(q) + ')\n\n'
            'This file is generated by `export_parameters.py` from the sole parameter source, `parameters.py`. Every table is complete, with no omitted rows. It records choices to be verified, not computation results or a record of successful verification.\n\n'
            + special + '\n\n' + window_note + '\n\n'
            '## Basic parameters\n\n' + markdown_table(['Parameter', 'Exact value'], parameter_rows)
            + '\n\nSee [parameters.csv](parameters.csv). A blank `terminal_limit` means that this algorithm does not use such a window.\n\n'
            '## Complete period table\n\n' + markdown_table(['i', 'B_i', 'ell_i'], period_rows)
            + '\n\nSee [period.csv](period.csv). Period positions `i` are numbered from 1.\n\n'
            '## All scale choices\n\n' + markdown_table(['Scale i', 'H', 'X', 'n', 'e'], scale_rows)
            + '\n\nSee [scales.csv](scales.csv). The exponent `n` normalizes the scale; `e` only selects a terminal table for positive bases. They are not interchangeable. The value `X` is always the exact, unnormalized radius.\n\n'
            '## All intervals and choices\n\n' + markdown_table(['Interval v', 'Left endpoint', 'Right endpoint', '(Scale i, exponent k), ordered by j'], interval_md_rows)
            + '\n\nSee [intervals.csv](intervals.csv), where each `i_j,k_j` specifies the choice for term `j` on that interval. '
            'The fully expanded [choices.csv](choices.csv) also lists one `(v,j,i,k)` per row. Intervals `v` and scales `i` are numbered from 1; summation terms `j` are numbered from 0.\n'
        )
        summary_rows.append([
            '(' + str(p) + ', ' + str(q) + ')', len(B), len(scale_rows),
            len(interval_rows), J, target, PAPER_TARGETS[pair],
            '[Complete tables](' + folder + '/README.md)',
        ])

    levels = []
    for m in range(1, POSITIVE_SETTINGS['levels'] + 1):
        levels.append([
            m,
            POSITIVE_SETTINGS['seed_exponent'] + POSITIVE_SETTINGS['width_increment'] * (m - 1),
            POSITIVE_SETTINGS['seed_span'] + POSITIVE_SETTINGS['span_increment'] * (m - 1),
            POSITIVE_SETTINGS['height_caps'][m - 1],
        ])
    files['positive_levels.csv'] = csv_text(['m', 'nu_m', 'W_m', 'K_m'], levels)
    files['README.md'] = (
        '# Complete parameter supplement\n\n'
        'The sole parameter source is `parameters.py` in the program root. This directory is generated by `python export_parameters.py`. '
        'The verification program reads `parameters.py`, not the exported tables in this directory or height distributions, terminal counts, or density results from previous runs.\n\n'
        'There are ' + str(len(PARAMETERS)) + ' base pairs, ' + str(total_periods) + ' period positions, '
        + str(total_scales) + ' scales, ' + str(total_intervals) + ' contiguous intervals, and '
        + str(total_choices) + ' interval choices. The case pages and CSV files below list everything, including every endpoint and every `(i,k)`.\n\n'
        '## Case index\n\n'
        + markdown_table(['Bases', 'Period length', 'Scales', 'Intervals', 'J', 'Internal target', 'Paper target', 'Complete data'], summary_rows)
        + '\n\nBoth internal and paper targets are thresholds to be exceeded strictly, not results of a computation. The two targets differ only for `(5,2)`.\n\n'
        '## Common computation settings\n\n'
        + markdown_table(['Parameter', 'Exact value'], global_rows)
        + '\n\nSee [common_settings.csv](common_settings.csv). The setting `retained_heights=(15,28)` means that lifting for positive bases uses only contributions from nonzero residues at heights in this range; other heights are discarded.\n\n'
        '## Five-step settings for positive bases\n\n'
        + markdown_table(['m', 'nu_m', 'W_m', 'K_m'], levels)
        + '\n\nSee [positive_levels.csv](positive_levels.csv). These satisfy `nu_m=13+11*(m-1)` and `W_m=400+2*(m-1)`.\n\n'
        '## Format and consistency checks\n\n'
        'Integers are written in decimal notation, and rational values as reduced fractions `numerator/denominator`. Values are never converted to floating point, and fractions are not replaced by decimal approximations. '
        'CSV files use UTF-8 encoding and comma separators. Import columns containing large integers or fractions as text in editors or spreadsheet software to prevent automatic conversion to approximate decimals.\n\n'
        'Run `python export_parameters.py --check` from the program root to check every exported file against the current parameter source character for character. '
        'This checks parameter transcription only; it does not establish that residue, terminal, or density verification has passed. For the role of each parameter and where it is checked, see the [parameter documentation](../docs/01_parameters.md).\n\n'
        'The SHA-256 of `parameters.py` is\n\n`' + source_hash + '`\n'
    )
    summary = {'pairs': len(PARAMETERS), 'period_entries': total_periods,
               'scales': total_scales, 'intervals': total_intervals,
               'choices': total_choices, 'files': len(files)}
    return files, summary


def export(destination, check=False):
    files, summary = build_files()
    if check:
        failures = []
        for relative, content in files.items():
            path = destination / relative
            if not path.is_file() or path.read_text(encoding='utf-8') != content:
                failures.append(relative)
        if failures:
            raise ValueError('Exported files are missing or differ from the current parameters:\n' + '\n'.join(failures))
    else:
        for relative, content in files.items():
            path = destination / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding='utf-8', newline='\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'supplement',
                        help='Output directory; defaults to supplement in the program directory')
    parser.add_argument('--check', action='store_true', help='Check existing exports without writing files')
    args = parser.parse_args()
    summary = export(args.output, check=args.check)
    print('Parameter export consistency check passed.' if args.check else 'The complete parameter supplement has been generated.')
    for name, value in summary.items():
        print(name + ': ' + str(value))


if __name__ == '__main__':
    main()
