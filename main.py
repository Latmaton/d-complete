"""Construct and verify eight fixed cases from scratch. Example: python main.py -2 -3."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import time

from parameters import (PARAMETERS, POSITIVE_SETTINGS, TERMINAL_LIMITS,
                        CHUNK, PAPER_TARGETS)
from arithmetic import read_inputs, constants, tails, scale_coefficients, density
from residues import residue_counts
from terminals import terminal_counts
from reporting import (write_json, source_hashes, environment, residue_report,
                       tail_report, terminal_report, scale_report, endpoint_report)


def resource_notes(pair):
    """Estimate the main arrays only, not peak process memory."""
    p, q = pair
    if pair == (5, 2):
        settings = POSITIVE_SETTINGS
        base = 5 ** settings['seed_exponent']
        return {'odd_seed_count': base,
                'two_seed_height_arrays_bytes': 2 * base,
                'residue_modulus': 5 ** PARAMETERS[pair][0],
                'note': 'This large positive case also requires joint counts and temporary arrays, excluded here.'}
    if pair in TERMINAL_LIMITS:
        limit = TERMINAL_LIMITS[pair]
        length = (2 * limit) // 3 + 1 if pair == (-2, 3) else 2 * limit + 1
        return {'terminal_window_radius': limit,
                'two_terminal_boolean_arrays_bytes_at_most': 2 * length,
                'note': 'Residue searches and temporary arrays require additional memory.'}
    d = PARAMETERS[pair][0]
    return {'two_residue_boolean_tables_bytes':
            2 * sum(abs(p) ** e for e in range(d + 1)),
            'note': 'The terminal lower bound follows from two single-term representations.'}


def parameter_check(pair, chunk=CHUNK):
    params, scales, intervals = read_inputs(*pair)
    if not isinstance(chunk, int) or isinstance(chunk, bool) or chunk < 1:
        raise ValueError('chunk must be a positive integer')
    params['chunk'] = chunk
    # Run these inexpensive checks before allocating large arrays.
    const = constants(params)
    return params, scales, intervals, const


def run_case(pair, output='results', chunk=CHUNK, verbose=True, save_arrays=False):
    params, scales, intervals, const = parameter_check(pair, chunk)
    pair = params['p'], params['q']
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    folder = Path(output) / (f'{pair[0]}_{pair[1]}_' + stamp)
    folder.mkdir(parents=True, exist_ok=False)
    hashes = source_hashes()
    summary = {'status': 'running', 'pair': pair,
               'started_utc': stamp, 'current_stage': 'residues', 'source_sha256': hashes,
               'environment': environment(), 'paper_target': params['paper_target'],
               'verification_target': params['target'], 'stages_seconds': {},
               'resources': resource_notes(pair)}
    write_json(folder / 'summary.json', summary)
    write_json(folder / 'parameters.json', {'parameters': params,
               'scales': scales, 'intervals': intervals,
               'positive_settings': POSITIVE_SETTINGS if pair == (5, 2) else None})
    c, A, kappa, theta = const
    write_json(folder / 'constants.json', {'c': c, 'A': A, 'kappa': kappa, 'theta': theta})
    start = time.perf_counter()
    stage = 'residues'
    try:
        print(f'\nBases {pair}, output directory {folder}', flush=True)
        print('1/4 Constructing residue families from scratch', flush=True)
        now = time.perf_counter()
        hist, joint = residue_counts(params, verbose=verbose)
        summary['stages_seconds'][stage] = time.perf_counter() - now
        residue_report(folder, hist, joint, save_arrays)
        write_json(folder / 'summary.json', summary)

        stage = 'tails'
        summary['current_stage'] = stage
        write_json(folder / 'summary.json', summary)
        print('2/4 Checking exponential moments and computing tail bounds', flush=True)
        now = time.perf_counter()
        w, g, Lambda = tails(params, hist)
        summary['stages_seconds'][stage] = time.perf_counter() - now
        tail_report(folder, w, g, Lambda)
        write_json(folder / 'summary.json', summary)

        stage = 'terminals'
        summary['current_stage'] = stage
        write_json(folder / 'summary.json', summary)
        print('3/4 Constructing terminal lower bounds from scratch', flush=True)
        now = time.perf_counter()
        terminal = terminal_counts(params, scales, joint, verbose=verbose)
        summary['stages_seconds'][stage] = time.perf_counter() - now
        terminal_report(folder, terminal)
        # The joint tables are no longer needed, so release their memory.
        del joint
        write_json(folder / 'summary.json', summary)

        stage = 'density'
        summary['current_stage'] = stage
        write_json(folder / 'summary.json', summary)
        print('4/4 Computing scale coefficients and checking every endpoint', flush=True)
        now = time.perf_counter()
        scale_coefficients(params, scales, w, terminal, A)
        records = []
        lower = density(params, scales, intervals, kappa, theta, records=records)
        summary['stages_seconds'][stage] = time.perf_counter() - now
        scale_report(folder, scales)
        endpoint_report(folder, records)
        if source_hashes() != hashes:
            raise RuntimeError('Core source files changed during the run; rerun with a fixed version')
        summary.update(status='passed', current_stage='complete', lower_bound=lower,
                       terminal_request_count=len(terminal),
                       scale_count=len(scales), interval_count=len(intervals),
                       endpoint_count=len(records),
                       elapsed_seconds=time.perf_counter() - start)
        write_json(folder / 'summary.json', summary)
        print('All checks passed. Lower density >=', lower, '>', params['target'], flush=True)
        print('Decimal for reference only:', float(lower), flush=True)
        return summary, folder
    except BaseException as error:
        summary.update(status='interrupted' if isinstance(error, KeyboardInterrupt) else 'failed',
                       failed_stage=stage, error=type(error).__name__ + ': ' + str(error),
                       elapsed_seconds=time.perf_counter() - start)
        write_json(folder / 'summary.json', summary)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('p', type=int, nargs='?')
    parser.add_argument('q', type=int, nargs='?')
    parser.add_argument('--all', action='store_true', help='verify all eight cases sequentially from scratch')
    parser.add_argument('--list', action='store_true', help='list supported cases')
    parser.add_argument('--check-parameters', action='store_true', help='check only the parameters and inexpensive constant conditions')
    parser.add_argument('--plan', action='store_true', help='show the cases and main array sizes without starting a full run')
    parser.add_argument('--output', default='results', help='parent directory for new run results')
    parser.add_argument('--chunk', type=int, default=CHUNK, help='positions per chunk; affects temporary memory only')
    parser.add_argument('--quiet', action='store_true', help='show only the main stages')
    parser.add_argument('--save-arrays', action='store_true', help='also save compressed joint counts for the positive case')
    args = parser.parse_args(argv)
    if args.chunk < 1:
        parser.error('--chunk must be a positive integer')
    if (args.p is None) != (args.q is None):
        parser.error('supply both p and q')
    if args.all and args.p is not None:
        parser.error('--all cannot be combined with a single pair of bases')
    if args.list:
        for pair in PARAMETERS:
            print(pair, 'paper target', PAPER_TARGETS[pair], 'verification target', PARAMETERS[pair][-1])
        return 0
    if args.all or ((args.check_parameters or args.plan) and args.p is None):
        pairs = list(PARAMETERS)
    elif args.p is not None:
        pair = args.p, args.q
        if pair not in PARAMETERS and pair[::-1] in PARAMETERS:
            pair = pair[::-1]
        if pair not in PARAMETERS:
            parser.error('no fixed parameters for this pair; run --list')
        pairs = [pair]
    else:
        parser.print_help()
        return 0
    if args.check_parameters or args.plan:
        for pair in pairs:
            parameter_check(pair, args.chunk)
            print(pair, 'parameter and constant checks passed')
            if args.plan:
                print(resource_notes(pair))
        print('This mode does not verify residue families, terminal counts, or the final density bound.')
        return 0
    try:
        for pair in pairs:
            run_case(pair, args.output, args.chunk, not args.quiet, args.save_arrays)
    except KeyboardInterrupt:
        print('\nCalculation interrupted; completed stage records have been saved.')
        return 130
    except (ValueError, RuntimeError, OverflowError, MemoryError) as error:
        print('Verification did not finish:', error)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
