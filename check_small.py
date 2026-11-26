"""Independent checks with small parameters. All checks also run under python -O.

Expected values come from direct enumeration of exponent chains, sets, or interval subsets, independently of the earlier program.
These small tests do not replace running main.py with the full parameters in the paper.
"""
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import subprocess
import sys
import argparse
import numpy as np
from residues import (signed_families, boolean_families, odd_heights,
                      positive_families, verify_family)
from terminals import (signed_terminal_counts, seed_intervals,
                       next_intervals, pack_intervals)


def enumerate_families(p, q, d, H, bounds):
    """Enumerate every admissible exponent chain to find the minimum height of each residue."""
    P = abs(p)**d
    heights = {B: {0: 0} for B in bounds}
    for length in range(1, min(d, H) + 1):
        for aa in combinations(range(d), length):
            for gg in combinations(range(H), length):
                value = sum((Fraction(p**a, q**g) for a, g in zip(aa, gg)), Fraction(0))
                r = value.numerator * pow(value.denominator, -1, P) % P
                for B in bounds:
                    if abs(value) <= B:
                        heights[B][r] = min(heights[B].get(r, H+1), gg[-1]+1)
    hist = {}
    for B, row in heights.items():
        if not (len(row) == P):
            raise AssertionError("Comparison with direct enumeration failed")
        hist[B] = {}
        for h in row.values():
            hist[B][h] = hist[B].get(h, 0) + 1
    return hist


def check_residues():
    for p, q in [(-2, 3), (-2, -3), (-5, 2), (5, -2), (-5, -2)]:
        params = {'p': p, 'q': q, 'd': 3, 'H_MAX': 16, 'B': [2, 3, 5]}
        expected = enumerate_families(p, q, 3, 16, [2, 3, 5])
        if not (signed_families(params, False) == expected):
            raise AssertionError("Comparison with direct enumeration failed")
    for p in [-7, 7]:
        params = {'p': p, 'q': -2, 'd': 2, 'H_MAX': 18, 'B': [7]}
        expected = enumerate_families(p, -2, 2, 18, [7])
        if not (boolean_families(params, False) == expected):
            raise AssertionError("Comparison with direct enumeration failed")
    # Enumerate integer antichains directly and compare with the minimum-height recurrence for odd integers.
    limit, width = 199, 3
    expected = {n: 255 for n in range(1, limit+1, 2)}
    for length in range(1, width+1):
        for aa in combinations(range(width), length):
            for bb in combinations(range(limit.bit_length()), length):
                n = sum(5**a * 2**b for a, b in zip(aa, reversed(bb)))
                if n in expected:
                    expected[n] = min(expected[n], max(bb)+1)
    actual = odd_heights(limit, width, verbose=False)
    if not (all(int(actual[(n-1)//2]) == h for n, h in expected.items())):
        raise AssertionError("Comparison with direct enumeration failed")
    print('Small tests passed for residue search, the Boolean recurrence, and the height recurrence for positive integers.')


def check_terminals():
    # Insert representations into a set one at a time, independently of the Boolean-array index formulas.
    for p, q in [(-2, 3), (-2, -3), (-5, 2), (5, -2), (-5, -2)]:
        requests = {(h, X) for h in range(1, 7) for X in (1, 7, 20, 50, 100)}
        got = signed_terminal_counts(p, q, requests, limit=1000, verbose=False)
        previous = {0}
        for h in range(1, 7):
            current = previous | {p*n + q**(h-1) for n in previous
                                  if abs(p*n + q**(h-1)) <= 1000}
            while True:
                more = {p*n for n in current if abs(p*n) <= 1000}
                if more <= current:
                    break
                current |= more
            for X in (1, 7, 20, 50, 100):
                expected = sum(abs(n) <= X and n % q != 0 for n in current)
                if not (got[h, X] == expected):
                    raise AssertionError("Comparison with direct enumeration failed")
            previous = current
    print('Signed terminal counts: 150 comparisons with set enumeration passed.')

    # Enumerate small antichains individually, without using the height recurrence to obtain expected values.
    possible = {}
    for size in range(1, 4):
        for aa in combinations(range(3), size):
            for bb in combinations(range(6), size):
                n = sum(5**a * 2**b for a, b in zip(aa, reversed(bb)))
                if n <= 50 and n % 2:
                    possible[n] = min(possible.get(n, 255), max(bb) + 1)
    got = seed_intervals(grid=25, denominator=1000, chunk=7,
                         seed_exponent=2, seed_width=3, max_height=10,
                         span=7, verbose=False)
    for h in range(11):
        for i in range(51):
            count = sum(hh <= h and i - 7 < n <= i for n, hh in possible.items())
            if not (int(got[h, i]) == count * 1000 // 25):
                raise AssertionError("Comparison with direct enumeration failed")
    print('Positive seeds: 561 comparisons with antichain enumeration passed.')

    # Enumerate all subsets of intervals and check disjointness directly.
    table = np.array([[0, 2, 0, 7, 1, 4, 3, 8]], dtype=np.uint32)
    for span in (1, 2, 3, 5):
        got = pack_intervals(table, span)
        for last in range(8):
            best = 0
            for size in range(last + 2):
                for chosen in combinations(range(last + 1), size):
                    if all(chosen[j+1] - chosen[j] >= span for j in range(len(chosen)-1)):
                        best = max(best, sum(int(table[0, i]) for i in chosen))
            if not (int(got[0, last]) == best):
                raise AssertionError("Comparison with direct enumeration failed")
    print('Interval packing: 32 comparisons with exhaustive subset enumeration passed.')


def check_positive_lift():
    """Evaluate the double sum in the appendix term by term and compare with grouped shifts and additions.

    The small arrays here are synthetic inputs used to check the formula implementation, not mathematical certificates.
    The checks include nonzero shifts, the zero residue, truncation of earlier heights, and maximization over two size bounds.
    """
    Q, R, D, P = 8, 16, 1000, 5000
    previous = np.zeros((41, 2*Q+1), dtype=np.uint32)
    for h in range(1, 41):
        for i in range(1, 2*Q+1):
            previous[h, i] = min(D, 7*h+13*i)
    joint = {}
    for B in (3, 4):
        counts = np.zeros((32, 65), dtype=np.uint32)
        counts[0, 0] = 1
        for a in range(15, 29):
            for b in range(65):
                counts[a, b] = (a+3*b+B) % 5
        joint[B] = counts
    actual = next_intervals(previous, joint, 2, grid=Q, size_grid=R,
                            denominator=D, P=P, height_caps=(40, 57),
                            verbose=False)
    expected = np.zeros(actual.shape, dtype=np.uint32)
    for h in range(1, 58):
        numerator = Q*2**(h-1)
        divisor = R*5**24
        for i in range(2, 2*Q+1):
            best = 0
            for counts in joint.values():
                total = int(previous[min(h, 40), i-2])
                for a in range(15, min(h, 28)+1):
                    for b in range(65):
                        shift = numerator*b // divisor
                        old_column = i-shift-2
                        if old_column >= 0:
                            total += int(counts[a, b])*int(
                                previous[min(h-a, 40), old_column])
                best = max(best, total // P)
            expected[h, i] = max(best, expected[h-1, i])
    if not (np.array_equal(actual, expected)):
        raise AssertionError("Comparison with direct enumeration failed")
    print('Positive interval lifting: 986 comparisons with direct double summation passed.')


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def expect_error(function, message, errors=(ValueError, RuntimeError)):
    try:
        function()
    except errors:
        return
    raise AssertionError('An invalid input was not rejected: ' + message)


def direct_odd_heights(limit, width):
    """The positive sum is at most limit, so every exponent of 2 is less than limit.bit_length()."""
    result = np.full((limit + 1) // 2, 255, dtype=np.uint8)
    for size in range(1, min(width, limit.bit_length()) + 1):
        for aa in combinations(range(width), size):
            for bb in combinations(range(limit.bit_length()), size):
                n = sum(5**a * 2**b for a, b in zip(aa, reversed(bb)))
                if n <= limit and n % 2:
                    i = (n - 1) // 2
                    result[i] = min(int(result[i]), bb[-1] + 1)
    return result


def check_height_chunks():
    for limit in (1, 2, 199, 200):
        for width in (1, 2, 3, 4):
            expected = direct_odd_heights(limit, width)
            for chunk in (1, 2, 7, 1000):
                actual = odd_heights(limit, width, chunk=chunk, verbose=False)
                require(np.array_equal(actual, expected),
                        'Minimum height for odd integers differs: ' + str((limit, width, chunk)))
    # Also check chunk boundaries in the Boolean recurrence and rational size bounds in the signed search.
    for p, q, d in [(-7, -2, 2), (7, -2, 2), (-2, 3, 3)]:
        B = abs(p)**(d - 1)
        params = dict(p=p, q=q, d=d, H_MAX=18, B=[B, B+1])
        expected = enumerate_families(p, q, d, 18, params['B'])
        for chunk in (1, 7, 1000):
            require(boolean_families(params, False, chunk=chunk) == expected,
                    'The Boolean-family results differ between chunks')
    bounds = [Fraction(3, 2), 2]
    params = dict(p=-5, q=2, d=3, H_MAX=16, B=bounds)
    require(signed_families(params, False) ==
            enumerate_families(-5, 2, 3, 16, bounds), 'The rational size-bound check failed')
    print('Two-array heights: 64 endpoint/width/chunk checks passed, as did the Boolean chunk and rational size-bound checks.')


def check_positive_families():
    """Enumerate rational chains directly and choose the representative minimizing (h,sigma) for each residue."""
    cap, d, height, size_grid = 999, 2, 16, 7
    bounds = (1, 2, 5)
    chosen = {B: {0: (0, Fraction(0))} for B in bounds}
    for size in range(1, d + 1):
        for aa in combinations(range(d), size):
            for gg in combinations(range(height), size):
                value = sum((Fraction(5**a, 2**g)
                             for a, g in zip(aa, gg)), Fraction(0))
                h = gg[-1] + 1
                numerator = value * 2**(h - 1)
                require(numerator.denominator == 1, 'The numerator after reflection must be an integer')
                if numerator > cap:
                    continue
                residue = value.numerator * pow(value.denominator, -1, 5**d) % (5**d)
                for B, row in chosen.items():
                    if value <= B:
                        candidate = h, value
                        if residue not in row or candidate < row[residue]:
                            row[residue] = candidate
    expected_hist, expected_joint = {}, {}
    for B, row in chosen.items():
        require(len(row) == 5**d, 'The positive test fixture does not cover every residue')
        expected_hist[B] = {}
        counts = np.zeros((height + 1, B * size_grid + 1), dtype=np.uint32)
        for h, value in row.values():
            scaled = value * size_grid
            b = scaled.numerator // scaled.denominator
            counts[h, b] += 1
            expected_hist[B][h] = expected_hist[B].get(h, 0) + 1
        expected_joint[B] = counts
    for chunk in (1, 7, 1000):
        hist, joint = positive_families(cap=cap, d=d, bounds=bounds,
                                       hmax=height, size_grid=size_grid,
                                       chunk=chunk, verbose=False)
        require(hist == expected_hist, 'The positive height distribution differs')
        for B in bounds:
            require(np.array_equal(joint[B], expected_joint[B]),
                    'The positive joint distribution differs: ' + str((B, chunk)))
    print('Positive residue families: joint distributions for 3 size bounds and 3 chunk sizes agree with direct rational-chain enumeration.')


def check_seed_boundaries():
    base, width, height = 25, 3, 10
    expected = direct_odd_heights(2 * base, width)
    checks = 0
    for grid, span, chunk in [(17, 7, 1), (31, 1, 1000)]:
        got = seed_intervals(grid=grid, denominator=1000, chunk=chunk,
                             seed_exponent=2, seed_width=width,
                             max_height=height, span=span, verbose=False)
        for h in range(height + 1):
            for i in range(2 * grid + 1):
                # Compare the rational endpoints directly, without copying the rounded-index formula.
                left = Fraction((i - span) * base, grid)
                right = Fraction(i * base, grid)
                count = sum(int(hh) <= h and left < 2*j + 1 <= right
                            for j, hh in enumerate(expected))
                require(int(got[h, i]) == count * 1000 // base,
                        'The seed interval boundary check failed: ' + str((grid, span, h, i)))
                checks += 1
    # A smaller height cap should omit representations above that cap, never count them prematurely.
    small = seed_intervals(grid=17, denominator=1000, chunk=7,
                           seed_exponent=2, seed_width=width,
                           max_height=2, span=7, verbose=False)
    large = seed_intervals(grid=17, denominator=1000, chunk=7,
                           seed_exponent=2, seed_width=width,
                           max_height=height, span=7, verbose=False)
    require(np.array_equal(small, large[:3]), 'Seed height truncation differs')
    print('Positive seeds: a further ' + str(checks) + ' nonintegral-endpoint checks and the height-truncation check passed.')


def check_tail_bound():
    """Enumerate finite height sequences individually and check upper bounds on the true tail probabilities at the ends of complete periods."""
    from arithmetic import tails
    params = dict(P=3, B=[1, 1], ell=[1, 3], z=Fraction(2),
                  C=Fraction(1), D0=1000, t_max=5, m0=2)
    w, g, Lambda = tails(params, {1: {0: 1, 1: 1, 3: 1}})
    require(g == Fraction(121, 144) and Lambda == 1, 'The exponential moments for the small example differ')
    for length in (4, 6):
        exceed = [0] * 6
        for heights in product((0, 1, 3), repeat=length):
            excess = 0
            for i, h in enumerate(heights):
                excess = max(0, excess + h - params['ell'][i % 2])
            for t in range(1, 6):
                exceed[t] += excess >= t
        for t in range(1, 6):
            require(Fraction(exceed[t], 3**length) <= w[t],
                    'A true tail probability exceeds the computed upper bound')
    print('Tail bounds: all 810 height sequences of lengths 4 and 6 passed the checks against 5 tail-probability bounds.')


def check_errors():
    """These checks also run in a separate python -O subprocess."""
    from arithmetic import constants, read_inputs, tails
    require(verify_family(-2, 3, 3, 0, (), 0) == 0, 'The empty witness must have height 0')
    require(verify_family(-2, 3, 3, 1, ((0, 0),), 1) == 1,
            'The one-term witness must have height 1')
    require(signed_terminal_counts(-2, 3, set(), verbose=False) == {},
            'Empty terminal requests must return an empty dictionary')
    requests = {(0, 0), (0, 5), (1, 0)}
    require(signed_terminal_counts(-2, 3, requests, limit=5, verbose=False) ==
            {key: 0 for key in requests}, 'Counts at height 0 and at X=0 must be 0')
    cases = [
        ('incorrect residue', lambda: verify_family(-2, 3, 3, 2, ((0, 0),), 1)),
        ('size-bound failure', lambda: verify_family(-2, 3, 3, 1, ((0, 0),), 0)),
        ('not an antichain', lambda: verify_family(-2, 3, 3, 1, ((0, 1), (1, 0)), 10)),
        ('noncoprime bases', lambda: signed_families(dict(p=6, q=3, d=2, H_MAX=3, B=[2]), False)),
        ('empty size bounds', lambda: signed_families(dict(p=-2, q=3, d=2, H_MAX=3, B=[]), False)),
        ('incomplete coverage by the signed family', lambda: signed_families(dict(p=-5, q=2, d=2, H_MAX=1, B=[1]), False)),
        ('the Boolean family cannot omit the size bound', lambda: boolean_families(dict(p=7, q=-2, d=2, H_MAX=2, B=[1]), False)),
        ('incomplete coverage by the Boolean family', lambda: boolean_families(dict(p=7, q=-2, d=2, H_MAX=1, B=[7]), False)),
        ('empty positive size bounds', lambda: positive_families(cap=9, d=1, bounds=(), hmax=3, verbose=False)),
        ('incomplete positive coverage', lambda: positive_families(cap=1, d=2, bounds=(1,), hmax=3, verbose=False)),
        ('empty height range', lambda: odd_heights(0, 1, verbose=False)),
        ('zero width', lambda: odd_heights(9, 0, verbose=False)),
        ('zero chunk size', lambda: odd_heights(9, 1, chunk=0, verbose=False)),
        ('negative terminal height', lambda: signed_terminal_counts(-2, 3, {(-1, 1)}, verbose=False)),
        ('noninteger terminal boundary', lambda: signed_terminal_counts(-2, 3, {(1, Fraction(1, 2))}, verbose=False)),
        ('terminal window too small', lambda: signed_terminal_counts(-2, 3, {(1, 5)}, limit=4, verbose=False)),
        ('noncoprime terminal bases', lambda: signed_terminal_counts(6, 3, {(1, 5)}, verbose=False)),
        ('zero seed grid size', lambda: seed_intervals(grid=0, seed_exponent=1, seed_width=1, verbose=False)),
        ('zero interval width', lambda: pack_intervals(np.zeros((1, 3), dtype=np.uint32), 0)),
        ('empty interval table', lambda: pack_intervals(np.zeros((0, 0), dtype=np.uint32), 1)),
        ('unknown fixed parameters', lambda: read_inputs(4, 6)),
    ]
    previous = np.zeros((41, 17), dtype=np.uint32)
    counts = np.zeros((32, 65), dtype=np.uint32)
    counts[0, 0] = 1
    options = dict(grid=8, size_grid=16, denominator=1000, P=5000,
                   height_caps=(40, 57), verbose=False)
    cases.extend([
        ('empty joint distribution', lambda: next_intervals(previous, {}, 2, **options)),
        ('incorrect interval-table shape', lambda: next_intervals(previous[:, :-1], {4: counts}, 2, **options)),
        ('interval count exceeds the denominator', lambda: next_intervals(previous + 1001, {4: counts}, 2, **options)),
    ])
    params, _, _ = read_inputs(-5, 2)
    bad_constants = dict(params, P=2)
    bad_hist = {B: {0: 1} for B in params['B']}
    cases.extend([
        ('contraction-condition failure', lambda: constants(bad_constants)),
        ('incorrect height-distribution total', lambda: tails(params, bad_hist)),
    ])
    for label, function in cases:
        expect_error(function, label)
    print('Input and certificate validation: ' + str(len(cases)) + ' exception cases and the empty/one-term checks passed.')
    return len(cases)


def check_optimized_mode():
    result = subprocess.run([sys.executable, '-O', str(Path(__file__).resolve()),
                             '--errors-only'], capture_output=True, text=True)
    require(result.returncode == 0,
            'Validation under python -O failed: \n' + result.stdout + result.stderr)
    print('python -O: ' + result.stdout.strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--errors-only', action='store_true',
                        help='check only input and certificate validation (used by the -O subprocess)')
    args = parser.parse_args()
    if args.errors_only:
        check_errors()
        return
    check_residues()
    check_terminals()
    check_positive_lift()
    check_height_chunks()
    check_positive_families()
    check_seed_boundaries()
    check_tail_bound()
    check_errors()
    check_optimized_mode()
    print('All small-parameter tests passed.')


if __name__ == '__main__':
    main()
