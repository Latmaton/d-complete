"""Exact rational verification: fixed inputs → constants → tails → scales → endpoints.

This module uses fixed parameters and recomputes the required counts. All decisive
comparisons use integers or Fraction. Explicit exceptions keep verification active
even when Python is run with -O.
"""
from fractions import Fraction as F
from math import gcd
from numbers import Integral
from parameters import (PARAMETERS, SCALES, INTERVALS, HEIGHT_LIMITS,
                        PAPER_TARGETS, CHUNK, TERMINAL_LIMITS, POSITIVE_SETTINGS)


def _integer(value, name, minimum=None):
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f'{name} must be an integer; got {value!r}')
    value = int(value)
    if minimum is not None and value < minimum:
        raise ValueError(f'{name} must be at least {minimum}; got {value}')
    return value


def _rational(value, name, minimum=None):
    if isinstance(value, bool) or not isinstance(value, (Integral, F)):
        raise TypeError(f'{name} must be an integer or Fraction; floating-point values are not accepted')
    value = F(int(value)) if isinstance(value, Integral) else value
    if minimum is not None and value < minimum:
        raise ValueError(f'{name} must be at least {minimum}; got {value}')
    return value


def _intervals(rows, scale_count, J, p):
    """Copy and check the full partition without modifying the tables in parameters.py."""
    if not rows:
        raise ValueError('The interval table must not be empty')
    result = []
    end = None
    for v, row in enumerate(rows, 1):
        if len(row) != 3:
            raise ValueError(f'Interval {v} must contain left, right, choices')
        left, right, choices = row
        left = _rational(left, f'left endpoint of interval {v}')
        right = _rational(right, f'right endpoint of interval {v}')
        if not 0 < left < right or (end is not None and left != end):
            raise ValueError(f'Interval {v} must have positive increasing endpoints and meet the preceding interval')
        if len(choices) != J:
            raise ValueError(f'Interval {v} must have exactly J={J} choices')
        checked = []
        for j, choice in enumerate(choices):
            if len(choice) != 2:
                raise ValueError(f'Interval {v}, choice j={j} must be a pair (i,k)')
            i = _integer(choice[0], f'i for interval {v}, j={j}', 1)
            k = _integer(choice[1], f'k for interval {v}, j={j}')
            if i > scale_count:
                raise ValueError(f'Interval {v}, j={j}: scale index i={i} is out of range')
            checked.append((i, k))
        result.append((left, right, checked))
        end = right
    if end != p * result[0][0]:
        raise ValueError('The interval table must cover exactly one full multiplicative period of |p|')
    return result


def read_inputs(p, q):
    """Check all fixed inputs before large searches; supported base pairs may be reversed."""
    p, q = _integer(p, 'p'), _integer(q, 'q')
    if (p, q) not in PARAMETERS and (q, p) in PARAMETERS:
        p, q = q, p
    pair = p, q
    if pair not in PARAMETERS:
        raise ValueError(f'No fixed inputs are available for the base pair {pair}')
    if abs(p) <= 1 or abs(q) <= 1 or gcd(p, q) != 1:
        raise ValueError('The bases must be coprime and have absolute values greater than 1')
    if len(PARAMETERS[pair]) != 9:
        raise ValueError('Each PARAMETERS entry must contain exactly nine items')
    d, B, ell, z, C, m0, tmax, J, target = PARAMETERS[pair]
    d = _integer(d, 'd', 1)
    if not B or len(B) != len(ell):
        raise ValueError('B and ell must be nonempty lists of the same length')
    B = [_integer(b, f'B[{i}]', 1) for i, b in enumerate(B, 1)]
    ell = [_integer(t, f'ell[{i}]', 1) for i, t in enumerate(ell, 1)]
    z, C = _rational(z, 'z'), _rational(C, 'C', 1)
    if z <= 1:
        raise ValueError('z must be strictly greater than 1')
    m0, tmax, J = (_integer(m0, 'm0', 1), _integer(tmax, 't_max', 1),
                   _integer(J, 'J', 1))
    target = _rational(target, 'verification target', 0)
    paper_target = _rational(PAPER_TARGETS[pair], 'paper target', 0)
    if not paper_target <= target < 1:
        raise ValueError('The targets must satisfy 0 <= paper target <= verification target < 1')
    H_MAX = _integer(HEIGHT_LIMITS[pair], 'H_MAX', 1)
    chunk = _integer(CHUNK, 'CHUNK', 1)
    D0 = _integer(POSITIVE_SETTINGS['denominator'], 'D0', 1)
    params = dict(p=p, q=q, d=d, B=B, ell=ell, z=z, C=C,
                  m0=m0, t_max=tmax, J=J, target=target, paper_target=paper_target,
                  D0=D0, P=abs(p)**d, s=len(B), H_MAX=H_MAX, chunk=chunk)
    if pair not in SCALES or not SCALES[pair]:
        raise ValueError(f'No scale table is available for the base pair {pair}')
    scales = []
    for index, source in enumerate(SCALES[pair], 1):
        if len(source) != 4:
            raise ValueError(f'Scale {index} must contain H, X, n, e')
        H, X, n, e = source
        H = _integer(H, f'H of scale {index}', 1)
        X = _rational(X, f'X of scale {index}')
        n = _integer(n, f'n of scale {index}', 0)
        e = _integer(e, f'e of scale {index}', 0)
        if X <= 0:
            raise ValueError(f'X of scale {index} must be positive')
        if pair == (5, 2):
            settings = POSITIVE_SETTINGS
            levels = _integer(settings['levels'], 'levels', 1)
            initial = _integer(settings['seed_exponent'], 'seed_exponent', 1)
            increment = _integer(settings['width_increment'], 'width_increment', 1)
            grid = _integer(settings['grid'], 'grid', 1)
            if e not in range(initial, initial + levels * increment, increment):
                raise ValueError(f'e of scale {index} has no corresponding positive-case interval table')
            column = X * grid / abs(p)**e
            if column.denominator != 1 or not 0 < column <= 2 * grid:
                raise ValueError(f'X of scale {index} is outside the valid grid of its interval table')
        else:
            if e != 0 or X.denominator != 1:
                raise ValueError(f'Signed-case scale {index} must have e=0 and an integer X')
            if pair in TERMINAL_LIMITS:
                limit = _integer(TERMINAL_LIMITS[pair], 'terminal window', 1)
                if X > limit:
                    raise ValueError(f'X of scale {index} exceeds the terminal window')
            elif pair in ((7, -2), (-7, -2)):
                if H != 1 or X != 7:
                    raise ValueError('The base-seven cases provide only the two single-term representations for H=1, X=7')
            else:
                raise ValueError(f'No terminal counting method is available for the base pair {pair}')
        scales.append(dict(H=H, X=X, n=n, count_exponent=e,
                           tail_terms=min(H, tmax)))
    if pair not in INTERVALS:
        raise ValueError(f'No interval table is available for the base pair {pair}')
    intervals = _intervals(INTERVALS[pair], len(scales), J, abs(p))
    constants(params)  # Check contraction and multiplicative order before large searches.
    return params, scales, intervals


def constants(params):
    """Use only fixed inputs; run before residue searches and terminal counting."""
    p, q = _integer(params['p'], 'p'), _integer(params['q'], 'q')
    if abs(p) <= 1 or abs(q) <= 1 or gcd(p, q) != 1:
        raise ValueError('The constants require coprime bases with absolute values greater than 1')
    P = _integer(params['P'], 'P', 2)
    B = [_integer(b, 'B_i', 1) for b in params['B']]
    ell = [_integer(t, 'ell_i', 1) for t in params['ell']]
    if not B or len(B) != len(ell):
        raise ValueError('B and ell must be nonempty lists of the same length')
    s = len(B)
    c = F(abs(q)**sum(ell), P**s)
    if c >= 1:
        raise ValueError(f'The contraction condition c < 1 fails: c={c}')
    A = sum((F(B[i]*abs(q)**(sum(ell[i:])-1), P**(s-i))
             for i in range(s)), F(0)) / (1-c)
    kappa = 1 if p > 0 and q > 0 else 2
    for order in range(1, abs(q) + 1):
        if pow(p, order, abs(q)) == 1:
            break
    else:
        raise ValueError('Could not determine the multiplicative order of p modulo |q|')
    theta = F(kappa*order, abs(q))
    return c, A, kappa, theta


def tails(params, hist):
    """Start with the exponential-moment bound, then update full periods with upward rounding."""
    P, B, ell = params['P'], params['B'], params['ell']
    z, C, D0 = params['z'], params['C'], params['D0']
    s = len(B)
    moments = []
    checked_hist = {}
    for i in range(s):
        if B[i] not in hist:
            raise ValueError(f'Missing height histogram for B={B[i]}')
        row = {_integer(h, 'histogram height', 0): _integer(number, 'histogram count', 0)
               for h, number in hist[B[i]].items()}
        if sum(row.values()) != P or row.get(0) != 1:
            raise ValueError(f'The histogram for B={B[i]} must cover P residues, with exactly one at height 0')
        checked_hist[B[i]] = row
        moments.append(sum((number*z**(h-ell[i])
                            for h, number in row.items()), F(0))/P)
    g = F(1)
    for value in moments:
        g *= value
    partial = F(1)
    Lambda = F(1)
    for i in range(s-1, 0, -1):
        partial *= moments[i]
        Lambda = max(Lambda, partial)
    if g >= 1 or Lambda > C:
        raise ValueError(f'Exponential-moment conditions fail: g < 1 and Lambda_s <= C are required; g={g}, Lambda_s={Lambda}')

    # The largest accessed index is t+ell[i]; retain the extra entries needed for all periods.
    length = params['t_max'] + params['m0']*sum(ell)
    w = [D0]
    for t in range(1, length+1):
        value = D0*min(F(1), C/z**t)
        w.append(-(-value.numerator // value.denominator))
    for step in range(params['m0']*s):
        i = step % s
        length -= ell[i]
        new = [D0]
        for t in range(1, length+1):
            total = 0
            for h, number in checked_hist[B[i]].items():
                u = t+ell[i]-h
                total += number*(D0 if u <= 0 else w[u])
            new.append(-(-total // P))
        w = new
    w = [F(number, D0) for number in w]
    if len(w) != params['t_max']+1 or w[0] != 1:
        raise RuntimeError('Incorrect length or initial value in the tail recurrence')
    if any(not 0 <= w[t+1] <= w[t] <= 1 for t in range(len(w)-1)):
        raise RuntimeError('Tail bounds must be nonincreasing and lie in [0,1]')
    return w, g, Lambda


def scale_coefficients(params, scales, w, terminal, A):
    """Derive valid rational scales from terminal and tail bounds; fill a,b in place."""
    p, q, D0 = abs(params['p']), abs(params['q']), params['D0']
    A = _rational(A, 'A', 0)
    w = [_rational(value, f'w[{t}]', 0) for t, value in enumerate(w)]
    if not w or w[0] != 1 or any(w[t+1] > w[t] for t in range(len(w)-1)):
        raise ValueError('Tail bounds must start at 1 and be nonincreasing')
    for index, row in enumerate(scales, 1):
        H, X, n = row['H'], row['X'], row['n']
        terms = _integer(row['tail_terms'], 'tail_terms', 1)
        if terms >= len(w) or terms > H:
            raise ValueError(f'Scale {index} has insufficient tail entries or an invalid height')
        boundary = (A*q**H+X)/p**n
        scaled = D0*boundary
        # Use the strictly higher grid point; add 1 even when boundary is already on the grid.
        row['a'] = F(scaled.numerator//scaled.denominator+1, D0)
        beta = F(0)
        for v in range(terms):
            weight = w[v]-w[v+1]
            if weight > 0:
                key = H-v, X
                if key not in terminal:
                    raise ValueError(f'Scale {index} is missing the terminal count {key}')
                lower = _rational(terminal[key], f'terminal count {key}', 0)
                beta += weight*lower
        scaled = D0**2*beta/p**n
        row['b'] = F(scaled.numerator//scaled.denominator, D0**2)
        if row['a'] <= boundary or row['b'] <= 0:
            raise ValueError(f'Scale {index} requires a strictly above the boundary and a positive lower bound b')


def density(params, scales, intervals, kappa, theta, records=None):
    """Check both endpoints of each interval exactly; optionally record every contribution."""
    if records is not None and not isinstance(records, list):
        raise TypeError('records must be a list or None')
    p, q, J = abs(params['p']), abs(params['q']), params['J']
    kappa = _integer(kappa, 'kappa', 1)
    theta = _rational(theta, 'theta', 0)
    intervals = _intervals(intervals, len(scales), J, p)
    minimum = None
    for v, (left, right, choices) in enumerate(intervals, 1):
        for y in (left, right):
            total = F(0)
            contributions = []
            for j, (i, k) in enumerate(choices):
                a = _rational(scales[i-1]['a'], f'a of scale {i}', 0)
                b = _rational(scales[i-1]['b'], f'b of scale {i}', 0)
                factor = F(p)**k
                penalty = theta*max(F(0), factor*a-y/q**j)
                contribution = factor*b-penalty
                total += contribution
                if records is not None:
                    contributions.append(dict(j=j, i=i, k=k, factor=factor,
                                              a=a, b=b, penalty=penalty,
                                              contribution=contribution))
            value = total/(kappa*(1-F(1, q**J))*y)
            if records is not None:
                records.append(dict(interval_index=v, y=y, left=left, right=right,
                                    value=value, numerator=total,
                                    contributions=contributions))
            minimum = value if minimum is None else min(minimum, value)
    if minimum <= params['target']:
        raise ValueError(f'The minimum endpoint value {minimum} does not strictly exceed the verification target {params["target"]}')
    return minimum
