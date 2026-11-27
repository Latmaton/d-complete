"""Construct terminal counting lower bounds from the proved recurrences.

The positive case is explained in docs/03_positive_terminals.md, and the signed
cases in docs/03_signed_terminals.md. Integer checks do not rely on assert.
"""
from fractions import Fraction as F
from math import gcd
import numpy as np
from parameters import POSITIVE_SETTINGS, TERMINAL_LIMITS, CHUNK
from residues import odd_heights


INT64_MAX = 2**63 - 1
UINT32_MAX = 2**32 - 1
UINT64_MAX = 2**64 - 1


def integer(value, name, minimum=None):
    """Accept integral values and exact rationals with denominator 1, rejecting nonintegral values and booleans."""
    try:
        result = int(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(name + ' must be an integer') from error
    if isinstance(value, (bool, np.bool_)) or result != value:
        raise ValueError(name + ' must be an integer')
    if minimum is not None and result < minimum:
        raise ValueError(name + ' must be at least ' + str(minimum))
    return result


def integer_table(table, name):
    if not isinstance(table, np.ndarray) or table.ndim != 2 or 0 in table.shape:
        raise ValueError(name + ' must be a nonempty two-dimensional array')
    if not np.issubdtype(table.dtype, np.integer):
        raise ValueError(name + ' must contain integers')
    if np.issubdtype(table.dtype, np.signedinteger) and np.any(table < 0):
        raise ValueError(name + ' must not contain negative values')


def ceil_div(a, b):
    if b <= 0:
        raise ValueError('The divisor in ceiling division must be positive')
    return -((-a) // b)


def add_image(source, target, p, shift, first, step, lo, hi, chunk=CHUNK):
    """Mark images p*n+shift of sources n in the window; position i represents first+step*i."""
    chunk = integer(chunk, 'chunk', 1)
    if source.ndim != 1 or target.shape != source.shape:
        raise ValueError('Source and target must be one-dimensional arrays of equal length')
    if source.dtype != np.bool_ or target.dtype != np.bool_:
        raise ValueError('Source and target must be Boolean arrays')
    if step <= 0:
        raise ValueError('The common difference of the represented arithmetic progression must be positive')
    if len(source) == 0 or lo > hi:
        return
    last = first + step * (len(source) - 1)
    lo, hi = max(lo, first), min(hi, last)
    if lo > hi:
        return
    image_lo, image_hi = sorted((p * lo + shift, p * hi + shift))
    # Use Python integers to exclude images outside the window before passing large q powers to NumPy.
    if image_hi < first or image_lo > last:
        return
    bound = max(abs(first), abs(last))
    if abs(p) * bound + abs(shift) > INT64_MAX or 2 * bound > INT64_MAX:
        raise ValueError('The affine image calculation exceeds the int64 range')
    left = max(0, ceil_div(lo - first, step))
    right = min(len(source), (hi - first) // step + 1)
    for start in range(left, right, chunk):
        stop = min(start + chunk, right)
        indices = np.flatnonzero(source[start:stop]) + start
        integers = first + step * indices
        images = p * integers + shift
        images = images[(images >= first) & (images <= last)]
        if np.any((images - first) % step != 0):
            raise ValueError('The affine image does not belong to the represented arithmetic progression')
        target[(images - first) // step] = True


def signed_terminal_counts(p, q, requests, limit=None, verbose=True, chunk=CHUNK):
    """Retain valid constructions in a fixed window, giving a lower bound without claiming completeness."""
    p, q = integer(p, 'p'), integer(q, 'q')
    chunk = integer(chunk, 'chunk', 1)
    if abs(p) <= 1 or abs(q) <= 1 or gcd(p, q) != 1:
        raise ValueError('The bases must be coprime and have absolute values greater than 1')
    requests = list(requests)
    if not requests:
        return {}
    points = {}
    for h, X in requests:
        h, X = integer(h, 'terminal height', 0), integer(X, 'terminal endpoint', 0)
        points.setdefault(h, set()).add(X)
    largest_X = max(X for values in points.values() for X in values)
    if limit is None:
        limit = largest_X
    limit = integer(limit, 'window radius', 0)
    if limit < largest_X:
        raise ValueError('The window radius is smaller than a required terminal endpoint')
    if 2 * limit > INT64_MAX or abs(p) * limit > INT64_MAX:
        raise ValueError('The window or base exceeds the int64 implementation range')
    result = {(h, X): 0 for h, values in points.items() for X in values}
    if limit == 0 or max(points) == 0:
        return result

    # For (-2,3), represented integers not divisible by 3 are 1 mod 3; the construction preserves this class.
    primitive = (p, q) == (-2, 3)
    step = 3 if primitive else 1
    first = -limit + (1 + limit) % 3 if primitive else -limit
    previous = np.zeros((limit - first) // step + 1, dtype=bool)
    if not primitive:
        previous[-first] = True  # The empty sum 0 assists the construction; terminal counts exclude all multiples of q.
    for h in range(1, max(points) + 1):
        current = previous.copy()
        if primitive:
            current[(1 - first) // 3] = True
        if not primitive or h >= 2:
            add_image(previous, current, p, q**(h - 1),
                      first, step, -limit, limit, chunk=chunk)

        # Multiplication by p preserves q exponents; use increasing absolute-value ranges to complete the closure.
        bound = 1
        while bound <= limit // abs(p):
            upper = min(abs(p) * bound - 1, limit // abs(p))
            add_image(current, current, p, 0, first, step,
                      bound, upper, chunk=chunk)
            add_image(current, current, p, 0, first, step,
                      -upper, -bound, chunk=chunk)
            bound *= abs(p)

        for X in points.get(h, []):
            left = max(0, ceil_div(-X - first, step))
            right = min(len(current), (X - first) // step + 1)
            count = int(np.count_nonzero(current[left:right]))
            if not primitive:
                first_multiple = left + (-(first + left)) % abs(q)
                count -= int(np.count_nonzero(current[first_multiple:right:abs(q)]))
            result[h, X] = count
        previous = current
        if verbose:
            print('  terminal height', h, '/', max(points), flush=True)
    return result


def seed_intervals(grid=POSITIVE_SETTINGS['grid'],
                   denominator=POSITIVE_SETTINGS['denominator'], chunk=CHUNK,
                   seed_exponent=POSITIVE_SETTINGS['seed_exponent'],
                   seed_width=POSITIVE_SETTINGS['seed_width'],
                   max_height=POSITIVE_SETTINGS['height_caps'][0],
                   span=POSITIVE_SETTINGS['seed_span'], verbose=True):
    """Lemma 1 and Appendix A: count minimum heights, take two cumulative sums, then subtract interval endpoints."""
    grid = integer(grid, 'grid', 1)
    denominator = integer(denominator, 'denominator', 1)
    seed_exponent = integer(seed_exponent, 'seed_exponent', 0)
    seed_width = integer(seed_width, 'seed_width', 1)
    max_height = integer(max_height, 'max_height', 0)
    span, chunk = integer(span, 'span', 1), integer(chunk, 'chunk', 1)
    base, last = 5**seed_exponent, 2 * grid
    if denominator > UINT32_MAX or base * denominator > INT64_MAX:
        raise ValueError('Scaling the initial counts exceeds the integer implementation range')
    if (2 * base - 1) * grid + base - 1 > INT64_MAX:
        raise ValueError('The initial interval index calculation exceeds the int64 range')
    heights = odd_heights(2 * base, seed_width, chunk=chunk, verbose=verbose)
    raw = np.zeros((max_height + 1, last + 1), dtype=np.int64)
    for left in range(0, len(heights), chunk):
        right = min(left + chunk, len(heights))
        integers = 2 * np.arange(left, right, dtype=np.int64) + 1
        h = heights[left:right].astype(np.int64)
        keep = (h != 255) & (h <= max_height)
        columns = (integers[keep] * grid + base - 1) // base
        np.add.at(raw, (h[keep], columns), 1)
    raw = np.cumsum(raw, axis=0)
    raw = np.cumsum(raw, axis=1)
    left_columns = np.maximum(np.arange(last + 1) - span, 0)
    counts = raw - raw[:, left_columns]
    return (counts * denominator // base).astype(np.uint32)


def next_intervals(previous, joint, level,
                   grid=POSITIVE_SETTINGS['grid'],
                   size_grid=POSITIVE_SETTINGS['size_grid'],
                   denominator=POSITIVE_SETTINGS['denominator'],
                   P=5**POSITIVE_SETTINGS['width_increment'],
                   height_caps=POSITIVE_SETTINGS['height_caps'], verbose=True,
                   seed_exponent=POSITIVE_SETTINGS['seed_exponent'],
                   width_increment=POSITIVE_SETTINGS['width_increment'],
                   retained_heights=POSITIVE_SETTINGS['retained_heights'],
                   chunk=CHUNK):
    """Proposition 1 and Appendix B: group equal shifts, then translate and add."""
    grid, size_grid = integer(grid, 'grid', 1), integer(size_grid, 'size_grid', 1)
    denominator, P = integer(denominator, 'denominator', 1), integer(P, 'P', 1)
    level, chunk = integer(level, 'level', 2), integer(chunk, 'chunk', 1)
    seed_exponent = integer(seed_exponent, 'seed_exponent', 0)
    width_increment = integer(width_increment, 'width_increment', 1)
    height_caps = tuple(integer(x, 'height_caps', 0) for x in height_caps)
    if level > len(height_caps):
        raise ValueError('level is outside the height_caps range')
    if len(retained_heights) != 2:
        raise ValueError('retained_heights requires a lower and an upper bound')
    first_a, last_a = (integer(x, 'retained height', 1) for x in retained_heights)
    if first_a > last_a:
        raise ValueError('The lower retained height exceeds the upper bound')
    integer_table(previous, 'previous interval table')
    last = 2 * grid
    if previous.shape != (height_caps[level - 2] + 1, last + 1):
        raise ValueError('The previous interval table shape does not match the height and grid')
    if np.any(previous[0]) or np.any(previous[:, 0]):
        raise ValueError('Height 0 and interval index 0 must give zero lower bounds')
    if int(np.max(previous)) > denominator:
        raise ValueError('A previous interval lower bound exceeds denominator')
    if denominator > UINT32_MAX or P > UINT32_MAX or P * denominator > UINT64_MAX:
        raise ValueError('The interval recurrence exceeds the uint32 or uint64 range')
    if not joint:
        raise ValueError('Missing joint residue counts')
    exponent = seed_exponent + width_increment * (level - 1)
    max_height = height_caps[level - 1]
    divisor = size_grid * 5**exponent
    if max_height > 0 and grid * 2**(max_height - 1) > divisor:
        raise ValueError('The condition eta <= 1 required for interval transfer fails')
    table = np.zeros((max_height + 1, last + 1), dtype=np.uint32)
    prefix = {}
    for B, counts in joint.items():
        integer_table(counts, 'joint residue counts')
        if int(np.max(counts)) > P or counts.size * P > UINT64_MAX:
            raise ValueError('The joint counts exceed the safe summation range')
        if int(counts.sum(dtype=np.uint64)) > P:
            raise ValueError('The sum of the joint counts exceeds the modulus P')
        if counts[0, 0] != 1 or np.any(counts[0, 1:]):
            raise ValueError('Height zero must contain exactly one empty chain')
        prefix[B] = np.zeros((len(counts), counts.shape[1] + 1), dtype=np.uint32)
        prefix[B][:, 1:] = np.cumsum(counts, axis=1, dtype=np.uint32)

    for h in range(1, max_height + 1):
        numerator = grid * 2**(h - 1)
        for B, counts in prefix.items():
            edges = [min(counts.shape[1] - 1, ceil_div(j * divisor, numerator))
                     for j in range(last)]
            old_h = min(h, len(previous) - 1)
            total = previous[old_h, :last - 1].astype(np.uint64)
            for a in range(first_a, min(h, last_a, len(counts) - 1) + 1):
                old_h = min(h - a, len(previous) - 1)
                source = previous[old_h, :last - 1].astype(np.uint64)
                if not np.any(source):
                    continue
                weights = counts[a, edges[1:]] - counts[a, edges[:-1]]
                for j in np.flatnonzero(weights):
                    # Processing in chunks only limits the length of temporary multiplication arrays.
                    for left in range(int(j), len(total), chunk):
                        right = min(left + chunk, len(total))
                        total[left:right] += np.uint64(weights[j]) * source[left-j:right-j]
            if np.any(total > P * denominator):
                raise ArithmeticError('The interval recurrence numerator exceeds the proved integer bound')
            table[h, 2:] = np.maximum(table[h, 2:], total // P)
        table[h] = np.maximum(table[h], table[h - 1])
        if verbose and (h % 10 == 0 or h == max_height):
            print('  Step', level, 'interval table, height', h, '/', max_height, flush=True)
    return table


def pack_intervals(table, span):
    """Lemma 5: maximize total weight over disjoint intervals with right endpoint at most i."""
    integer_table(table, 'interval table')
    span = integer(span, 'span', 1)
    if np.any(table[:, 0]):
        raise ValueError('The lower bound at interval index 0 must be zero')
    if int(np.max(table)) * table.shape[1] > UINT64_MAX:
        raise ValueError('Interval packing may exceed the uint64 range')
    U = np.zeros(table.shape, dtype=np.uint64)
    for i in range(table.shape[1]):
        skip = U[:, i - 1] if i > 0 else 0
        before = U[:, i - span] if i >= span else 0
        use = table[:, i].astype(np.uint64) + before
        U[:, i] = np.maximum(skip, use)
    return U


def positive_terminals(scales, joint, verbose=True, chunk=CHUNK):
    settings = POSITIVE_SETTINGS
    chunk = integer(chunk, 'chunk', 1)
    if not scales:
        return {}
    levels = settings['levels']
    if levels != len(settings['height_caps']):
        raise ValueError('The number of positive-case steps does not match the number of height bounds')
    exponent_levels = {settings['seed_exponent'] + settings['width_increment'] * i: i + 1
                       for i in range(levels)}
    rows_by_level = {}
    for row in scales:
        exponent = integer(row['count_exponent'], 'count_exponent', 0)
        if exponent not in exponent_levels:
            raise ValueError('The scale count_exponent does not correspond to a configured step')
        H = integer(row['H'], 'H', 1)
        tail_terms = integer(row['tail_terms'], 'tail_terms', 1)
        if tail_terms > H:
            raise ValueError('The scale record requests a nonpositive terminal height')
        column = F(row['X']) * settings['grid'] / 5**exponent
        if column.denominator != 1 or not 0 <= column <= 2 * settings['grid']:
            raise ValueError('The terminal endpoint does not lie on the configured interval grid')
        rows_by_level.setdefault(exponent_levels[exponent], []).append((row, int(column)))
    previous = seed_intervals(
        grid=settings['grid'], denominator=settings['denominator'], chunk=chunk,
        seed_exponent=settings['seed_exponent'], seed_width=settings['seed_width'],
        max_height=settings['height_caps'][0], span=settings['seed_span'], verbose=verbose)
    result = {}
    for level in range(1, max(rows_by_level) + 1):
        if level > 1:
            previous = next_intervals(
                previous, joint, level, grid=settings['grid'], size_grid=settings['size_grid'],
                denominator=settings['denominator'], P=5**settings['width_increment'],
                height_caps=settings['height_caps'], verbose=verbose,
                seed_exponent=settings['seed_exponent'],
                width_increment=settings['width_increment'],
                retained_heights=settings['retained_heights'], chunk=chunk)
        used = rows_by_level.get(level, [])
        if not used:
            continue
        exponent = settings['seed_exponent'] + settings['width_increment'] * (level - 1)
        counts = pack_intervals(previous, settings['seed_span'] + settings['span_increment'] * (level - 1))
        for row, column in used:
            for v in range(row['tail_terms']):
                h = row['H'] - v
                value = int(counts[min(h, len(counts) - 1), column])
                result[h, row['X']] = F(value * 5**exponent, settings['denominator'])
    return result


def terminal_counts(params, scales, joint=None, verbose=True):
    chunk = integer(params.get('chunk', CHUNK), 'chunk', 1)
    requests = set()
    for row in scales:
        H = integer(row['H'], 'H', 1)
        tail_terms = integer(row['tail_terms'], 'tail_terms', 1)
        if tail_terms > H:
            raise ValueError('The scale record requests a nonpositive terminal height')
        requests.update((H - v, row['X']) for v in range(tail_terms))
    pair = params['p'], params['q']
    if pair == (5, 2):
        if joint is None:
            raise ValueError('Joint counts are missing for the positive terminal calculation')
        return positive_terminals(scales, joint, verbose=verbose, chunk=chunk)
    if pair in ((7, -2), (-7, -2)):
        if requests != {(1, 7)}:
            raise ValueError('The fixed analytic bound for the two cases involving 7 applies only to (h,X)=(1,7)')
        return {(1, 7): 2}  # The distinct singleton representations 1 and p are both indivisible by q.
    if pair not in TERMINAL_LIMITS:
        raise ValueError('No signed terminal window is configured for these bases')
    return signed_terminal_counts(*pair, requests, limit=TERMINAL_LIMITS[pair],
                                  verbose=verbose, chunk=chunk)
