"""Residue families: verify representations and retain the counts needed later.

Read verify_family, then signed_families. The two larger constructions use
a Boolean recurrence and a positive-integer height recurrence. Their NumPy
operations process batches of integers simultaneously.
"""
from fractions import Fraction as F
from math import gcd
import numpy as np
from parameters import CHUNK, POSITIVE_SETTINGS


# These control the computation, rather than the mathematical construction.
PROGRESS_INTERVAL = 4096  # Report progress after this many residues
UNREPRESENTED = 255       # Height marker for an unrepresented integer
INT64_MAX = 2**63 - 1
UINT32_MAX = 2**32 - 1


def check_integer(value, name, minimum=None):
    """Require integer parameters; do not silently truncate floating-point values."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(name + ' must be an integer')
    value = int(value)
    if minimum is not None and value < minimum:
        raise ValueError(name + ' is outside the allowed range')
    return value


def check_bases(p, q):
    p = check_integer(p, 'p')
    q = check_integer(q, 'q')
    if abs(p) < 2 or abs(q) < 2 or gcd(p, q) != 1:
        raise ValueError('The bases must be coprime and have absolute value at least 2')
    return p, q


def check_bounds(bounds):
    """Use integers or Fraction values to keep size bounds exact."""
    bounds = list(bounds)
    if not bounds or any(isinstance(B, bool) or not isinstance(B, (int, F))
                         or B < 0 for B in bounds):
        raise ValueError('Size bounds must be nonnegative integers or exact rational numbers')
    return sorted(set(bounds))


def verify_family(p, q, d, r, terms, B):
    """Check a rational representation and return its height. It need not be retained.

    terms is a sequence of (alpha, gamma) pairs for p**alpha / q**gamma.
    Both exponent sequences increase strictly, alpha < d, and the absolute
    value of the sum is at most B.
    """
    p, q = check_bases(p, q)
    d = check_integer(d, 'd', 1)
    r = check_integer(r, 'r', 0)
    check_bounds([B])
    P = abs(p)**d
    if r >= P:
        raise ValueError('The residue must lie between 0 and P-1')
    terms = tuple(terms)
    value = F(0)                  # The exact rational sum
    residue = 0                   # The same sum modulo P, reduced at the end
    last_a = -1                   # Previous alpha
    last_b = -1                   # Previous gamma
    for a, b in terms:
        a = check_integer(a, 'alpha', 0)
        b = check_integer(b, 'gamma', 0)
        if not (last_a < a < d and last_b < b):
            raise ValueError('The exponent range or strict ordering condition fails')
        value += F(p**a, q**b)
        # pow(q, -b, P) is the inverse of q**b modulo P, not a floating-point power.
        residue += pow(p, a, P) * pow(q, -b, P)
        last_a = a
        last_b = b
    if abs(value) > B or residue % P != r or (r == 0 and terms):
        raise ValueError('The residue witness failed verification')
    return last_b + 1             # Return 0 for the empty chain, otherwise max(gamma)+1


def signed_families(params, verbose=True):
    """Search and verify each residue. hist[B][h] counts residues of height h."""
    p, q = check_bases(params['p'], params['q'])
    d = check_integer(params['d'], 'd', 1)
    Hmax = check_integer(params['H_MAX'], 'H_MAX', 1)
    bounds = check_bounds(params['B'])  # Process distinct bounds B in increasing order
    P = abs(p)**d
    modulus = [abs(p)**e for e in range(d + 1)]
    inverse = [[0] * Hmax]       # At e=0 the modulus is 1, so no inverse is needed
    for e in range(1, d + 1):
        inverse.append([pow(q, -gamma_min, modulus[e]) for gamma_min in range(Hmax)])
    geometry = {}                # Reuse integer-sum bounds for each H

    def divide_interval(lo, hi):
        """Return the closed interval of integers x satisfying lo <= p*x <= hi."""
        if p < 0:
            lo, hi = -hi, -lo
        lower = -((-lo) // abs(p))  # ceil(lo / abs(p))
        upper = hi // abs(p)       # floor(hi / abs(p))
        return lower, upper

    def search(r, B, H):
        """Find a representation with gamma < H; return None on failure."""
        if r == 0:
            return ()
        if H not in geometry:
            # Multiply the rational sum by q**(H-1) to obtain an integer sum.
            powers = [q**(H - 1 - gamma_min) for gamma_min in range(H)]
            # low[e][gamma_min] and high[e][gamma_min]: integer-sum extrema without a congruence condition.
            # The exponent ranges are 0 <= alpha < e and gamma_min <= gamma < H.
            # Include the empty sum 0. No terms are available when e=0 or g=H.
            low = [[0] * (H + 1) for _ in range(d + 1)]
            high = [[0] * (H + 1) for _ in range(d + 1)]
            for e in range(1, d + 1):
                for gamma_min in range(H - 1, -1, -1):
                    # Case 1: no alpha=0 term; divide all terms by p.
                    without_low, without_high = sorted(
                        (p * low[e - 1][gamma_min], p * high[e - 1][gamma_min]))
                    # Case 2: (alpha,gamma)=(0,gamma_min) occurs; remove it and divide by p.
                    with_low, with_high = sorted(
                        (p * low[e - 1][gamma_min + 1], p * high[e - 1][gamma_min + 1]))
                    with_low += powers[gamma_min]
                    with_high += powers[gamma_min]
                    # Case 3: all gamma are greater than g.
                    low[e][gamma_min] = min(without_low, with_low, low[e][gamma_min + 1])
                    high[e][gamma_min] = max(without_high, with_high, high[e][gamma_min + 1])
            geometry[H] = powers, low, high
        powers, low, high = geometry[H]
        failures = set()         # States proved infeasible in this search

        def dfs(e, residue, gamma_min, lo, hi):
            """Search the remaining e p-exponents, with gamma >= g and sum in [lo,hi]."""
            lo = max(lo, low[e][gamma_min])
            hi = min(hi, high[e][gamma_min])
            if lo > hi:
                return None
            if residue == 0 and lo <= 0 <= hi:
                return ()       # The empty chain solves the remaining problem
            if e == 0:
                return None
            key = (e, residue, gamma_min, lo, hi)
            if key in failures:
                return None
            if residue % abs(p) == 0:
                # A residue of 0 modulo p excludes alpha=0: such a term is a unit modulo p.
                lower, upper = divide_interval(lo, hi)
                smaller_residue = (residue // p) % modulus[e - 1]
                answer = dfs(e - 1, smaller_residue, gamma_min, lower, upper)
                if answer is not None:
                    return answer
            else:
                # Exactly one alpha=0 term is required; try each denominator exponent gamma.
                for gamma in range(gamma_min, H):
                    difference = (residue - inverse[e][gamma]) % modulus[e]
                    if difference % abs(p):
                        continue
                    lower, upper = divide_interval(
                        lo - powers[gamma], hi - powers[gamma])
                    smaller_residue = (difference // p) % modulus[e - 1]
                    answer = dfs(e - 1, smaller_residue, gamma + 1, lower, upper)
                    if answer is not None:
                        # We have removed d-e factors of p, so the original alpha is d-e.
                        return ((d - e, gamma),) + answer
            failures.add(key)
            return None

        # The sum is integral, so a nonintegral size bound can be rounded down.
        # This remains exact for Fraction bounds; no float is used.
        size = F(B) * abs(q)**(H - 1)
        integer_bound = size.numerator // size.denominator
        return dfs(d, r, 0, -integer_bound, integer_bound)

    hist = {B: {0: 1} for B in bounds}  # Residue zero uses the empty chain of height 0
    for r in range(1, P):
        best = None
        for B in bounds:
            # The previous chain remains valid as B increases. For the first B, try Hmax.
            if best is None:
                best = search(r, B, Hmax)
                if best is None:
                    raise RuntimeError('A residue has no representation within the search range')
            # Bisect to find the least feasible height: lo is infeasible, hi is feasible.
            lo = 0
            hi = best[-1][1] + 1
            while hi - lo > 1:
                mid = (lo + hi) // 2
                answer = search(r, B, mid)
                if answer is None:
                    lo = mid
                else:
                    best = answer
                    hi = answer[-1][1] + 1
            h = verify_family(p, q, d, r, best, B)
            hist[B][h] = hist[B].get(h, 0) + 1
        # Retain only the height count, not this residue and its full chain.
        if verbose and (r % PROGRESS_INTERVAL == 0 or r == P - 1):
            print('  Checked', r + 1, '/', P, flush=True)
    if any(sum(row.values()) != P for row in hist.values()):
        raise RuntimeError('The height counts have an incorrect total')
    return hist


def boolean_families(params, verbose=True, chunk=CHUNK):
    """When the size bound is automatic, record existence without storing chains."""
    p, q = check_bases(params['p'], params['q'])
    d = check_integer(params['d'], 'd', 1)
    Hmax = check_integer(params['H_MAX'], 'H_MAX', 1)
    chunk = check_integer(chunk, 'chunk', 1)
    bounds = check_bounds(params['B'])
    if min(bounds) < abs(p)**(d - 1):
        raise ValueError('The absolute-value constraint cannot be omitted for these bounds')
    P = abs(p)**d
    if max(abs(p), abs(q)) * P > INT64_MAX:
        raise ValueError('The Boolean recurrence indices may exceed the int64 range')
    # previous[e] stores existence at height H-1; current[e] stores height H.
    previous = [np.zeros(abs(p)**e, dtype=bool) for e in range(d + 1)]
    current = [row.copy() for row in previous]
    for row in previous + current:
        row[0] = True           # The empty chain is allowed at every width and height
    hist = {0: 1}
    total = 1
    for H in range(1, Hmax + 1):
        for e in range(1, d + 1):
            m = abs(p)**e       # Modulus for the current width
            small = abs(p)**(e - 1)
            for left in range(0, m, chunk):
                right = min(left + chunk, m)
                r = np.arange(left, right, dtype=np.int64)
                # Case 1: all gamma>=1; multiplying by q reduces the height by 1.
                answer = previous[e][(q * r) % m].copy()
                # Case 2: all alpha>=1; dividing by p reduces the width by 1.
                divisible = r % abs(p) == 0
                answer[divisible] |= current[e - 1][(r[divisible] // p) % small]
                # Case 3: (0,0) occurs; remove 1, divide by p, then multiply by q.
                has_one = r % abs(p) == 1
                smaller_residue = (q * ((r[has_one] - 1) // p)) % small
                answer[has_one] |= previous[e - 1][smaller_residue]
                current[e][left:right] = answer
        count = int(np.count_nonzero(current[d]))
        if count < total:
            raise RuntimeError('The cumulative height count decreased')
        if count > total:
            hist[H] = count - total
        total = count
        if verbose:
            print('  H =', H, ', covered', total, '/', P, flush=True)
        previous, current = current, previous  # Swap arrays to reuse their storage
        if total == P:
            return {B: dict(hist) for B in bounds}
    raise RuntimeError('Not all residues are covered within the height limit')


def odd_heights(limit, width, chunk=CHUNK, verbose=True):
    """Compute the minimum heights h_min,k(n) of odd n, returning k=width.

    previous[i] = h_min,k-1(2*i+1) is read without modification.
    current[i]  = h_min,k(2*i+1) receives the current results.
    The marker 255 denotes no representation; finite heights must be below 255.
    """
    limit = check_integer(limit, 'limit', 1)
    width = check_integer(width, 'width', 1)
    chunk = check_integer(chunk, 'chunk', 1)
    if limit >= 2**32 or width > 20:
        raise ValueError('The finite height table is outside the supported range')
    # Each array stores one round of heights at one byte per entry; this uses one more table than the earlier version.
    previous = np.full((limit + 1) // 2, UNREPRESENTED, dtype=np.uint8)
    current = np.empty_like(previous)

    for k in range(1, width + 1):
        power = 5**(k - 1)
        # Without the new term 5**(k-1), h_min,k(n) is at most h_min,k-1(n).
        current[:] = previous

        # Read previous and write current, so increasing traversal needs no reversal.
        for left in range((power - 1) // 2, len(previous), chunk):
            right = min(left + chunk, len(previous))
            n = 2 * np.arange(left, right, dtype=np.int64) + 1
            parent = (n - power) // 2

            # parent=2**exponent*odd; handle parent=0 separately below.
            odd = parent.copy()
            exponent = np.zeros(len(parent), dtype=np.int64)
            even = (odd > 0) & (odd % 2 == 0)
            while np.any(even):
                odd[even] //= 2
                exponent[even] += 1
                even = (odd > 0) & (odd % 2 == 0)

            indices = np.maximum(0, (odd - 1) // 2)
            # Read the preceding round: 1 + h_min,k-1((n-5**(k-1))/2).
            candidate = previous[indices].astype(np.int64) + exponent + 1
            candidate[parent == 0] = 1
            candidate = np.minimum(candidate, UNREPRESENTED)
            # Take the minimum with the old height and write the current array.
            current[left:right] = np.minimum(previous[left:right], candidate)

        # Reuse the same two buffers in the next round; swap references without copying.
        previous, current = current, previous
        if verbose:
            print('  Finite height table k =', k, '/', width, flush=True)
    return previous


def positive_families(cap=POSITIVE_SETTINGS['cap'], d=11,
                      bounds=(10,11,12,14,22), hmax=31,
                      size_grid=POSITIVE_SETTINGS['size_grid'],
                      chunk=CHUNK, verbose=True):
    """Search S/2**(h-1) for (5,2), retaining height and joint height/size counts.

    Scan h, then odd S, in increasing order and select the first occurrence
    of each residue. All batches share seen, so every residue is counted once.
    """
    cap = check_integer(cap, 'cap', 1)
    d = check_integer(d, 'd', 1)
    hmax = check_integer(hmax, 'hmax', 1)
    size_grid = check_integer(size_grid, 'size_grid', 1)
    chunk = check_integer(chunk, 'chunk', 1)
    if cap >= 2**32 or d > 20 or hmax >= UNREPRESENTED:
        raise ValueError('Invalid width, height, grid, or chunk size')
    bounds = sorted({check_integer(B, 'B', 1) for B in bounds})
    if not bounds:
        raise ValueError('The positive case requires at least one size bound')
    P = 5**d
    # NumPy multiplication uses fixed-width integers; check bounds before allocation.
    if P > UINT32_MAX or cap * (P - 1) > INT64_MAX or cap * size_grid > INT64_MAX:
        raise ValueError('These parameters may overflow the residue calculation or joint counts')
    if max(bounds) * 2**(hmax - 1) > INT64_MAX:
        raise ValueError('The size comparison exceeds the supported int64 range')
    if chunk * 2 >= P:
        chunk = max(1, (P - 1) // 2)  # Ensure distinct S in a batch give distinct residues
    heights = odd_heights(cap, d, chunk=chunk, verbose=verbose)
    histogram = {}
    joint = {}
    for B in bounds:
        seen = np.zeros(P, dtype=bool)
        seen[0] = True
        # counts[h,b] counts residues of height h with floor(size_grid*sigma)=b.
        counts = np.zeros((hmax + 1, B * size_grid + 1), dtype=np.uint32)
        counts[0, 0] = 1
        covered = 1
        for h in range(1, hmax + 1):
            denominator = 2**(h - 1)
            inverse = pow(denominator, -1, P)
            last = (min(cap, B * denominator) + 1) // 2  # Number of odd integers to inspect
            for left in range(0, last, chunk):
                right = min(left + chunk, last)
                S = 2 * np.arange(left, right, dtype=np.int64) + 1
                residues = (S * inverse) % P
                keep = (heights[left:right] <= h) & ~seen[residues]
                if np.any(keep):
                    values = S[keep]
                    accepted = residues[keep]
                    if np.any(values > B * denominator):
                        raise RuntimeError('The residue family violates its size bound')
                    seen[accepted] = True
                    bins = values * size_grid // denominator
                    added = np.bincount(bins)
                    counts[h, :len(added)] += added.astype(np.uint32)
                    covered += len(accepted)
                if covered == P:
                    break
            if verbose:
                print('  B =', B, ', h =', h, ', covered', covered, '/', P, flush=True)
            if covered == P:
                break
        if covered != P or not np.all(seen):
            raise RuntimeError('Some residues are uncovered; density verification cannot continue')
        row = {}
        for h in range(hmax + 1):
            number = int(counts[h].sum(dtype=np.uint64))
            if number:
                row[h] = number
        if sum(row.values()) != P or row.get(0) != 1:
            raise RuntimeError('The residue family counts are inconsistent')
        histogram[B] = row
        joint[B] = counts
    return histogram, joint


def residue_counts(params, verbose=True):
    """Select an applicable algorithm and return (height counts, joint counts or None)."""
    pair = check_bases(params['p'], params['q'])
    chunk = check_integer(params.get('chunk', CHUNK), 'chunk', 1)
    d = check_integer(params['d'], 'd', 1)
    bounds = check_bounds(params['B'])
    if pair == (5, 2):
        hist, joint = positive_families(
            cap=POSITIVE_SETTINGS['cap'], d=d,
            bounds=bounds, hmax=params['H_MAX'],
            size_grid=POSITIVE_SETTINGS['size_grid'], chunk=chunk, verbose=verbose)
        return hist, joint
    if min(bounds) >= abs(pair[0])**(d - 1):
        return boolean_families(params, verbose=verbose, chunk=chunk), None
    return signed_families(params, verbose=verbose), None
