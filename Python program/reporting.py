"""Write exact results and readable tables. These files are outputs, not assumptions used in the proof."""
from fractions import Fraction
from pathlib import Path
import csv
import hashlib
import json
import platform
import sys

import numpy as np


CORE_FILES = ('main.py', 'parameters.py', 'residues.py',
              'terminals.py', 'arithmetic.py', 'reporting.py')


def exact(value):
    """Keep large integers exact in JSON and encode fractions as exact a/b strings."""
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError('Cannot serialize this result type: ' + type(value).__name__)


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('w', encoding='utf-8') as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, default=exact)
        handle.write('\n')
    temporary.replace(path)


def write_csv(path, fields, rows):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerows(rows)


def source_hashes():
    root = Path(__file__).resolve().parent
    return {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
            for name in CORE_FILES}


def environment():
    return {'python': sys.version, 'numpy': np.__version__,
            'platform': platform.platform()}


def residue_report(folder, hist, joint=None, save_arrays=False):
    rows = [(B, h, number) for B in sorted(hist)
            for h, number in sorted(hist[B].items())]
    write_csv(folder / 'residue_heights.csv', ['B', 'height', 'count'], rows)
    write_json(folder / 'residue_heights.json', hist)
    if save_arrays and joint is not None:
        np.savez_compressed(folder / 'joint_counts.npz',
                            **{'B_' + str(B): table for B, table in joint.items()})


def tail_report(folder, w, g, Lambda):
    write_csv(folder / 'tails.csv', ['t', 'upper_bound'], enumerate(w))
    write_json(folder / 'moments.json', {'g': g, 'Lambda': Lambda})


def terminal_report(folder, terminal):
    rows = [(h, X, terminal[h, X]) for h, X in sorted(terminal)]
    write_csv(folder / 'terminal_lower_bounds.csv', ['h', 'X', 'lower_bound'], rows)


def scale_report(folder, scales):
    fields = ['index', 'H', 'X', 'n', 'count_exponent', 'tail_terms', 'a', 'b']
    write_csv(folder / 'scales.csv', fields,
              [(i, *(row[key] for key in fields[1:]))
               for i, row in enumerate(scales, 1)])


def endpoint_report(folder, records):
    # Retain each contribution for an independent check of the final finite inequalities.
    write_json(folder / 'endpoints.json', records)
    fields = ['interval_index', 'left', 'right', 'y', 'value']
    write_csv(folder / 'endpoints.csv', fields,
              [[row[key] for key in fields] for row in records])
