#!/usr/bin/env python3
"""Exact verification of the finite exact-repair certificates through n=7.

Python 3.10+; standard library only. No LP, floating point, or search imports.
Run `python verify.py`, or add `--json verification_results.json`.

The certificate proves scale * (p D_parent - q D_child) = sum weight * CMI.
Node labels in JSON are 0,...,n-1. Entropy coordinates are message-set masks;
bit positions enumerate ordered pairs (i,j), i != j, lexicographically.

Every coordinate identification is checked by repair closure and an explicit
node permutation. Labels need not be lexicographically minimal: their explicit
permutation witnesses are sufficient. The target is recomputed from n and r.
The supplied `claimed_sum` is only cross-checked, never used as the target.

This program checks the base reconstruction sets and the dependency coverage.
It is an arithmetic/functional-dependency checker, not a proof-assistant kernel.
"""

import argparse
import hashlib
import json
import sys
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations
from pathlib import Path


class VerificationError(ValueError):
    """A certificate or its stated search restriction failed a check."""


def require(condition, message):
    # Explicit checks remain enabled under `python -O`.
    if not condition:
        raise VerificationError(message)


def nonzero(vector):
    return {k: v for k, v in vector.items() if v}


class RepairModel:
    def __init__(self, n):
        require(type(n) is int and 2 <= n <= 7, "unsupported node count")
        self.n = n
        self.edges = [(i, j) for i in range(n) for j in range(n) if i != j]
        self.index = {edge: k for k, edge in enumerate(self.edges)}
        self.full = (1 << len(self.edges)) - 1
        self.rows = [sum(1 << k for k, (i, j) in enumerate(self.edges) if i == v)
                     for v in range(n)]
        self.columns = [sum(1 << k for k, (i, j) in enumerate(self.edges) if j == v)
                        for v in range(n)]

    def edge(self, i, j):
        return 1 << self.index[i, j]

    def check_mask(self, mask):
        require(type(mask) is int and 0 <= mask <= self.full, "invalid message mask")

    @lru_cache(None)
    def close(self, mask):
        self.check_mask(mask)
        while True:
            previous = mask
            for incoming, outgoing in zip(self.columns, self.rows):
                if mask & incoming == incoming:
                    mask |= outgoing
            if previous == mask:
                return mask

    def permute(self, mask, permutation):
        self.check_mask(mask)
        require(isinstance(permutation, list)
                and all(type(v) is int for v in permutation)
                and sorted(permutation) == list(range(self.n)),
                "invalid node permutation")
        return sum(self.edge(permutation[i], permutation[j])
                   for k, (i, j) in enumerate(self.edges) if mask >> k & 1)

    def residual(self, conditioning):
        closed = self.close(conditioning)
        core = 0
        for row in self.rows:
            if closed & row == row:
                core |= row
        return closed & ~core

    def has_directed_cycle(self, mask):
        adjacency = [[] for _ in range(self.n)]
        for k, (i, j) in enumerate(self.edges):
            if mask >> k & 1:
                adjacency[i].append(j)
        color = [0] * self.n

        def visit(v):
            color[v] = 1
            for w in adjacency[v]:
                if color[w] == 1 or (color[w] == 0 and visit(w)):
                    return True
            color[v] = 2
            return False

        return any(color[v] == 0 and visit(v) for v in range(self.n))

    @lru_cache(None)
    def layer_blocks(self, size):
        # One independent single-parity-check block per `size`-node subset.
        return tuple(tuple(sum(self.edge(i, j) for j in nodes if j != i)
                           for i in nodes)
                     for nodes in combinations(range(self.n), size))

    @lru_cache(None)
    def layer_entropy(self, mask, size):
        # A helper message exposes its source symbol in every block containing
        # both endpoints. The entropy of t distinct symbols is min(t, size-1).
        return sum(min(size - 1, sum(bool(mask & source_mask)
                                    for source_mask in block))
                   for block in self.layer_blocks(size))

    def target(self, r, p, q):
        n = self.n
        c, d = r * (r - 1), r * (r + 1)
        # H(empty)=0 is padded to make the literal coefficient sum zero.
        parent = [(self.edge(0, 1), n * (n - 1)),
                  (self.rows[0], n * c), (self.full, -d),
                  (0, -n * (n - 1) - n * c + d)]
        incoming = self.columns[0]
        edge_coefficient = (n - 1) * (n - 2)
        row_coefficient = (n - 1) * c
        child = [(incoming | self.edge(1, 2), edge_coefficient),
                 (incoming | self.rows[1], row_coefficient),
                 (self.full, -d),
                 (incoming, d - edge_coefficient - row_coefficient)]
        return ([(mask, p * coefficient) for mask, coefficient in parent]
                + [(mask, -q * coefficient) for mask, coefficient in child])


def verify_certificate(path):
    certificate = json.loads(path.read_text())
    n, r = certificate['n'], certificate['r']
    require(type(r) is int and 2 <= r <= n - 2, "not an interior certificate")
    model = RepairModel(n)
    p, q, scale = (certificate[k] for k in ('p', 'q', 'scale'))
    require(all(type(v) is int and v > 0 for v in (p, q, scale)),
            "p, q, and scale must be positive integers")
    require(isinstance(certificate['terms'], list) and certificate['terms'],
            "empty certificate")
    total = defaultdict(int)
    depths, arguments = Counter(), Counter()
    cyclic_terms = conditional_entropy_terms = 0

    for t_index, term in enumerate(certificate['terms']):
        weight = term['weight']
        require(type(weight) is int and weight > 0, "nonpositive CMI weight")
        R, X, Y = (term[k] for k in ('R', 'X', 'Y'))
        for mask in (R, X, Y):
            model.check_mask(mask)
        labels, permutations = term['labels'], term['permutations']
        require(len(labels) == len(permutations) == 4, "four CMI corners required")
        corners = [R, R | X | Y, R | X, R | Y]
        for raw, label, permutation, sign in zip(
                corners, labels, permutations, (-1, -1, 1, 1)):
            model.check_mask(label)
            require(model.permute(model.close(raw), permutation) == label,
                    f"incorrect closure/permutation at term {t_index}")
            total[label] += sign * weight

        residual = model.residual(R)
        depth = residual.bit_count()
        require(depth <= n - 2, "conditioning depth exceeds stated search family")
        depths[depth] += 1
        cyclic_terms += model.has_directed_cycle(residual)
        conditional_entropy_terms += X == Y
        for mask in (X, Y):
            kind = ('row' if mask in model.rows else
                    'message' if mask.bit_count() == 1 else 'other')
            require(kind != 'other', "CMI argument is neither a message nor a full row")
            arguments[kind] += 1
        for size in (r, r + 1):
            h = lambda mask: model.layer_entropy(mask, size)
            require(h(R | X) + h(R | Y) - h(R) - h(R | X | Y) == 0,
                    f"term {t_index} is not tight on the size-{size} layer")

    # The archived JSON uses the opposite sign convention: -CMI.
    claimed = {int(k): -v for k, v in certificate['claimed_sum'].items() if v}
    require(nonzero(total) == claimed, "archived claimed_sum differs from recomputation")
    raw_target = model.target(r, p, q)
    witnesses = certificate['target_witnesses']
    require(len(witnesses) == len(raw_target), "wrong number of target witnesses")
    expected = defaultdict(int)
    for (mask, coefficient), witness in zip(raw_target, witnesses):
        require(witness['mask'] == mask and witness['coefficient'] == -coefficient,
                "target coefficient differs from n,r,p,q formula")
        label = witness['label']
        model.check_mask(label)
        require(model.permute(model.close(mask), witness['permutation']) == label,
                "incorrect target closure/permutation")
        expected[label] += scale * coefficient
    require(nonzero(total) == nonzero(expected), "exact target identity failed")
    return {'file': path.name, 'n': n, 'r': r, 'p': p, 'q': q, 'scale': scale,
            'terms': len(certificate['terms']), 'max_conditioning_depth': max(depths),
            'depth_histogram': dict(sorted(depths.items())),
            'cyclic_residual_terms': cyclic_terms,
            'argument_counts': dict(arguments),
            'conditional_entropy_terms': conditional_entropy_terms,
            'both_adjacent_layers_tight': True, 'exact_identity_verified': True}


def verify_bases_and_coverage(records):
    pairs = {(item['n'], item['r']) for item in records}
    require(len(pairs) == len(records), "duplicate (n,r) certificate")
    expected = {(n, r) for n in range(4, 8) for r in range(2, n - 1)}
    require(pairs == expected, "incomplete finite induction: expected all ten cases")
    known = set()
    for n in range(3, 8):
        model = RepairModel(n)
        terminal = model.edge(n - 1, n - 2)
        storage = 0
        for row in model.rows[:n - 2]:
            terminal |= row
        for row in model.rows[:n - 1]:
            storage |= row
        triangle = sum(model.edge(i, j) for i in range(n) for j in range(i))
        require(model.close(terminal) == model.full, "terminal cutset reconstruction failed")
        require(model.close(storage) == model.full, "storage cutset reconstruction failed")
        require(model.close(triangle) == model.full, "bandwidth cutset reconstruction failed")
        require(triangle.bit_count() == n * (n - 1) // 2, "wrong triangle edge count")
        known.add((n, n - 1))
    for n, r in sorted(pairs):
        require((n - 1, r) in known, f"missing child/base for ({n},{r})")
        known.add((n, r))
    for n in (6, 7):
        require(all((n, r) in known for r in range(2, n)), "missing region facet")
    return {'cutset_reconstruction_checked_for_n': list(range(3, 8)),
            'complete_regions': [4, 5, 6, 7],
            'explanation': 'Positive recursion coefficients, every smaller child/base, '
                           'and all endpoint reconstruction sets are present.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', type=Path, help='also save the exact verification report')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / 'manifest.json').read_text())
    records = []
    for item in manifest['certificates']:
        relative = Path(item['file'])
        require(not relative.is_absolute() and '..' not in relative.parts,
                "invalid manifest path")
        path = root / relative
        require(hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'],
                f"certificate checksum mismatch: {relative}")
        record = verify_certificate(path)
        for key in ('n', 'r', 'p', 'q', 'terms', 'scale'):
            require(record[key] == item[key], f"manifest mismatch: {key}")
        records.append(record)
        print(f"PASS n={record['n']} r={record['r']}: "
              f"{record['p']} D_parent >= {record['q']} D_child; "
              f"{record['terms']} positive CMIs; depth <= {record['max_conditioning_depth']}")
    coverage = verify_bases_and_coverage(records)
    print('PASS cutset bases and endpoints, and all finite induction dependencies.')
    print('PASS complete six- and seven-node regions (also n=4,5).')
    print('PASS conditioning grammar and tightness on both adjacent layers.')
    print('All checks used exact integers; no LP solver was used.')
    if args.json:
        args.json.write_text(json.dumps({'certificates': records, 'coverage': coverage},
                                       indent=2) + '\n')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (VerificationError, OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        print(f'FAIL: {exc}', file=sys.stderr)
        sys.exit(1)
