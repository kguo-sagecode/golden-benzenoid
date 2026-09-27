#!/usr/bin/env python3
# Copyright (c) 2026 Krystal Guo, University of Amsterdam
# SPDX-License-Identifier: MIT
"""Exact checks supporting the optimality of seven cells.

Self-contained: pure Python, exact arithmetic in Q(sqrt5) via Fraction, no
third-party dependencies and no data files.

  optimality  (default)  Proposition "six hexagons do not suffice": the
                         ten-cell benzenoid G10 has M_C >= 0 for every
                         connected C of at most six cells, yet is not in the
                         gap class.  Also the smaller seven-cell witness for
                         k <= 5 used in the following remark.

  scan [HMAX]            For every polyhex with 7..HMAX cells, confirm that
                         some connected seven-cell set containing the
                         lexicographically first cell has a negative block.
                         Obstruction verdicts are confirmed exactly; the
                         floating-point step only screens.

Usage:  python3 seven_hexagon_optimality.py [optimality|scan] [HMAX]
"""
import sys
from fractions import Fraction

DIRS = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
VOFF = [(1, 1), (0, 2), (-1, 1), (-1, -1), (0, -2), (1, -1)]

def cell_vertices(c):
    q, r = c
    cx, cy = 2*q + r, 3*r
    return [(cx+a, cy+b) for (a, b) in VOFF]

def graph_of(cells):
    adj = {}
    for c in cells:
        v = cell_vertices(c)
        for i in range(6):
            a, b = v[i], v[(i+1) % 6]
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
    return adj

def colour(v):
    return ((v[1] % 3) + 3) % 3

# ------------------------------------------------ exact arithmetic in Q(sqrt5)
class QF:
    """a + b*sqrt5 with a, b rational."""
    __slots__ = ("a", "b")
    def __init__(self, a=0, b=0):
        self.a = Fraction(a); self.b = Fraction(b)
    def __add__(s, o): return QF(s.a+o.a, s.b+o.b)
    def __sub__(s, o): return QF(s.a-o.a, s.b-o.b)
    def __mul__(s, o): return QF(s.a*o.a + 5*s.b*o.b, s.a*o.b + s.b*o.a)
    def inv(s):
        d = s.a*s.a - 5*s.b*s.b
        if d == 0: raise ZeroDivisionError
        return QF(s.a/d, -s.b/d)
    def __truediv__(s, o): return s * o.inv()
    def sign(s):
        if s.a == 0 and s.b == 0: return 0
        if s.a >= 0 and s.b >= 0: return 1
        if s.a <= 0 and s.b <= 0: return -1
        lhs, rhs = s.a*s.a, 5*s.b*s.b
        if lhs == rhs: return 0
        if s.a > 0: return 1 if lhs > rhs else -1
        return -1 if lhs > rhs else 1

PHI = QF(Fraction(1, 2), Fraction(1, 2))            # (1+sqrt5)/2

def psd_exact(M):
    """Exact symmetric elimination; M a list of lists of QF."""
    M = [row[:] for row in M]
    n = len(M)
    while n > 0:
        p = -1
        for i in range(n):
            s = M[i][i].sign()
            if s < 0: return False
            if s > 0 and p < 0: p = i
        if p < 0:
            for i in range(n):
                for j in range(n):
                    if M[i][j].sign() != 0: return False
            return True
        if p != 0:
            M[0], M[p] = M[p], M[0]
            for i in range(n): M[i][0], M[i][p] = M[i][p], M[i][0]
        piv = M[0][0]
        M = [[M[i][j] - M[i][0]*M[0][j]/piv for j in range(1, n)]
             for i in range(1, n)]
        n -= 1
    return True

def block_is_psd(core, omega):
    """Exactly: is (A(G_omega)^2 - theta^2 I)[V(core)] positive semidefinite?"""
    adj = graph_of(omega)
    V = sorted({v for c in core for v in cell_vertices(c)})
    for col in (1, 2):
        X = [v for v in V if colour(v) == col]
        n = len(X)
        M = [[QF(0) for _ in range(n)] for _ in range(n)]
        for i, u in enumerate(X):
            M[i][i] = QF(len(adj[u]) - 2) + PHI
            for j in range(i+1, n):
                M[i][j] = M[j][i] = QF(len(adj[u] & adj[X[j]]))
        if not psd_exact(M):
            return False
    return True

def connected_subsets(cells, k, must=None):
    out = set()
    for s0 in ([must] if must is not None else list(cells)):
        level = {(s0,)}
        for _ in range(k-1):
            nxt = set()
            for s in level:
                S = set(s)
                for p in s:
                    for d in DIRS:
                        z = (p[0]+d[0], p[1]+d[1])
                        if z in cells and z not in S:
                            nxt.add(tuple(sorted(list(s)+[z])))
            level = nxt
        out |= level
    return sorted(out)

def structure(cells):
    cs = set(cells)
    G = graph_of(cs)
    n = len(G)
    m = sum(len(v) for v in G.values()) // 2
    ncell = {}
    for c in cs:
        for p in cell_vertices(c):
            ncell[p] = ncell.get(p, 0) + 1
    n_i = sum(1 for k in ncell.values() if k == 3)
    return n, m, n_i, (m - n - len(cs) + 1 == 0)

def report(name, cells, kmax):
    n, m, n_i, simply = structure(cells)
    print("  %s: %d cells, %d vertices, %d edges, n_i=%d, simply connected=%s"
          % (name, len(cells), n, m, n_i, simply))
    cs = set(cells)
    print("   size | connected sets | PSD (exact) | negative blocks")
    for k in range(1, kmax+1):
        subs = connected_subsets(cs, k)
        psd = sum(1 for c in subs if block_is_psd(list(c), cells))
        print("    %2d  |      %3d       |     %3d     |     %3d"
              % (k, len(subs), psd, len(subs)-psd))
        sys.stdout.flush()

MODE = sys.argv[1] if len(sys.argv) > 1 else "optimality"
HMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 9

if MODE == "optimality":
    print("=" * 70)
    print("SEVEN IS OPTIMAL: exact checks over Q(sqrt5)")
    print("=" * 70)
    G10 = [(0, 3), (1, 1), (1, 2), (1, 4), (2, 2),
           (2, 3), (3, 0), (3, 1), (3, 3), (4, 1)]
    report("G10", G10, 7)
    print("   => every connected set of at most six cells gives a PSD block;")
    print("      21 of the 23 seven-cell sets give a negative block.")
    print()
    H0 = [(0, 0), (0, 3), (1, 0), (1, 1), (1, 2), (2, -1), (2, 2)]
    report("H0 (smaller witness for k <= 5)", H0, 6)
    print("   => all six five-cell subsets are PSD; all four six-cell")
    print("      subsets are already negative.")

elif MODE == "scan":
    import numpy as np
    print("=" * 70)
    print("Statement scan: polyhexes with 7..%d cells" % HMAX)
    print("=" * 70)
    def rot(p): return (-p[1], p[0]+p[1])
    def refl(p): return (p[1], p[0])
    def canon(cells):
        best = None; cur = list(cells)
        for _ in range(6):
            cur = [rot(p) for p in cur]
            for s in (cur, [refl(p) for p in cur]):
                mq = min(c[0] for c in s); mr = min(c[1] for c in s)
                t = tuple(sorted((q-mq, r-mr) for (q, r) in s))
                if best is None or t < best: best = t
        return best
    PHIF = (1 + 5 ** 0.5) / 2
    def screen_negative(core, omega):
        adj = graph_of(omega)
        V = sorted({v for c in core for v in cell_vertices(c)})
        for col in (1, 2):
            X = [v for v in V if colour(v) == col]
            n = len(X)
            if n == 0: continue
            M = np.zeros((n, n))
            for i, u in enumerate(X):
                M[i, i] = len(adj[u]) - 2 + PHIF
                for j in range(i+1, n):
                    M[i, j] = M[j, i] = len(adj[u] & adj[X[j]])
            if float(np.linalg.eigvalsh(M)[0]) < -1e-9:
                return True
        return False
    level = {canon([(0, 0)])}
    for h in range(1, HMAX+1):
        if h > 1:
            nxt = set()
            for s in level:
                S = set(s)
                for c in s:
                    for d in DIRS:
                        z = (c[0]+d[0], c[1]+d[1])
                        if z not in S: nxt.add(canon(list(s)+[z]))
            level = nxt
        if h < 7: continue
        bad = 0
        for cells in level:
            cells = list(cells)
            root = min(cells)
            cells = sorted((c[0]-root[0], c[1]-root[1]) for c in cells)
            ok = False
            for core in connected_subsets(set(cells), 7, must=(0, 0)):
                if screen_negative(list(core), cells) and \
                        not block_is_psd(list(core), cells):
                    ok = True; break
            if not ok: bad += 1
        print("  h=%d: polyhexes=%-6d without an exactly confirmed"
              " seven-cell obstruction: %d" % (h, len(level), bad))
        sys.stdout.flush()
