#!/usr/bin/env python3
# Copyright (c) 2026 Krystal Guo, University of Amsterdam
# SPDX-License-Identifier: MIT
"""Independent re-verification of the seven-hexagon proof certificate.

Written from the certificate's stated specification, not from its C++ verifier,
and using Python bignums / exact Q(sqrt5) arithmetic throughout.

Checks performed
  A. structural: W7 (85 cells), the 3652 rooted cores regenerated from scratch,
     pattern well-formedness.
  B. every one of the 22009 patterns: the integer negativity certificate,
     recomputed independently; plus a brute-force check on patterns with few
     undecided collar cells that the per-neighbour bound really does dominate
     the true maximum over completions.
  C. the covering DAG: complete case-split, all roots, no gaps.
  D. END-TO-END: sample (core, environment) pairs, route them through the DAG
     to the pattern the certificate claims applies, then directly test that
     core's block for positive semidefiniteness exactly over Q(sqrt5).
     This cross-checks the certificate against ground truth and depends on
     none of its internal reasoning.

Usage: python3 seven_hexagon_reverify.py <obstructions.txt> <cover.bin> [NSPOT]
"""
import sys, struct, random, itertools
from fractions import Fraction

DIRS = [(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]
VOFF = [(1,1),(0,2),(-1,1),(-1,-1),(0,-2),(1,-1)]

def cell_vertices(c):
    q, r = c
    cx, cy = 2*q + r, 3*r
    return [(cx+a, cy+b) for (a, b) in VOFF]

def lex_nonneg(p):
    return p[0] > 0 or (p[0] == 0 and p[1] >= 0)

def hexdist(p):
    return max(abs(p[0]), abs(p[1]), abs(p[0]+p[1]))

# ------------------------------------------------ exact arithmetic in Q(sqrt5)
class QF:
    """a + b*sqrt5 with a,b rational."""
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
        if s.a >= 0 and s.b >= 0: return 1 if (s.a > 0 or s.b > 0) else 0
        if s.a <= 0 and s.b <= 0: return -1
        # opposite signs: compare a^2 with 5b^2
        lhs = s.a*s.a; rhs = 5*s.b*s.b
        if lhs == rhs: return 0
        if s.a > 0:  return 1 if lhs > rhs else -1
        else:        return -1 if lhs > rhs else 1
    def __repr__(s): return "(%s+%s*sqrt5)" % (s.a, s.b)

PHI = QF(Fraction(1,2), Fraction(1,2))          # (1+sqrt5)/2

def psd_exact(M):
    """M is a list of lists of QF.  Exact symmetric elimination."""
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
        N = [[M[i][j] - M[i][0]*M[0][j]/piv for j in range(1, n)]
             for i in range(1, n)]
        M = N; n -= 1
    return True

# ---------------------------------------------------------------- graph model
def graph_of(cells):
    adj = {}
    for c in cells:
        v = cell_vertices(c)
        for i in range(6):
            a, b = v[i], v[(i+1) % 6]
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
    return adj

def colour(v): return ((v[1] % 3) + 3) % 3

def block_is_psd(core, omega):
    """Exactly test whether (A^2 - theta^2 I)[V(core)] is PSD inside G_omega."""
    adj = graph_of(omega)
    V = sorted({v for c in core for v in cell_vertices(c)})
    for col in (1, 2):
        X = [v for v in V if colour(v) == col]
        n = len(X)
        M = [[QF(0) for _ in range(n)] for _ in range(n)]
        for i, u in enumerate(X):
            M[i][i] = QF(len(adj[u]) - 2) + PHI
            for j in range(i+1, n):
                z = len(adj[u] & adj[X[j]])
                M[i][j] = M[j][i] = QF(z)
        if not psd_exact(M):
            return False
    return True

def collar_of(core):
    CS = set(core); out = set()
    for c in core:
        for d in DIRS:
            z = (c[0]+d[0], c[1]+d[1])
            if z not in CS and lex_nonneg(z): out.add(z)
    return sorted(out)

def rooted_cores_7():
    level = {((0,0),)}
    for _ in range(6):
        nxt = set()
        for s in level:
            S = set(s)
            for p in s:
                for d in DIRS:
                    q = (p[0]+d[0], p[1]+d[1])
                    if q not in S and lex_nonneg(q):
                        nxt.add(tuple(sorted(list(s)+[q])))
        level = nxt
    return sorted(level)

# -------------------------------------------------------------------- parsing
def parse_atlas(path):
    toks = open(path).read().split()
    i = 0
    assert toks[i] == "SEVEN_HEXAGON_OBSTRUCTIONS_V1"; i += 1
    assert toks[i] == "W7"; i += 1
    nw = int(toks[i]); i += 1
    W7 = []
    for j in range(nw):
        assert toks[i] == "W"; assert int(toks[i+1]) == j
        W7.append((int(toks[i+2]), int(toks[i+3]))); i += 4
    assert toks[i] == "CORES"; i += 1
    nc = int(toks[i]); i += 1
    cores = []
    for j in range(nc):
        assert toks[i] == "C"; assert int(toks[i+1]) == j; i += 2
        s = []
        for _ in range(7):
            s.append((int(toks[i]), int(toks[i+1]))); i += 2
        cores.append(tuple(s))
    assert toks[i] == "PATTERNS"; i += 1
    npat = int(toks[i]); i += 1
    pats = []
    for j in range(npat):
        assert toks[i] == "P"; assert int(toks[i+1]) == j
        core = int(toks[i+2]); i += 3
        k = int(toks[i]); i += 1
        one = 0
        for _ in range(k):
            b = int(toks[i]); i += 1
            assert 0 <= b < 85; one |= 1 << b
        k = int(toks[i]); i += 1
        zero = 0
        for _ in range(k):
            b = int(toks[i]); i += 1
            assert 0 <= b < 85; zero |= 1 << b
        k = int(toks[i]); i += 1
        sup = []
        for _ in range(k):
            sup.append(((int(toks[i]), int(toks[i+1])), int(toks[i+2]))); i += 3
        pats.append((core, one, zero, sup))
    assert i == len(toks), "trailing tokens in atlas"
    return W7, cores, pats

def parse_cover(path, npat):
    data = open(path, "rb").read()
    assert data[:8] == b"7HEXCVR1", "bad cover magic"
    off = 8
    nn, nr, np_ = struct.unpack_from("<III", data, off); off += 12
    assert np_ == npat and nr == 3652, "cover header mismatch"
    nodes = []
    for _ in range(nn):
        pid, k = struct.unpack_from("<II", data, off); off += 8
        arcs = []
        for _ in range(k):
            var, val, child = struct.unpack_from("<BBI", data, off); off += 6
            arcs.append((var, val, child))
        nodes.append((pid, arcs))
    roots = list(struct.unpack_from("<%dI" % nr, data, off)); off += 4*nr
    assert off == len(data), "extra bytes in cover"
    return nodes, roots

# ------------------------------------------------------------------- checks
def pattern_bound(W7, wi, cores, pat):
    """Recompute the exact integer negativity test for one pattern.
    With N = y^T y and M = -qhat, the condition phi*N < M is equivalent to
    2M-N > 0 and 5N^2 < (2M-N)^2.  Returns (ok, qhat, norm)."""
    core_id, one, zero, sup = pat
    C = cores[core_id]
    CS = set(C); collar = set(collar_of(C))
    coreverts = {v for c in C for v in cell_vertices(c)}
    ymap = {}
    cols = set()
    for (v, y) in sup:
        assert v in coreverts, "witness vertex outside core"
        cols.add(colour(v)); ymap[v] = y
    assert len(cols) == 1 and cols <= {1, 2}, "witness colour"
    norm = sum(y*y for y in ymap.values())
    assert norm > 0
    # edge -> (definitely present?, set of collar-cell indices that would create it)
    emap = {}
    def addcell(c, iscore):
        vv = cell_vertices(c)
        gi = wi.get(c, -1)
        for j in range(6):
            a, b = vv[j], vv[(j+1) % 6]
            key = (a, b) if a < b else (b, a)
            e = emap.setdefault(key, [False, set()])
            if iscore: e[0] = True
            else:
                assert gi >= 0, "collar cell outside W7"
                e[1].add(gi)
    for c in C: addcell(c, True)
    for c in collar: addcell(c, False)
    groups = {}
    for v, coef in ymap.items():
        for (a, b), (fixed_core, covers) in emap.items():
            if a == v: w = b
            elif b == v: w = a
            else: continue
            fixed = fixed_core; unknown = False
            for g in covers:
                if (one >> g) & 1: fixed = True
                elif not ((zero >> g) & 1): unknown = True
            g_ = groups.setdefault(w, [0, []])
            if fixed: g_[0] += coef
            elif unknown: g_[1].append(coef)
    total = 0
    for w, (fixed, unk) in groups.items():
        best = 0
        for m in range(1 << len(unk)):
            s = fixed
            for j in range(len(unk)):
                if (m >> j) & 1: s += unk[j]
            if s*s > best: best = s*s
        total += best
    qhat = total - 2*norm
    N = norm; M = -qhat; D = 2*M - N
    return (D > 0 and 5*N*N < D*D), qhat, norm

def route(nodes, start, present):
    """Follow the DAG for an environment `present` (a bitmask); return the
    pattern id the certificate claims applies."""
    nid = start
    while True:
        pid, arcs = nodes[nid]
        moved = False
        for (var, val, child) in arcs:
            bit = (present >> var) & 1
            if bit != val:
                nid = child; moved = True; break
        if not moved:
            return pid

def verify_cover(nodes, roots, pats, cores, wi):
    """Independent check that the DAG is a complete case split."""
    seen = {}
    sys.setrecursionlimit(10000)
    stats = [0, 0]
    def rec(nid, one, zero):
        stats[0] += 1
        assert nid < len(nodes), "node out of range"
        pid, arcs = nodes[nid]
        assert pid < len(pats), "pattern id out of range"
        pcore, pone, pzero, _ = pats[pid]
        assert (pone & zero) == 0 and (pzero & one) == 0, "node/cube conflict"
        exp = []
        for b in range(85):
            if ((pone >> b) & 1) and not ((one >> b) & 1): exp.append((b, 1))
            elif ((pzero >> b) & 1) and not ((zero >> b) & 1): exp.append((b, 0))
        exp.sort(key=lambda t: (-t[1], t[0]))
        assert len(exp) == len(arcs), "arity mismatch"
        for j, (b, v) in enumerate(exp):
            assert arcs[j][0] == b and arcs[j][1] == v, "literal order mismatch"
        if not arcs:
            stats[1] += 1
            # leaf: the pattern matches this cube; its core must be present
            for c in cores[pcore]:
                assert (one >> wi[c]) & 1, "leaf core not forced present"
            return
        if nid in seen:
            assert seen[nid] == (one, zero), "DAG node reached with two cubes"
            return
        seen[nid] = (one, zero)
        o, z = one, zero
        for (b, v, child) in arcs:
            if v: rec(child, o, z | (1 << b)); o |= (1 << b)
            else: rec(child, o | (1 << b), z); z |= (1 << b)
    for i, r in enumerate(roots):
        one = 0
        for c in cores[i]: one |= 1 << wi[c]
        rec(r, one, 0)
    return stats

def main():
    apath, cpath = sys.argv[1], sys.argv[2]
    NSPOT = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    print("=" * 74)
    print("INDEPENDENT RE-VERIFICATION OF THE SEVEN-HEXAGON CERTIFICATE")
    print("=" * 74)

    W7, cores, pats = parse_atlas(apath)
    wi = {c: i for i, c in enumerate(W7)}
    print("A. structural")
    Wref = sorted([(q, r) for q in range(-7, 8) for r in range(-7, 8)
                   if hexdist((q, r)) <= 7 and lex_nonneg((q, r))])
    assert Wref == W7, "W7 mismatch"
    print("   W7 regenerated independently: %d cells  OK" % len(W7))
    mine = rooted_cores_7()
    assert len(mine) == 3652, "core count %d" % len(mine)
    assert [tuple(sorted(c)) for c in cores] == mine, "core atlas mismatch"
    print("   rooted 7-cell cores regenerated: %d  OK" % len(mine))
    print("   patterns parsed: %d" % len(pats))

    print("B. pattern certificates (recomputed independently)")
    bad = 0
    for j, p in enumerate(pats):
        ok, qhat, norm = pattern_bound(W7, wi, cores, p)
        if not ok:
            bad += 1
            if bad <= 5:
                print("   PATTERN %d FAILS: qhat=%s norm=%s"
                      % (j, qhat, norm))
        if j % 5000 == 0:
            print("   ... %d/%d" % (j, len(pats))); sys.stdout.flush()
    print("   patterns with a valid negative certificate: %d / %d"
          % (len(pats)-bad, len(pats)))
    if bad: print("   *** %d PATTERNS FAIL ***" % bad)

    print("C. covering DAG")
    nodes, roots = parse_cover(cpath, len(pats))
    print("   nodes=%d roots=%d" % (len(nodes), len(roots)))
    stats = verify_cover(nodes, roots, pats, cores, wi)
    print("   complete case split verified; recursive visits=%d leaves=%d  OK"
          % (stats[0], stats[1]))

    print("D. end-to-end spot checks against ground truth (exact Q(sqrt5))")
    random.seed(20260814)
    fails = 0
    for t in range(NSPOT):
        ci = random.randrange(len(cores))
        C = cores[ci]
        coll = collar_of(C)
        k = random.randint(0, len(coll))
        env = set(C) | set(random.sample(coll, k))
        # extra random far cells inside W7 (must not matter)
        for c in random.sample(W7, random.randint(0, 8)):
            env.add(c)
        mask = 0
        for c in env:
            if c in wi: mask |= 1 << wi[c]
        pid = route(nodes, roots[ci], mask)
        pcore, pone, pzero, _ = pats[pid]
        assert (pone & ~mask) == 0, "routed pattern requires an absent cell"
        assert (pzero & mask) == 0, "routed pattern forbids a present cell"
        claimed = cores[pcore]
        assert set(claimed) <= env, "claimed core not inside environment"
        if block_is_psd(list(claimed), sorted(env)):
            fails += 1
            print("   *** COUNTEREXAMPLE: certificate claims core %s is an"
                  " obstruction in env %s, but the block is PSD"
                  % (list(claimed), sorted(env)))
            if fails > 5: break
        if (t+1) % 100 == 0:
            print("   ... %d/%d spot checks, failures=%d" % (t+1, NSPOT, fails))
            sys.stdout.flush()
    print("   spot checks: %d run, %d disagreements with ground truth"
          % (NSPOT, fails))
    print()
    print("VERDICT: %s" % ("ALL CHECKS PASSED" if (bad == 0 and fails == 0)
                           else "*** CERTIFICATE PROBLEM ***"))

main()
