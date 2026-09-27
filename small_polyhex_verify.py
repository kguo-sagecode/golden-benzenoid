#!/usr/bin/env python3
# Copyright (c) 2026 Krystal Guo, University of Amsterdam
# SPDX-License-Identifier: MIT
"""Independent verifier for small_polyhex_witnesses.txt.

Checks, from scratch and without trusting the search that produced the file:

  (a) the 116 polyhexes with at most six cells are regenerated here, up to
      translation and the dihedral symmetries of the tiling, the counts
      1, 1, 3, 7, 22, 82 are confirmed, the ring of six cells around a hole
      is confirmed to be amongst them, and the list in the file is compared,
      in order, with the list generated here;

  (b) for each of the 113 polyhexes carrying a witness, the graph is rebuilt
      from the cell list, the witness is checked to be a nonzero integer
      vector supported on vertices of one colour, the integers
      N = y^T y  and  M = -y^T (A^2 - 2I) y  are recomputed from the
      adjacency matrix in exact integer arithmetic and compared with the
      file, and the two strict integer inequalities

              2M - N > 0        and        5 N^2 < (2M - N)^2

      are checked.  Together these say phi N < M, that is
      y^T (A^2 - theta^2 I) y = phi N - M < 0, so the graph has an
      eigenvalue in (-theta, theta) and the polyhex is not in the golden
      gap class.  As a cross-check the same quadratic form is evaluated
      again in exact Q(sqrt5) arithmetic and its sign confirmed negative;

  (c) for each of the three remaining polyhexes, both colour blocks of
      A^2 - theta^2 I are shown to be positive semidefinite by exact
      symmetric elimination over Q(sqrt5), the characteristic polynomial is
      recomputed over Z, the factorization recorded in the file is verified
      by expansion, and a Sturm sequence evaluated exactly at the algebraic
      points -theta and theta confirms that no factor has a root in the open
      interval (-theta, theta); the factor carrying the smallest positive
      root is localized exactly;

  (d) the a priori integer bound is recomputed and reported.

Pure Python 3, standard library only.  No floating point is used anywhere in
the verification path.

Conventions, as everywhere else in this project: the cell (q, r) in axial
coordinates has centre (2q + r, 3r) and vertices at the six offsets
(1,1), (0,2), (-1,1), (-1,-1), (0,-2), (1,-1); two cells are adjacent when
they differ by one of (1,0), (-1,0), (0,1), (0,-1), (1,-1), (-1,1); the
colour of a vertex (x, y) is y mod 3, which takes only the values 1 and 2.

Usage:  python3 small_polyhex_verify.py [small_polyhex_witnesses.txt]
"""
import sys
from fractions import Fraction

DATA = sys.argv[1] if len(sys.argv) > 1 else "small_polyhex_witnesses.txt"

DIRS = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
VOFF = [(1, 1), (0, 2), (-1, 1), (-1, -1), (0, -2), (1, -1)]

FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("  FAIL: " + msg)
    return cond


# --------------------------------------------------------------- the tiling
def cell_vertices(c):
    q, r = c
    cx, cy = 2 * q + r, 3 * r
    return [(cx + a, cy + b) for (a, b) in VOFF]


def graph_of(cells):
    adj = {}
    for c in cells:
        v = cell_vertices(c)
        for i in range(6):
            a, b = v[i], v[(i + 1) % 6]
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
    return adj


def colour(v):
    return ((v[1] % 3) + 3) % 3


def rot(p):
    q, r = p
    return (-r, q + r)


def refl(p):
    q, r = p
    return (r, q)


def normalise(s):
    mq = min(q for q, r in s)
    mr = min(r for q, r in s)
    return tuple(sorted((q - mq, r - mr) for (q, r) in s))


def canon(cells):
    """Lexicographically least of the twelve dihedral images, each
    translated so that the least q and the least r are both zero."""
    best = None
    cur = list(cells)
    for _ in range(6):
        cur = [rot(p) for p in cur]
        for s in (cur, [refl(p) for p in cur]):
            t = normalise(s)
            if best is None or t < best:
                best = t
    return best


def enumerate_small(hmax=6):
    """All polyhexes with at most hmax cells, up to translation and D6."""
    out = {}
    level = {canon([(0, 0)])}
    out[1] = sorted(level)
    for h in range(2, hmax + 1):
        nxt = set()
        for s in level:
            S = set(s)
            for c in s:
                for d in DIRS:
                    z = (c[0] + d[0], c[1] + d[1])
                    if z not in S:
                        nxt.add(canon(list(s) + [z]))
        level = nxt
        out[h] = sorted(level)
    return out


def structure(cells):
    """(vertices, edges, number of holes) of the graph of a polyhex."""
    G = graph_of(cells)
    n = len(G)
    m = sum(len(x) for x in G.values()) // 2
    # a plane graph has m - n + 1 bounded faces; h of them are cells
    return n, m, (m - n + 1) - len(cells)


def colour_block(cells, col):
    """The integer matrix (A^2)[S] for S the colour class col, and S."""
    adj = graph_of(cells)
    V = sorted(v for v in adj if colour(v) == col)
    n = len(V)
    B = [[0] * n for _ in range(n)]
    for i, u in enumerate(V):
        B[i][i] = len(adj[u])
        for j in range(i + 1, n):
            B[i][j] = B[j][i] = len(adj[u] & adj[V[j]])
    return V, B


# ------------------------------------------------------- exact in Q(sqrt5)
class QF(object):
    """a + b sqrt5 with a, b rational."""
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(s, o):
        return QF(s.a + o.a, s.b + o.b)

    def __sub__(s, o):
        return QF(s.a - o.a, s.b - o.b)

    def __mul__(s, o):
        return QF(s.a * o.a + 5 * s.b * o.b, s.a * o.b + s.b * o.a)

    def inv(s):
        d = s.a * s.a - 5 * s.b * s.b
        if d == 0:
            raise ZeroDivisionError
        return QF(s.a / d, -s.b / d)

    def __truediv__(s, o):
        return s * o.inv()

    def sign(s):
        if s.a == 0 and s.b == 0:
            return 0
        if s.a >= 0 and s.b >= 0:
            return 1
        if s.a <= 0 and s.b <= 0:
            return -1
        lhs, rhs = s.a * s.a, 5 * s.b * s.b
        if lhs == rhs:
            return 0
        if s.a > 0:
            return 1 if lhs > rhs else -1
        return -1 if lhs > rhs else 1

    def phi_form(s):
        """Return (a, b) with s = a + b phi, or None if not in Z[phi]."""
        b = 2 * s.b
        a = s.a - s.b
        if a.denominator != 1 or b.denominator != 1:
            return None
        return (int(a), int(b))

    def phi_str(s):
        pf = s.phi_form()
        if pf is None:
            return "%s + %s sqrt5" % (s.a, s.b)
        a, b = pf
        if b == 0:
            return "%d" % a
        t = ("phi" if b == 1 else ("-phi" if b == -1 else "%d phi" % b))
        if a == 0:
            return t
        return "%s %+d" % (t, a)


PHI = QF(Fraction(1, 2), Fraction(1, 2))       # phi   = (1 + sqrt5)/2
THETA = PHI - QF(1)                            # theta = phi - 1 = (sqrt5-1)/2
MTHETA = QF(0) - THETA


def psd_exact(M):
    """Symmetric elimination, pivoting on positive diagonal entries."""
    M = [row[:] for row in M]
    n = len(M)
    while n > 0:
        p = -1
        for i in range(n):
            s = M[i][i].sign()
            if s < 0:
                return False
            if s > 0 and p < 0:
                p = i
        if p < 0:
            for i in range(n):
                for j in range(n):
                    if M[i][j].sign() != 0:
                        return False
            return True
        if p != 0:
            M[0], M[p] = M[p], M[0]
            for i in range(n):
                M[i][0], M[i][p] = M[i][p], M[i][0]
        piv = M[0][0]
        M = [[M[i][j] - M[i][0] * M[0][j] / piv for j in range(1, n)]
             for i in range(1, n)]
        n -= 1
    return True


def block_psd_exact(cells, col):
    """Is (A^2 - theta^2 I)[S] positive semidefinite?   theta^2 = 2 - phi."""
    V, B = colour_block(cells, col)
    n = len(V)
    M = [[QF(B[i][j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        M[i][i] = QF(B[i][i] - 2) + PHI
    return psd_exact(M)


# ------------------------------------------------- polynomials over Q, Sturm
def ptrim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def pmul(a, b):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return ptrim(r)


def pdivmod(a, b):
    a = [Fraction(x) for x in a]
    b = [Fraction(x) for x in b]
    q = [Fraction(0)] * max(len(a) - len(b) + 1, 0)
    while len(a) >= len(b) and a:
        c = a[-1] / b[-1]
        d = len(a) - len(b)
        q[d] = c
        for i, x in enumerate(b):
            a[d + i] -= c * x
        ptrim(a)
        while a and a[-1] == 0:
            a.pop()
    return ptrim(q), a


def pderiv(a):
    return ptrim([a[i] * i for i in range(1, len(a))])


def pgcd(a, b):
    a = [Fraction(x) for x in a]
    b = [Fraction(x) for x in b]
    while b:
        a, b = b, pdivmod(a, b)[1]
    if a:
        a = [x / a[-1] for x in a]
    return a


def sturm_chain(p):
    g = pgcd(p, pderiv(p))
    p = [Fraction(x) for x in p] if len(g) <= 1 else pdivmod(p, g)[0]
    chain = [p, pderiv(p)]
    while len(chain[-1]) > 1:
        r = pdivmod(chain[-2], chain[-1])[1]
        if not r:
            break
        chain.append([-x for x in r])
    return chain


def peval_qf(p, x):
    acc = QF(0)
    for c in reversed(p):
        acc = acc * x + QF(c)
    return acc


def variations(chain, x):
    signs = []
    for p in chain:
        if not p:
            continue
        s = peval_qf(p, x).sign()
        if s:
            signs.append(s)
    return sum(1 for i in range(len(signs) - 1) if signs[i] != signs[i + 1])


def roots_in(p, a, b):
    """Number of distinct real roots of p in the half-open interval (a, b]."""
    ch = sturm_chain(p)
    return variations(ch, a) - variations(ch, b)


def charpoly(A):
    """det(xI - A) over Z, ascending coefficients (Faddeev-LeVerrier)."""
    n = len(A)

    def mul(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)]
                for i in range(n)]

    c = [Fraction(0)] * (n + 1)
    c[n] = Fraction(1)
    M = None
    for k in range(1, n + 1):
        if k == 1:
            M = [[Fraction(1 if i == j else 0) for j in range(n)]
                 for i in range(n)]
        else:
            AM = mul(A, M)
            M = [[AM[i][j] + (c[n - k + 1] if i == j else 0)
                  for j in range(n)] for i in range(n)]
        AM = mul(A, M)
        c[n - k] = -sum(AM[i][i] for i in range(n)) / k
    out = []
    for x in c:
        if x.denominator != 1:
            raise ValueError("characteristic polynomial not integral")
        out.append(int(x))
    return out


def adjacency(cells):
    adj = graph_of(cells)
    V = sorted(adj)
    idx = {v: i for i, v in enumerate(V)}
    n = len(V)
    A = [[0] * n for _ in range(n)]
    for u in V:
        for w in adj[u]:
            A[idx[u]][idx[w]] = 1
    return V, A


# ------------------------------------------------------------- read the file
def read_records(path):
    hdr = {}
    polys = []
    certs = {}
    order = []
    with open(path) as f:
        lines = f.readlines()
    if not lines or lines[0].strip() != "SMALL_POLYHEX_WITNESSES_V1":
        raise SystemExit("bad magic line in %s" % path)
    cur = None
    for ln in lines[1:]:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        t = s.split()
        k = t[0]
        if k in ("COUNTS", "TOTAL", "MEMBERS", "WITNESSED", "MAXENTRY",
                 "MAXN", "MAXM", "MAXINT", "MAXSUPPORT"):
            hdr[k] = [int(x) for x in t[1:]]
        elif k == "POLYHEX":
            pid, h, n, m = (int(x) for x in t[1:5])
            rest = [int(x) for x in t[5:]]
            cells = tuple(sorted((rest[2 * i], rest[2 * i + 1])
                                 for i in range(h)))
            cur = {"id": pid, "h": h, "n": n, "m": m, "cells": cells}
            polys.append(cur)
            order.append(cells)
        elif k == "WITNESS":
            col = int(t[1])
            kk = int(t[2])
            body = [int(x) for x in t[3:]]
            ent = [(body[3 * i], body[3 * i + 1], body[3 * i + 2])
                   for i in range(kk)]
            N, M = body[3 * kk], body[3 * kk + 1]
            cur["witness"] = (col, ent, N, M)
        elif k == "MEMBER":
            cur["member"] = t[1]
        elif k == "MEMBERCERT":
            certs[t[1]] = {"charpoly": [int(x) for x in t[3:]],
                           "factors": [], "minroot": None}
        elif k == "MEMBERFACTOR":
            certs[t[1]]["factors"].append(
                (int(t[2]), [int(x) for x in t[4:]]))
        elif k == "MEMBERMINROOT":
            certs[t[1]]["minroot"] = (int(t[2]), int(t[3]), int(t[4]))
        elif k == "END":
            break
        else:
            raise SystemExit("unknown record %r" % k)
    return hdr, polys, certs, order


# ------------------------------------------------------------------ the run
def main():
    print("=" * 72)
    print("VERIFYING %s" % DATA)
    print("=" * 72)

    hdr, polys, certs, order = read_records(DATA)

    # ---- (a) independent regeneration ------------------------------------
    print()
    print("(a) polyhexes with at most six cells, regenerated independently")
    allp = enumerate_small(6)
    counts = [len(allp[h]) for h in range(1, 7)]
    print("    counts by number of cells: %s   total %d"
          % (" ".join(str(c) for c in counts), sum(counts)))
    check(counts == [1, 1, 3, 7, 22, 82], "counts are not 1 1 3 7 22 82")
    check(sum(counts) == 116, "total is not 116")
    check(hdr.get("COUNTS") == counts, "COUNTS line disagrees")
    check(hdr.get("TOTAL") == [116], "TOTAL line disagrees")

    holed = [c for h in range(1, 7) for c in allp[h] if structure(c)[2] > 0]
    print("    polyhexes with a hole: %d   %s"
          % (len(holed), " ".join(str(c) for c in holed)))
    check(len(holed) == 1 and len(holed[0]) == 6,
          "the six-cell ring around a hole was not found")

    generated = [c for h in range(1, 7) for c in allp[h]]
    check(order == generated,
          "the polyhexes listed in the file differ from those generated here")
    print("    the 116 records of the file agree, in order, with the "
          "regenerated list")

    # ---- (b) integer witnesses -------------------------------------------
    print()
    print("(b) integer witnesses")
    nwit = 0
    max_entry = max_N = max_M = max_int = 0
    max_support = 0
    hist = {}
    worst = None
    for rec in polys:
        cells = rec["cells"]
        n, m, holes = structure(cells)
        check(n == rec["n"] and m == rec["m"],
              "vertex or edge count wrong for %s" % (cells,))
        if "witness" not in rec:
            continue
        nwit += 1
        col, ent, N, M = rec["witness"]
        adj = graph_of(cells)
        V = sorted(adj)
        idx = {v: i for i, v in enumerate(V)}
        # the witness is a nonzero integer vector on one colour class
        seen = set()
        ok = True
        for (x, y, c) in ent:
            v = (x, y)
            ok &= check(v in adj, "witness vertex %s not in the graph of %s"
                        % (v, (cells,)))
            ok &= check(colour(v) == col,
                        "witness vertex %s is not of colour %d" % (v, col))
            ok &= check(v not in seen, "witness vertex %s repeated" % (v,))
            ok &= check(c != 0, "witness entry at %s is zero" % (v,))
            seen.add(v)
        if not ok:
            continue
        check(len(ent) > 0, "empty witness for %s" % (cells,))
        # recompute N and M in exact integer arithmetic from the adjacency
        # matrix:  (A^2)_{uv} = |N(u) cap N(v)|, and (A^2)_{uu} = deg(u)
        yv = {(x, y): c for (x, y, c) in ent}
        NN = sum(c * c for c in yv.values())
        QQ = 0
        keys = sorted(yv)
        for i, u in enumerate(keys):
            QQ += len(adj[u]) * yv[u] * yv[u]
            for j in range(i + 1, len(keys)):
                w = keys[j]
                QQ += 2 * len(adj[u] & adj[w]) * yv[u] * yv[w]
        MM = 2 * NN - QQ
        check(NN == N, "N disagrees for %s: file %d, recomputed %d"
              % (cells, N, NN))
        check(MM == M, "M disagrees for %s: file %d, recomputed %d"
              % (cells, M, MM))
        D = 2 * MM - NN
        check(NN > 0, "N is not positive for %s" % (cells,))
        check(D > 0, "2M - N is not positive for %s" % (cells,))
        check(5 * NN * NN < D * D,
              "5N^2 < (2M-N)^2 fails for %s" % (cells,))
        # cross-check the same form in exact Q(sqrt5):  phi N - M < 0
        val = PHI * QF(NN) - QF(MM)
        check(val.sign() < 0,
              "y^T(A^2 - theta^2 I)y is not negative for %s" % (cells,))
        e = max(abs(c) for c in yv.values())
        if e > max_entry or (e == max_entry and NN > max_N):
            worst = (cells, col, ent, NN, MM)
        max_entry = max(max_entry, e)
        max_N = max(max_N, NN)
        max_M = max(max_M, MM)
        max_support = max(max_support, len(ent))
        hist[e] = hist.get(e, 0) + 1
        max_int = max(max_int, D * D, 5 * NN * NN)
    print("    %d witnesses checked by exact integer arithmetic" % nwit)
    check(nwit == 113, "there are not 113 witnesses")
    check(hdr.get("WITNESSED") == [113], "WITNESSED line disagrees")
    print("    max |y(v)| = %d,  max support = %d,  max N = %d,  max M = %d"
          % (max_entry, max_support, max_N, max_M))
    print("    largest integer occurring in any comparison: %d  (< 2^%d)"
          % (max_int, max_int.bit_length()))
    print("    distribution: %s"
          % ", ".join("%d witnesses with max |y(v)| = %d" % (v, k)
                      for k, v in sorted(hist.items())))
    check(hdr.get("MAXENTRY") == [max_entry], "MAXENTRY line disagrees")
    check(hdr.get("MAXSUPPORT") == [max_support], "MAXSUPPORT line disagrees")
    check(hdr.get("MAXN") == [max_N], "MAXN line disagrees")
    check(hdr.get("MAXM") == [max_M], "MAXM line disagrees")
    check(hdr.get("MAXINT") == [max_int], "MAXINT line disagrees")
    if worst is not None:
        cells, col, ent, NN, MM = worst
        print("    largest witness: cells %s, colour %d, N = %d, M = %d"
              % (" ".join("(%d,%d)" % c for c in cells), col, NN, MM))
        print("        y = " + " ".join("y(%d,%d)=%d" % e for e in ent))

    # a priori bound, from the support size, the largest entry and the
    # degree bound three
    k, c = max_support, max_entry
    bnd_N = k * c * c
    bnd_M = 11 * k * c * c
    bnd = max((2 * bnd_M + bnd_N) ** 2, 5 * bnd_N * bnd_N)
    print("    a priori bound.  The certificate declares support <= %d and"
          " |y(v)| <= %d," % (k, c))
    print("      both checked above.  The graph is subcubic, so every row sum"
          " of A^2")
    print("      is at most 9, whence N <= %d*%d^2 = %d and"
          " |M| <= (2+9)*%d*%d^2 = %d;" % (k, c, bnd_N, k, c, bnd_M))
    print("      every integer compared is then at most %d < 2^%d, far"
          " inside 64 bits." % (bnd, bnd.bit_length()))

    # ---- (c) the three members -------------------------------------------
    print()
    print("(c) the polyhexes in the golden gap class")
    members = [rec for rec in polys if "member" in rec]
    check(len(members) == 3, "there are not exactly 3 members")
    check(hdr.get("MEMBERS") == [3], "MEMBERS line disagrees")
    expected = {
        "benzene": ((0, 0),),
        "naphthalene": ((0, 0), (0, 1)),
        "triphenylene": ((0, 1), (1, 1), (1, 2), (2, 0)),
    }
    for rec in members:
        nm = rec["member"]
        cells = rec["cells"]
        print("    %s: cells %s, %d vertices"
              % (nm, "".join("(%d,%d)" % c for c in cells), rec["n"]))
        check(nm in expected and expected[nm] == cells,
              "%s does not have the expected cells" % nm)
        # exact positive semidefiniteness of both colour blocks
        for col in (1, 2):
            check(block_psd_exact(cells, col),
                  "colour block %d of %s is not PSD" % (col, nm))
        print("        both colour blocks of A^2 - theta^2 I are positive"
              " semidefinite (exact elimination over Q(sqrt5))")
        # the characteristic polynomial and its factorization over Z
        V, A = adjacency(cells)
        cp = charpoly(A)
        cert = certs.get(nm)
        check(cert is not None, "no certificate block for %s" % nm)
        if cert is None:
            continue
        check(cp == cert["charpoly"],
              "characteristic polynomial of %s disagrees with the file" % nm)
        prod = [1]
        for (mult, f) in cert["factors"]:
            for _ in range(mult):
                prod = pmul(prod, f)
        check(prod == cp, "the factorization of %s does not expand to the"
              " characteristic polynomial" % nm)
        print("        chi(x) = " + polystr(cp))
        print("               = " + "".join(
            ("(" + polystr(f) + ")" + ("^%d" % mult if mult > 1 else ""))
            for (mult, f) in cert["factors"]))
        # no factor has a root in the open interval (-theta, theta)
        total_open = 0
        for (mult, f) in cert["factors"]:
            half = roots_in(f, MTHETA, THETA)     # roots in (-theta, theta]
            at = peval_qf(f, THETA)
            isroot = 1 if at.sign() == 0 else 0
            openc = half - isroot
            total_open += openc
            print("        factor %-24s roots in (-theta,theta): %d,"
                  "   value at theta = %s"
                  % (polystr(f), openc, at.phi_str()))
            check(openc == 0, "a factor of %s has a root in (-theta,theta)"
                  % nm)
        check(total_open == 0,
              "%s has an eigenvalue in (-theta,theta)" % nm)
        # the factor carrying the smallest positive root
        fi, hn, hd = cert["minroot"]
        hi = QF(Fraction(hn, hd))
        for j, (mult, f) in enumerate(cert["factors"]):
            c0 = roots_in(f, QF(0), hi)
            if j == fi:
                check(c0 == 1, "the named factor of %s does not have exactly"
                      " one root in (0, %d/%d]" % (nm, hn, hd))
            else:
                check(c0 == 0, "another factor of %s has a root in (0, %d/%d]"
                      % (nm, hn, hd))
        f = cert["factors"][fi][1]
        at = peval_qf(f, THETA)
        # how many roots does the whole characteristic polynomial have in
        # (0, theta] ?
        in0 = sum(roots_in(g, QF(0), THETA) for (mult, g) in cert["factors"])
        if at.sign() == 0:
            check(in0 == 1, "%s should have exactly one root in (0, theta]"
                  % nm)
            print("        smallest positive root: the unique root of %s in"
                  " (0, %d/%d]; it equals theta exactly, and theta is the"
                  " only root of chi in (0, theta]"
                  % (polystr(f), hn, hd))
        else:
            check(in0 == 0, "%s has a root in (0, theta]" % nm)
            print("        smallest positive root: the unique root of %s in"
                  " (0, %d/%d]; chi has no root at all in (0, theta], so"
                  " that root exceeds theta" % (polystr(f), hn, hd))

    # ---- verdict ----------------------------------------------------------
    print()
    print("=" * 72)
    if FAIL:
        print("VERDICT: FAILED, %d check(s) did not pass" % len(FAIL))
        for f in FAIL:
            print("   " + f)
        return 1
    print("VERDICT: all checks passed.  Of the 116 polyhexes with at most "
          "six cells,")
    print("         113 carry a verified integer witness placing an "
          "eigenvalue in")
    print("         (-theta, theta), and exactly 3 lie in the golden gap "
          "class:")
    print("         benzene, naphthalene and triphenylene.  Every integer "
          "compared")
    print("         is at most %d, so 64-bit arithmetic suffices." % max_int)
    print("=" * 72)
    return 0


def polystr(c):
    terms = []
    for k in range(len(c) - 1, -1, -1):
        if c[k] == 0:
            continue
        if k == 0:
            s = "%+d" % c[k]
        elif abs(c[k]) == 1:
            s = "+" if c[k] > 0 else "-"
        else:
            s = "%+d" % c[k]
        terms.append(s + ("x^%d" % k if k > 1 else ("x" if k == 1 else "")))
    return "".join(terms).lstrip("+")


if __name__ == "__main__":
    sys.exit(main())
