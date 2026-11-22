"""Referee control (22 Nov 2026): find curves with Delta > 0 (c_oo = 2), rank 1 or 2 and a small supersingular prime
(a_p = 0, p in 5..13), then compare ref_m2ss.py's leading term with PARI's overconvergent ellpadicbsd (independent
algorithm). This tests the Omega_E = 2 Omega+ normalisation, which program A's controls (c_oo = 1) never exercised."""
import sys, json, itertools
import cypari2
pari = cypari2.Pari(size=2 * 10**8, sizemax=1_200_000_000)
mode = sys.argv[1]
if mode == "search":
    found = []
    for a1, a2, a3, a4, a6 in itertools.product([0, 1], [-1, 0, 1], [0, 1], range(-6, 7), range(-6, 7)):
        try:
            E = pari.ellinit([a1, a2, a3, a4, a6])
        except Exception:
            continue
        if len(E) == 0 or pari.sign(E.disc()) < 0: continue
        gr = pari.ellglobalred(E); N = int(gr[0])
        if N > 2000 or gr[1] != pari("[1,0,0,0]") or int(gr[2]) % 1: continue        # minimal model only
        sp = [p for p in [5, 7, 11, 13] if N % p and pari.ellap(E, p) == 0]
        if not sp: continue
        r = int(pari.ellanalyticrank(E)[0])
        if r == 0: continue
        found.append(dict(a=[a1, a2, a3, a4, a6], N=N, r=r, ss=sp, tam=int(gr[2]), tors=int(pari.elltors(E)[0])))
    found.sort(key=lambda d: (-d["r"], d["ss"][0], d["N"]))
    print(json.dumps(found[:12]))
else:   # mode "bsd": a-invariants, p, n, gens
    ai = [int(x) for x in sys.argv[2].split(",")]; p = int(sys.argv[3]); n = int(sys.argv[4])
    gens = [[pari(c) for c in P.split(",")] for P in sys.argv[5].split(";")]
    E = pari.ellinit(ai); r = len(gens)
    rr, Lp = pari.ellpadicbsd(E, p, n)
    R = pari.ellpadicregulator(E, p, n, gens)
    lg = pari.log(pari(f"1 + {p} + O({p}^{n + 4})"))
    lead = [Lp[0] / lg ** r, Lp[1] / lg ** r]
    print(json.dumps(dict(curve=ai, p=p, r_padic=int(rr), Lp_over_log_r=[str(x) for x in lead],
                          ratio_first=str(Lp[0] / R[0]), ratio_second=str(Lp[1] / R[1]))))
