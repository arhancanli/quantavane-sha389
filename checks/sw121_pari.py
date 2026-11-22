"""SW 2013 section 12.1 control (22 Nov 2026): E0: y^2 + xy = x^3 - x^2 - 4x + 4, p = 223 split multiplicative.
SW print: Reg_p = 153*223^2 + 125*223^3 + 124*223^4 + O(223^5), L = 179*223 + 85*223^2 + 30*223^3 + O(223^4).
Checks PARI's split-multiplicative ('Tate curve') canonical height, as given in the ellpadicheight docstring,
  h = ellpadicheight(E,p,n,P) * [1, -s2 + 1/log(q)/u2]~,
against SW's printed regulator, and the L-invariant log_p(q)/ord_p(q) against SW's printed L.
Usage: python -I sw121_pari.py [curve p gens]"""
import sys, json
import cypari2
pari = cypari2.Pari(size=4 * 10**8, sizemax=1_500_000_000)
ai = sys.argv[1] if len(sys.argv) > 1 else "1,-1,0,-4,4"
p = int(sys.argv[2]) if len(sys.argv) > 2 else 223
n = 6
E = pari.ellinit([int(x) for x in ai.split(",")])
N = int(pari.ellglobalred(E)[0]); ap = int(pari.ellap(E, p))
if len(sys.argv) > 3:
    gens = [[pari(c) for c in P.split(",")] for P in sys.argv[3].split(";")]
else:
    rk = pari.ellrank(E); gens = list(pari.ellsaturation(E, rk[3], 1000))
out = dict(curve=ai, N=N, p=p, ap=ap, pari=str(pari.version()), disc_sign=int(pari.sign(E.disc())),
           tamagawa=int(pari.ellglobalred(E)[2]), torsion=int(pari.elltors(E)[0]), gens=str(gens))
Ep = pari.ellinit(E[:5], pari(f"O({p}^{n + 4})"))
T = pari("(e)->e.tate")(Ep); u2, u, q = T[0], T[1], T[2]
lq = pari.log(q); vq = int(pari.valuation(q, p))
Linv = lq / vq
out["L_inv"] = str(Linv)
s2 = pari.ellpadics2(E, p, n)
col = pari.matrix(2, 1, [1, -s2 + 1 / lq / u2])
r = len(gens)
H = pari.matrix(r, r)
for i in range(r):
    for j in range(r):
        hv = pari.ellpadicheight(E, p, n, gens[i], gens[j])
        H[i, j] = (hv[0] * col[0, 0] + hv[1] * col[1, 0])
Reg_doc = pari.matdet(H)
out["Reg_p_docstring_formula"] = str(Reg_doc)
out["Reg_p_ellpadicregulator"] = str(pari.ellpadicregulator(E, p, n, gens))
lg = pari.log(pari(f"1 + {p} + O({p}^{n + 4})"))
out["Reg_gamma"] = str(Reg_doc / lg ** r)
out["L_over_loggamma"] = str(Linv / lg)
print(json.dumps(out), flush=True)
