"""Control (22 Nov 2026): 91b1, p = 7 split multiplicative, rank 1 (PARI ellpadicbsd docstring example).
Compares PARI's overconvergent ellpadicbsd (an independent algorithm) with split_mult.py's classical [T^2] L_p."""
import json, cypari2
pari = cypari2.Pari(size=4 * 10**8, sizemax=1_500_000_000)
E = pari.ellinit([0, 1, 1, -7, 5]); p = 7; n = 6   # Cremona 91b1 ([-1,3] on it, as in the PARI docstring)
r, Lp = pari.ellpadicbsd(E, p, n)
lg = pari.log(pari(f"1 + {p} + O({p}^{n + 4})"))
gens = [[-1, 3]]
print(json.dumps(dict(curve=str(E[:5]), p=p, bsd_r=int(r), bsd_Lp=str(Lp), Lp_over_loggamma_r=str(Lp / lg ** int(r)),
                      Reg=str(pari.ellpadicregulator(E, p, n, gens)), tamagawa=int(pari.ellglobalred(E)[2]),
                      disc_sign=int(pari.sign(E.disc())))))
# Cross-check with split_mult.py's classical coefficient (level 7^4, exact mod 7^3): [T^2] L_p == 281 (mod 7^3).
# PARI's Lp (s-variable, divided by the L-invariant and by prod c_v / #tors^2 = 1/9) should equal 9 * 281 * lg^2 / L.
Ep = pari.ellinit(E[:5], pari(f"O({p}^{n + 6})"))
q = list(pari("(e)->e.tate")(Ep))[2]
L = pari.log(q) / pari.valuation(q, p)
pred = 9 * 281 * lg ** 2 / L
print(json.dumps(dict(classical_T2_mod_343=281, predicted_pari_Lp_from_classical=str(pred + pari(f"O({p}^4)")),
                      pari_Lp=str(Lp), agree_mod_7_4=bool(pari.valuation(pred - Lp, p) >= 4))))
