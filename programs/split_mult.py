"""Split multiplicative p (a_p = +1): Sha(E/Q)[p^oo] = 0 test from classical modular symbols (22 Nov 2026).

Mathematics (written out in ROUTE.md, section 6):
  For p || N with a_p = 1, the Mazur-Tate-Teitelbaum measure is mu(a + p^n Z_p) = a_p^-n [a/p^n]+ = [a/p^n]+, so with
  <x> = (1+p)^s(x),  L_p(E,T) = int_{Z_p^x} (1+T)^s(x) dmu  ==  Theta_n(T) := sum_{a in (Z/p^n)^x} [a/p^n]+ (1+T)^s(a)
  modulo omega_{n-1}(T) = (1+T)^(p^(n-1)) - 1. The T^j coefficient of omega_{n-1} has valuation n-1-v_p(j) >= n-1 for
  1 <= j < p, so [T^j] L_p == [T^j] Theta_n (mod p^(n-1)) for j < p (the measure is p-integral: asserted below).
  Exceptional zero: ord_T L_p >= r + 1 (SW 2013, Thm 7.3 with T | d(T), and Thm 6.1).
  SW Thm 6.1 + Thm 7.3: if [T^(r+1)] L_p != 0 then Sha(E/Q)(p) is finite and
      ord_p #Sha(p) <= ord_p [T^(r+1)]L_p - ord_p(L/log_p(1+p)) - sum ord_p(c_v) - ord_p(Reg_gamma) + 2 ord_p #E(Q)(p),
  L = log_p(q_E)/ord_p(q_E) the L-invariant, Reg_gamma = Reg_p / log_p(1+p)^r with the split-multiplicative height of
  SW 4.2 (Schneider, Werner's correction) = PARI's height * [1, -s2 + 1/log(q)/u2] (checked against SW 12.1: sw121_pari.py).
Output: [T^j] L_p mod p^(n-1) for j = 0..4 in PARI's Omega+ normalisation (msfromell; Omega_E = c_oo Omega+), the
L-invariant, Reg_gamma and S = [T^(r+1)]L_p / ((L/log gamma) Reg_gamma) mod p^(n-1); S a unit gives Sha(p) = 0
(S == c_oo prod c_v #Sha mod p is the p-adic BSD prediction).
Usage: python -I split_mult.py a1,a2,a3,a4,a6 p n "x1,y1;x2,y2"   (n = 2 or 3: level p^n)"""
import sys, json, time, math
import cypari2
pari = cypari2.Pari(size=4 * 10**8, sizemax=1_500_000_000)
ai = [int(x) for x in sys.argv[1].split(",")]
p = int(sys.argv[2]); n = int(sys.argv[3])
gens = [[pari(c) for c in P.split(",")] for P in sys.argv[4].split(";")] if len(sys.argv) > 4 and sys.argv[4] else []
r = len(gens)
E = pari.ellinit(ai)
gr = pari.ellglobalred(E)
N = int(gr[0]); ap = int(pari.ellap(E, p))
assert N % p == 0 and N % (p * p) != 0 and ap == 1, ("need split multiplicative p", N, ap)
t0 = time.time()
M, xp = pari.msfromell(E, 1)
oo = pari("oo")
mod = p ** (n - 1)


def sym(num, den):
    v = pari.mseval(M, xp, [oo, pari(f"{num}/{den}")])
    q = int(pari.denominator(v))
    assert q % p != 0, ("symbol not p-integral", num, den)
    return int(pari.numerator(v)) * pow(q, -1, mod) % mod


pn = p ** n
tab, x = {}, 1
for s in range(p ** (n - 1)):          # <a> = (1+p)^s mod p^n
    tab[x] = s
    x = x * (1 + p) % pn
J = 5
coef = [0] * J
for a in range(1, pn):
    if a % p == 0:
        continue
    w = pow(a, p ** (n - 1), pn)       # Teichmueller representative of a mod p^n
    s = tab[a * pow(w, -1, pn) % pn]
    v = sym(a, pn)
    if v:
        for j in range(J):
            coef[j] = (coef[j] + v * math.comb(s, j)) % mod
res = dict(curve=ai, N=N, p=p, level=n, precision=f"mod {p}^{n-1}", r=r, pari=str(pari.version()),
           c_inf=2 if pari.sign(E.disc()) > 0 else 1, tamagawa=int(gr[2]), torsion=int(pari.elltors(E)[0]),
           coeffs_T0_T4_mod=coef, sec_symbols=round(time.time() - t0, 1))
if r:
    prec = 8
    Ep = pari.ellinit(E[:5], pari(f"O({p}^{prec + 4})"))
    u2, u, q = list(pari("(e)->e.tate")(Ep))[:3]
    lq = pari.log(q); vq = int(pari.valuation(q, p))
    Linv = lq / vq
    s2 = pari.ellpadics2(E, p, prec)
    corr = -s2 + 1 / lq / u2
    H = pari.matrix(r, r)
    for i in range(r):
        for j in range(r):
            hv = pari.ellpadicheight(E, p, prec, gens[i], gens[j])
            H[i, j] = hv[0] + hv[1] * corr
    Reg = pari.matdet(H)
    Reg_pari = pari.ellpadicregulator(E, p, prec, gens)
    lg = pari.log(pari(f"1 + {p} + O({p}^{prec + 4})"))
    Rg, Lg = Reg / lg ** r, Linv / lg
    res.update(L_inv=str(Linv), v_L_over_loggamma=int(pari.valuation(Lg, p)), Reg_gamma=str(Rg),
               v_Reg_gamma=int(pari.valuation(Rg, p)), Reg_matches_ellpadicregulator=bool(Reg == Reg_pari),
               rank_bounds_ellrank=str(pari.ellrank(E)[:2]))
    lead = coef[r + 1]
    if all(c == 0 for c in coef[:r + 1]):
        res["low_coeffs_vanish"] = True
    if lead % p and res["v_L_over_loggamma"] == 0 and res["v_Reg_gamma"] == 0:
        den = int(pari.lift(pari.Mod(pari.truncate(Lg * Rg + pari(f"O({p}^{n-1})")), mod)))
        res["S_mod"] = lead * pow(den, -1, mod) % mod
        res["verdict"] = "S unit: Sha(E/Q)(p) = 0 (given SW Thm 6.1, Thm 7.3)"
    else:
        res["verdict"] = "no conclusion at this precision"
    # red control: non-saturated basis [p P1, P2, ...] must raise v(Reg_gamma) by 2
    Rb = pari.ellpadicregulator(E, p, prec, [pari.ellmul(E, gens[0], p)] + gens[1:])
    res["red_nonsaturated_v_Reg_gamma"] = int(pari.valuation(Rb / lg ** r, p))
res["sec_total"] = round(time.time() - t0, 1)
print(json.dumps(res), flush=True)
