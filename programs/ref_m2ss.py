"""Referee re-computation (22 Nov 2026) of the M2-ss claim: Sha(E/Q)[p^oo] = 0 for E = 389a1 at supersingular p,
written independently of programs/pm_route.py (not imported, not copied).

Differences from program A (pm_route.py), on purpose:
  * symbols: own plus eigensymbol from msinit(N,2,1) + mssplit (Hecke match), scaled by TWISTED L-values
    (Birch: sum_a (D|a) x+(oo->a/D) = sqrt(D) L(E_D,1)/Omega+, D > 0) and then divided by c_oo, so values are SW's
    [r]+ = lambda+(r)/Omega_E (not msfromell, not Omega+). Expected S is then #Sha * (unit) == 1, not 2.
  * Mazur-Tate projection: (Z/p^n)^x enumerated as powers g^k of a primitive root; s(g^k) = k*s(g) mod p^(n-1)
    with s(g) = log_p(g^(p-1)) / ((p-1) log_p(1+p)) (p-adic log), instead of a table + Teichmueller.
  * T-expansion: Horner in F_p[T]/(T^D) instead of Lucas binomials.
  * extra check: [T^j] Theta_3 == 0 mod p for all j < p-1 (forced by Phi_p(1+T) L- == Theta_3 mod Phi_{p^2}).
  * regulator: assembled from ellpadicheight (pairwise [f,g]) via SW Lemma 4.2 / PARI's definition
    [omega,nu]^(r-1)[R,nu] = det<,>_nu, R = Y omega - X eta, converted with ellpadicfrobenius; cross-checked against
    ellpadicregulator. (1-phi)^-2 is applied as a p-adic matrix in the basis (omega, eta), not by a closed formula.
  * S is computed for three functionals: [.,F omega] (program A first coordinate), [.,eta] and [.,omega] (Perrin-Riou's
    component); proportionality (Perrin-Riou) forces all three to agree.
Usage: python -I ref_m2ss.py a1,a2,a3,a4,a6 p r "x1,y1;x2,y2" [levels=2,3] [scale_red=1]
"""
import sys, time, json, math
import cypari2

pari = cypari2.Pari(size=3 * 10**8, sizemax=1_500_000_000)
t00 = time.time()
ai = [int(x) for x in sys.argv[1].split(",")]
p = int(sys.argv[2]); r = int(sys.argv[3])
gens = [[pari(c) for c in P.split(",")] for P in sys.argv[4].split(";")] if r > 0 else []
levels = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else [2, 3]
scale_red = pari(sys.argv[6]) if len(sys.argv) > 6 else pari(1)     # red control: deliberately mis-scaled symbol
out = dict(curve=ai, p=p, r=r, pari=str(pari.version()))
refuse = []

E = pari.ellinit(ai)
N = int(pari.ellglobalred(E)[0]); ap = int(pari.ellap(E, p))
cinf = 2 if pari.sign(E.disc()) > 0 else 1
out.update(N=N, ap=ap, c_inf=cinf, tamagawa=int(pari.ellglobalred(E)[2]), torsion=int(pari.elltors(E)[0]))
if ap != 0: refuse.append("a_p != 0 (method needs a_p = 0)")
if N % p == 0: refuse.append("p | N")
if p < 5: refuse.append("p < 5")
for P in gens:
    if not pari.ellisoncurve(E, P): refuse.append(f"point {P} not on curve")

# ---------- 1. own plus eigensymbol ----------
M = pari.msinit(N, 2, 1)
spaces = pari.mssplit(M, pari.msnew(M))
ells = [l for l in [2, 3, 5, 7, 11, 13, 17, 19, 23] if N % l][:6]
phi = None
for Sp in spaces:
    B = Sp[0] if str(pari.type(Sp)) == "t_VEC" else Sp          # mssplit returns [basis, ...] in 2.17
    if int(pari.matsize(B)[1]) != 1: continue
    if all(pari.mshecke(M, l, Sp)[0, 0] == pari.ellap(E, l) for l in ells):
        phi = pari("(S)->S[,1]")(B)
        break
if phi is None: refuse.append("no rational eigensymbol matching E"); print(json.dumps(out | dict(refuse=refuse))); sys.exit()
den = pari.denominator(pari.content(phi)); phi = phi * den / pari.content(phi * den)   # primitive integral vector
oo = pari("oo")
ms = lambda a, m: pari.mseval(M, phi, [oo, pari(a) / m])
# scale via twisted L-values, Birch: sum_a chi(a) F(oo->a/D) = tau(chi) L(E,chi,1), F = 2 pi i int f dz
Om_plus = pari.real(E.omega()[0])
ratios = []
for D in [5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37, 40, 41, 44]:
    if math.gcd(D, N) != 1 or not pari.isfundamental(D): continue
    s = sum(int(pari.kronecker(D, a)) * ms(a, D) for a in range(1, D) if math.gcd(a, D) == 1)
    if s == 0: continue
    ED = pari.ellinit(pari.elltwist(E, D))
    LD = pari.lfun(ED, 1)
    ratios.append((D, pari.sqrt(D) * LD / Om_plus / s))
    if len(ratios) >= 3: break
if r == 0:
    s0 = ms(0, 1)
    if s0 != 0: ratios.append((1, pari.lfun(E, 1) / Om_plus / s0))
cs = [pari.bestappr(x, 10**6) for _, x in ratios]
out["scale_checks"] = [(D, str(x)[:22]) for D, x in ratios]
if len(cs) < 2 or any(c != cs[0] for c in cs) or any(abs(x - cs[0]) > 1e-15 for _, x in ratios):
    refuse.append("symbol normalisation not pinned consistently"); print(json.dumps(out | dict(refuse=refuse))); sys.exit()
scale = cs[0] / cinf * scale_red         # SW's [r]+ = lambda+(r)/Omega_E = x+ / c_oo
out["scale_xplus_over_phi"] = str(cs[0]); out["scale_used"] = str(scale)
sn, sd = int(pari.numerator(scale)), int(pari.denominator(scale))
if sd % p == 0: refuse.append("normalised symbols not p-integral")
scale_modp = sn * pow(sd, -1, p) % p if sd % p else None
out["sec_setup"] = round(time.time() - t00, 1)

# ---------- 2. Mazur-Tate projections and T-expansion ----------
def mt_poly(n, deg):
    pn, K = p ** n, p ** (n - 1)
    g = int(pari.lift(pari.znprimroot(pn)))
    sg = pari.log(pari(f"{g}^{p - 1} + O({p}^{n + 3})")) / ((p - 1) * pari.log(pari(f"1 + {p} + O({p}^{n + 3})")))
    sg = int(pari.lift(pari.Mod(pari.truncate(sg), K))) if K > 1 else 0
    assert sg % p != 0 or K == 1
    v = [0] * K; a = 1; s = 0
    for k in range((p - 1) * K):
        val = ms(a, pn)
        if val != 0:
            if pari.denominator(val) % p == 0: raise ValueError("raw symbol not p-integral")
            v[s] = (v[s] + int(pari.numerator(val)) * pow(int(pari.denominator(val)), -1, p)) % p
        a = a * g % pn; s = (s + sg) % K
    assert a == 1
    P = [0] * deg                       # Horner: sum_s v_s (1+T)^s in F_p[T]/(T^deg)
    for s in range(K - 1, -1, -1):
        for j in range(deg - 1, 0, -1):
            P[j] = (P[j] + P[j - 1]) % p
        P[0] = (P[0] + v[s]) % p
    return [x * scale_modp % p for x in P]

J = 4
cplus = cminus = None
if scale_modp is not None:
    t = time.time(); th2 = mt_poly(2, J); cplus = [(-x) % p for x in th2]
    out.update(cplus_modp=cplus, sec_level2=round(time.time() - t, 1))
    if 3 in levels:
        t = time.time(); th3 = mt_poly(3, p - 1 + J)
        cminus = th3[p - 1:p - 1 + J]
        low_ok = all(x == 0 for x in th3[:p - 1])
        out.update(cminus_modp=cminus, theta3_low_coeffs_vanish=low_ok, sec_level3=round(time.time() - t, 1))
        if not low_ok: refuse.append("Theta_3 low coefficients nonzero mod p (symbol/projection error)")
    for j in range(min(r, J)):
        if cplus[j] != 0 or (cminus and cminus[j] != 0): refuse.append(f"T^{j} coefficient nonzero mod p but ord_T L_p >= r")

# ---------- 3. D-valued regulator from heights (Lemma 4.2) ----------
prec = 6
Fr = pari.ellpadicfrobenius(E, p, prec + 2)            # matrix on (omega, eta)
lg = pari.log(pari(f"1 + {p} + O({p}^{prec + 4})"))
def to_omF(vec_oe, col=True):
    """coordinates (u,w) with vec = u omega + w F omega, given vec in (omega, eta)."""
    Fo = [Fr[0, 0], Fr[1, 0]] if col else [Fr[0, 0], Fr[0, 1]]
    w = vec_oe[1] / Fo[1]; u = vec_oe[0] - w * Fo[0]
    return [u, w]
if r > 0:
    H = [[pari.ellpadicheight(E, p, prec, gens[i], gens[j]) for j in range(r)] for i in range(r)]
    f = pari.matrix(r, r, [H[i][j][0] for i in range(r) for j in range(r)])
    g = pari.matrix(r, r, [H[i][j][1] for i in range(r) for j in range(r)])
    Y = pari.matdet(f)                     # Reg_eta   (h_eta = f)
    X = pari.matdet(f - g) - Y             # Reg_{omega+eta} - Reg_eta  (h_{omega+eta} = f - g)
    R_oe = [Y / lg ** r, -X / lg ** r]     # Reg_gamma in (omega, eta)
    Rpari = pari.ellpadicregulator(E, p, prec, gens)
    Rp = [Rpari[0] / lg ** r, Rpari[1] / lg ** r]
    for col in (True, False):
        mine = to_omF(R_oe, col)
        out[f"Reg_gamma_omF_mine_{'col' if col else 'row'}"] = [str(x) for x in mine]
    out["Reg_gamma_omF_ellpadicregulator"] = [str(x) for x in Rp]
    mine = to_omF(R_oe, True)
    agree = all(pari.valuation(mine[i] - Rp[i], p) >= 3 for i in range(2))
    out["regulator_assembly_agrees_with_ellpadicregulator(col)"] = bool(agree)
    out["v_Reg_gamma_omF"] = [int(pari.valuation(x, p)) for x in mine]
    out["v_Reg_gamma_omega_eta"] = [int(pari.valuation(x, p)) for x in R_oe]
else:
    R_oe = [pari(1), pari(0)]              # rank 0: SW 4.3, Reg = omega
    mine = [pari(1), pari(0)]

# ---------- 4. (1-phi)^-2 L_p^* as a matrix computation in (omega, eta), and the ratios ----------
def ratios_for(cp, cm):
    I2 = pari.matid(2); phim = Fr / p
    A = (I2 - phim) ** -2
    Fo = pari.matrix(2, 1, [Fr[0, 0], Fr[1, 0]])
    om = pari.matrix(2, 1, [1, 0])
    Lstar = om * pari(cp) / p - Fo * (pari(cm) / p if cm is not None else 0)
    v = A * Lstar
    v_oe = [v[0, 0], v[1, 0]]
    v_omF = to_omF(v_oe, True)
    res = {}
    def rat(x, y):
        if pari.valuation(y, p) != 0: return None
        q = x / y
        if pari.valuation(q, p) < 0: return "nonintegral"
        return int(pari.lift(pari.Mod(pari.truncate(q + pari(f"O({p})")), p)))
    res["S_[.,Fomega]_first_coord"] = rat(v_omF[0], mine[0])
    if r == 0:
        res["second_coord_vanishes_modp"] = pari.valuation(v_omF[1], p) >= 1
    elif cm is not None:
        res["S_[.,omega]_second_coord(PR)"] = rat(v_omF[1], mine[1])
        res["S_[.,eta]_omega_coord_in_(omega,eta)"] = rat(v_oe[0], R_oe[0])
    res["lead_omF_modp"] = [int(pari.lift(pari.Mod(pari.truncate(x + pari(f'O({p})')), p))) if pari.valuation(x, p) >= 0 else "neg" for x in v_omF]
    return res

if cplus is not None:
    k = r if r > 0 else 0
    out["main"] = ratios_for(cplus[k], cminus[k] if cminus else None)
    if cminus:
        out["red_flip_sign_of_H"] = ratios_for(cplus[k], (-cminus[k]) % p)
    S = out["main"]["S_[.,Fomega]_first_coord"]
    vals = [x for kk, x in out["main"].items() if kk.startswith("S_")]
    if r == 0 and not out["main"].get("second_coord_vanishes_modp", True): refuse.append("rank 0: second coordinate nonzero")
    if any(x is None or x == "nonintegral" or x == 0 for x in vals): refuse.append("S not a p-adic unit / regulator coordinate not a unit: no conclusion")
    if len(set(vals)) > 1: refuse.append(f"functionals disagree {vals}: proportionality violated -> error somewhere")
    out["S_modp"] = S
if r > 0:
    Rb = pari.ellpadicregulator(E, p, prec, [pari.ellmul(E, gens[0], p)] + gens[1:])
    out["red_nonsaturated_v_Reg"] = [int(pari.valuation(x / lg ** r, p)) for x in Rb]
out["verdict"] = "REFUSE: " + "; ".join(refuse) if refuse else "PASS: S unit, all functionals agree -> Sha(p) = 0 (given Kato+PR+Pollack)"
out["sec_total"] = round(time.time() - t00, 1)
print(json.dumps(out), flush=True)
