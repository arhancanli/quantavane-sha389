"""M2-ss second route (22 Nov 2026): the D-valued p-adic L-function at a supersingular prime (a_p = 0) from CLASSICAL
modular symbols, without PARI's overconvergent machinery (whose memory grows like 3 MB * p^2: p = 107 needs ~35 GB).

Mathematics (all exact identities; written out in ROUTE.md):
  MTT measure  mu_a(a + p^n Z_p) = alpha^-n [a/p^n]^+ - alpha^-(n+1) [a/p^(n-1)]^+,  alpha^2 = -p.
  For zeta of exact order p^m (m >= 1) the character x -> zeta^s(x) (<x> = gamma^s(x), gamma = 1+p) has conductor
  p^(m+1), so the alpha^-(n+1) term drops at level n = m+1 and
        L_alpha(zeta-1) = alpha^-(m+1) Theta_{m+1}(zeta-1),  Theta_n(T) = sum_{a in (Z/p^n)^x} [a/p^n]^+ (1+T)^s(a).
  Pollack (Duke 2003): L_alpha = log+ L+ + alpha log- L-, L+- in Lambda, log+(T) = (1/p) prod_k Phi_{p^2k}(1+T)/p,
  log-(T) = (1/p) prod_k Phi_{p^(2k-1)}(1+T)/p.  m = 1: log- vanishes, log+ = 1/p, so  L+ = -Theta_2  mod Phi_p(1+T);
  m = 2: log+ vanishes, log-(zeta-1) = Phi_p(zeta)/p^2, so  Phi_p(1+T) L- = Theta_3  mod Phi_{p^2}(1+T).
  Phi_{p^k}(1+T) == T^deg mod p (Eisenstein), hence, with c+-_j the T^j coefficients of L+-:
        c+_j == -[T^j] Theta_2 (mod p) for j <= p-2;   c-_j == [T^(p-1+j)] Theta_3 (mod p) for j <= p^2-2p.
  D-valued series (Stein-Wuthrich 2013, 3.5, a_p = 0): L_p = G omega - p H phi(omega), G = log+ L+, H = log- L-;
  if c+-_0 = c+-_1 = 0 the T^2 coefficient is (c+_2/p) omega - (c-_2/p) F omega (PARI basis (omega, F omega), F = p phi),
  and (1-phi)^-2 of it equals [(-c+ + p(c+ + 2c-)), (2c+ + c- - p c-)] / (1+p)^2   (exact; see the route note).
Output: c+-_j mod p for j = 0..3, the predicted leading digits, PARI's regulator (if requested) and the ratio mod p.
Usage: python -I pm_route.py CURVE P [reg] [gens]   CURVE = a1,a2,a3,a4,a6 ; gens like "-1,1;0,0"."""
import sys, time, json
import cypari2
pari = cypari2.Pari(size=2 * 10**8, sizemax=1_200_000_000)
ai = [int(x) for x in sys.argv[1].split(",")]
p = int(sys.argv[2])
want_reg = len(sys.argv) > 3 and sys.argv[3] == "reg"
gens = [[int(c) for c in P.split(",")] for P in sys.argv[4].split(";")] if len(sys.argv) > 4 else []
E = pari.ellinit(ai)
N = int(pari.ellglobalred(E)[0]); ap = int(pari.ellap(E, p))
assert N % p != 0 and ap == 0, (N, ap)
t0 = time.time()
M, xp = pari.msfromell(E, 1)
oo = pari("oo")


def sym_mod_p(num, den):
    v = pari.mseval(M, xp, [oo, pari(f"{num}/{den}")])
    q = pari.denominator(v)
    assert int(q) % p != 0, ("symbol not p-integral", num, den, v)
    return int(pari.numerator(v)) * pow(int(q), -1, p) % p


def slog_table(n):
    """s(a) in [0, p^(n-1)) with <a> = (1+p)^s(a) mod p^n, for a in (Z/p^n)^x."""
    pn, g, tab, x = p ** n, 1 + p, {}, 1
    for s in range(p ** (n - 1)):
        tab[x] = s
        x = x * g % pn
    return tab


def binom_mod_p(s, j):
    r = 1
    while s or j:
        a, b = s % p, j % p
        if b > a:
            return 0
        r = r * (__import__("math").comb(a, b)) % p
        s //= p; j //= p
    return r


def theta_coeffs(n, js):
    pn, tab = p ** n, slog_table(n)
    out = [0] * len(js)
    for a in range(1, pn):
        if a % p == 0:
            continue
        w = pow(a, p ** (n - 1), pn)                      # Teichmueller representative mod p^n
        s = tab[a * pow(w, -1, pn) % pn]
        v = sym_mod_p(a, pn)
        if v:
            for i, j in enumerate(js):
                out[i] = (out[i] + v * binom_mod_p(s, j)) % p
    return out


J = [0, 1, 2, 3]
th2 = theta_coeffs(2, J)
cplus = [(-x) % p for x in th2]
import os
SKIP3 = os.environ.get("PM_SKIP3") == "1"            # level p^3 (c-, the consistency check) costs p^3 symbols
cminus = theta_coeffs(3, [p - 1 + j for j in J]) if not SKIP3 else [None] * 4
res = dict(curve=ai, N=N, p=p, cplus_mod_p=cplus, cminus_mod_p=cminus, sec_symbols=round(time.time() - t0, 1))
cp, cm = cplus[2], cminus[2]
res["lead_digits_omega_Fomega"] = [(-cp) % p, (2 * cp + cm) % p if cm is not None else None]   # (1-phi)^-2 of T^2 coeff mod p
if want_reg:
    n = 4
    R = pari.ellpadicregulator(E, p, n, gens)
    lg = pari.log(pari(f"1 + {p} + O({p}^{n + 4})"))
    Rg = [R[0] / lg ** 2, R[1] / lg ** 2]                            # SW's Reg_gamma in PARI's basis (omega, F omega)
    res["Reg_gamma"] = [str(x) for x in Rg]
    vR = [int(pari.valuation(x, p)) for x in Rg]
    res["v_Reg_gamma"] = vR
    if vR == [0, 0] and cp % p:
        r0 = int(pari.lift(pari.Mod(pari.truncate(Rg[0] + pari(f"O({p})")), p)))
        r1 = int(pari.lift(pari.Mod(pari.truncate(Rg[1] + pari(f"O({p})")), p)))
        S = res["lead_digits_omega_Fomega"][0] * pow(r0, -1, p) % p
        res["S_mod_p"] = S
        if cm is not None:
            res["consistency_second_coord"] = (S * r1 - res["lead_digits_omega_Fomega"][1]) % p == 0
    # red controls: (R1) the non-saturated basis [p*P1, P2] must lose two powers of p in Reg_gamma (no conclusion);
    # (R2) flipping the sign of c-_2 must break the second-coordinate consistency.
    Rb = pari.ellpadicregulator(E, p, n, [pari.ellmul(E, gens[0], p), gens[1]])
    res["red_R1_v_Reg_nonsaturated"] = [int(pari.valuation(x / lg ** 2, p)) for x in Rb]
    if "S_mod_p" in res and cm is not None:
        flip = [(-cp) % p, (2 * cp - cm) % p]
        res["red_R2_flipped_sign_consistent"] = (res["S_mod_p"] * r1 - flip[1]) % p == 0 and (flip[0] - res["S_mod_p"] * r0) % p == 0
    res["rank_bounds_ellrank"] = str(pari.ellrank(E)[:2])
    res["c_infinity"] = 2 if E.disc() > 0 else 1
res["sec_total"] = round(time.time() - t0, 1)
print(json.dumps(res), flush=True)
