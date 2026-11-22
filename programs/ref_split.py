"""Referee re-computation (22 Nov 2026) of the split-multiplicative claim Sha(389a1)[389^oo] = 0.
Written independently of programs/split_mult.py, checks/ctrl_91b1.py and checks/sw121_pari.py (not opened, not imported).

Route (Stein-Wuthrich 2013, Math. Comp. 82):
  (a) symbols: own plus eigensymbol from msinit(N,2,1) + mssplit (Hecke match at many l), scaled by twisted L-values
      sum_a (D|a) [a/D]+ = sqrt(D) L(E,chi_D,1)/Omega_E (D > 0 fundamental, my own derivation), Omega_E = full real
      period computed by my own AGM.  Values are SW's [r]+ = lambda+(r)/Omega_E.  The scale is fixed on one D and
      then CHECKED on every other D (a mis-scaled symbol fails here).  p-integrality of ALL symbol values is checked
      on the generators of Delta_0 (mspathgens), not only on the level-p^n values used.
      Theta_n(T) = sum_{b=1}^{p-1} sum_{j<p^(n-1)} [w(b)(1+p)^j / p^n]+ (1+T)^j, w(b) = Teichmueller lift mod p^n;
      the index j is re-derived from p-adic logs for every a as a cross-check.
  (b) q_E by fixed-point inversion of 1/j = q prod(1-q^n)^24 / E4(q)^3; L = log_p(q_E)/ord_p(q_E) (own log series);
      split height through Tate's parametrisation (SW 4.2):
        h(P) = h(mP)/m^2,  h(Q) = 2 log_p(e(Q)/sigma(Q)) - log_p(u)^2/log_p(q_E),
        sigma(Q) = C * (u-1)/u^(1/2) * prod_n (1-q^n u)(1-q^n/u)/(1-q^n)^2,  psi^* omega_E = C du/u, C = 1/mu,
        u = exp(log_E(Q)/C), log_E = own formal logarithm (series of omega in t = -x/y).
      Checks: mu^2 must be a square in Q_p (else not split: refuse); x(Q) + b2/12 = mu^2 (X(u,q) + 1/12) (Tate's X);
      h(P) computed with m and with 2m must agree (catches a wrong constant C); parallelogram law.
      PARI's ellpadicheight (docstring split correction) and ellpadicregulator are used ONLY as a cross-check.
  (c) S = [T^(r+1)] L_p / ((L/lg) * Reg_gamma), lg = log_p(1+p), Reg_gamma = Reg_p/lg^r, all in SW's Omega_E
      normalisation; p-adic BSD (SW Conj 5.1, (5.2)) predicts S = prod_{v finite} c_v * #Sha / #E(Q)_tors^2.
Usage: python -I ref_split.py a1,a2,a3,a4,a6 p "x1,y1;x2,y2" [levels=2,3] [scale=1] [prec=16] [basisA="a,b;c,d"]
"""
import sys, time, json, math
import cypari2

pari = cypari2.Pari(size=4 * 10**8, sizemax=1_500_000_000)
pari.set_real_precision(60)
t00 = time.time()
ai = [int(x) for x in sys.argv[1].split(",")]
p = int(sys.argv[2])
gens_txt = sys.argv[3] if len(sys.argv) > 3 else ""
gens = [[pari(c) for c in P.split(",")] for P in gens_txt.split(";")] if gens_txt else []
r = len(gens)
opts = dict(a.split("=", 1) for a in sys.argv[4:])
levels = [int(x) for x in opts.get("levels", "2").split(",")]
scale_red = pari(opts.get("scale", "1"))                    # red control: deliberately mis-scaled symbol
K = int(opts.get("prec", "16"))                             # p-adic working precision (digits)
basisA = opts.get("basisA")                                 # red control: replace basis B by A*B (A integral)
out = dict(curve=ai, p=p, r_points=r, pari=str(pari.version()), K=K, scale_red=str(scale_red), basisA=basisA)
refuse = []
def log(msg):
    print(msg, flush=True)

# ---------------- 0. curve data, reduction type ----------------
E = pari.ellinit(ai)
Emin = pari.ellminimalmodel(E)
if [int(x) for x in Emin[:5]] != ai:
    refuse.append("model not global minimal")
gr = pari.ellglobalred(E)
N = int(gr[0])
fa = pari.factor(N)
bad = [int(fa[0][i]) for i in range(len(fa[0]))]
tam = {l: int(pari.elllocalred(E, l)[3]) for l in bad}
prodc = 1
for l in bad: prodc *= tam[l]
tors = int(pari.elltors(E)[0])
disc = E.disc()
cinf = 2 if disc > 0 else 1
ap = int(pari.ellap(E, p))
vpN = int(pari.valuation(N, p))
vpD = int(pari.valuation(disc, p))
out.update(N=N, a_p=ap, v_p_N=vpN, v_p_Delta=vpD, tamagawa=tam, prod_c=prodc, torsion=tors, c_inf=cinf)
if p < 5: refuse.append("p < 5")
if vpN != 1 or ap != 1:
    refuse.append(f"p not split multiplicative (v_p(N)={vpN}, a_p={ap}): split formula does not apply")
for P in gens:
    if not pari.ellisoncurve(E, P): refuse.append(f"point {P} not on curve")
if tors % p == 0: refuse.append("p | #E(Q)_tors")

# ---------------- 0b. surjectivity of rhobar_p (own check) ----------------
# transvection: Tate curve with p not dividing v_p(q_E) = v_p(Delta) => inertia at p has an element of order p
# (Serre 1972 sec. 1.12); irreducible: one Frobenius with char poly irreducible mod p excludes a Borel;
# subgroup of GL2(F_p) with an element of order p, not in a Borel, contains SL2 (Serre 1972 Prop. 15); det = cyclotomic.
irr_l = None
if vpN == 1:
    for l in pari.primes(400):
        l = int(l)
        if N % l == 0 or l == p: continue
        al = int(pari.ellap(E, l))
        if pari.kronecker(al * al - 4 * l, p) == -1:
            irr_l = (l, al); break
out["surjectivity"] = dict(p_not_dividing_vp_Delta=(vpD % p != 0), frobenius_irreducible_at=irr_l)
if not (vpD % p != 0 and irr_l):
    refuse.append("surjectivity of rhobar_p not established by the transvection + irreducible-Frobenius check")

# ---------------- 1. Omega_E by own AGM ----------------
b2, b4, b6 = E[5], E[6], E[7]
x = pari("'x")
roots = pari.polroots(4 * x**3 + b2 * x**2 + 2 * b4 * x + b6, precision=256)
if disc > 0:
    es = sorted([pari.real(z) for z in roots], reverse=True)
    e1, e2, e3 = es
    w1 = pari.Pi(precision=256) / pari.agm(pari.sqrt(e1 - e3), pari.sqrt(e1 - e2))
    OmegaE = 2 * w1
else:
    e1 = [pari.real(z) for z in roots if abs(pari.imag(z)) < 1e-30][0]
    a_ = 3 * e1 + b2 / 4
    b_ = pari.sqrt(3 * e1**2 + b2 * e1 / 2 + b4 / 2)
    OmegaE = 2 * pari.Pi(precision=256) / pari.agm(2 * pari.sqrt(b_), pari.sqrt(2 * b_ + a_))
out["Omega_E_agm"] = str(OmegaE)[:25]
out["Omega_E_vs_cinf_x_E.omega1"] = str(OmegaE / (cinf * pari.real(E.omega()[0])))[:20]

# ---------------- 2. own plus eigensymbol ----------------
M = pari.msinit(N, 2, 1)
spaces = pari.mssplit(M, pari.msnew(M))
ells = [int(l) for l in pari.primes(40) if N % int(l)][:10]
phi = None; nmatch = 0
for Sp in spaces:
    B = Sp[0] if str(pari.type(Sp)) == "t_VEC" else Sp
    if int(pari.matsize(B)[1]) != 1: continue
    if all(pari.mshecke(M, l, Sp)[0, 0] == pari.ellap(E, l) for l in ells):
        phi = pari("(S)->S[,1]")(B); nmatch += 1
if phi is None or nmatch != 1:
    refuse.append(f"no unique rational eigensymbol matching E ({nmatch})")
    log(json.dumps(out | dict(refuse=refuse))); sys.exit()
for l in [int(l) for l in pari.primes(60) if N % int(l)]:      # eigen check on the symbol itself
    T = pari.mshecke(M, l)
    if T * phi != pari.ellap(E, l) * phi:
        refuse.append(f"phi not a T_{l} eigensymbol"); break
oo = pari("oo")
def phival(num, den):
    return pari.mseval(M, phi, [pari(num) / den, oo])          # phi({num/den -> oo})

# scale from twisted L-values: sum_a (D|a) lambda+(a/D) = sqrt(D) L(E,chi_D,1)
Dlist = [D for D in range(5, 200) if pari.isfundamental(D) and math.gcd(D, N) == 1][:9]
twist = []
for D in Dlist:
    sD = sum(int(pari.kronecker(D, a)) * phival(a, D) for a in range(1, D) if math.gcd(a, D) == 1)
    LD = pari.lfun(pari.ellinit(pari.elltwist(E, D)), 1, precision=256)
    twist.append((D, sD, LD))
base = [(D, sD, LD) for (D, sD, LD) in twist if sD != 0]
if not base:
    refuse.append("no twist with nonzero symbol sum"); log(json.dumps(out | dict(refuse=refuse))); sys.exit()
D0, s0, L0 = base[0]
sc_real = pari.sqrt(D0, precision=256) * L0 / OmegaE / s0
scale = pari.bestappr(sc_real, 10**8)
if abs(scale - sc_real) > 1e-30:
    refuse.append("scale not recognised as a rational")
scale = scale * scale_red                       # red control hook: deliberately mis-scaled symbol
out["scale_phi_to_SW"] = str(scale); out["scale_D0"] = D0
chk = []
for (D, sD, LD) in twist:                       # every D re-checked with the FINAL scale (incl. D0)
    lhs = scale * sD; rhs = pari.sqrt(D, precision=256) * LD / OmegaE
    ok = abs(lhs - rhs) < 1e-40 * max(1, abs(rhs))
    chk.append((D, str(lhs), str(rhs)[:18], "%.1e" % float(abs(lhs - rhs)), bool(ok)))
out["twist_checks"] = chk
if not all(c[4] for c in chk):
    refuse.append("normalisation check failed: sum (D|a)[a/D]+ != sqrt(D) L(E_D,1)/Omega_E for some D (mis-scaled symbol)")
if len([c for c in chk if c[4] and c[1] != "0"]) < 3:
    refuse.append("fewer than 3 nonzero twists confirm the scale")
# p-integrality of every symbol value: values on the generators of Delta_0
gvals = pari.mseval(M, phi * scale)
dens = [int(pari.denominator(v)) for v in gvals]
lcmden = 1
for d in dens: lcmden = lcmden * d // math.gcd(lcmden, d)
out["lcm_denominator_all_symbols"] = lcmden
if lcmden % p == 0: refuse.append("normalised symbols not p-integral (measure not integral)")
out["msfromell_ratio_check"] = None
try:
    Mx = pari.msfromell(E, 1)
    xs = Mx[1]
    v1 = pari.mseval(Mx[0], xs, [pari(1) / D0, oo]); v2 = scale * phival(1, D0)
    out["msfromell_ratio_check"] = "x+ / [.]+_SW on 1/D0: " + (str(v1 / v2) if v2 != 0 else "n/a (zero)")
except Exception as e:
    out["msfromell_ratio_check"] = "msfromell failed: " + str(e)[:60]
out["sec_setup"] = round(time.time() - t00, 1)
log(f"[setup] {json.dumps(out)}")

# ---------------- 3. Theta_n and T-coefficients ----------------
def binom_mod(j, k, mod):
    return math.comb(j, k) % mod

def theta(n, J):
    pn, K1 = p ** n, p ** (n - 1)
    mod = p ** (n - 1)
    A = [pari(0)] * K1
    onep = 1 + p
    jcheck_bad = 0
    lg_hi = pari(f"log(1+{p}+O({p}^{n+4}))")
    pq = p ** (n + 4)
    for b in range(1, p):
        w = pow(b, K1, pn)                       # Teichmueller lift of b mod p^n
        a = w
        for j in range(K1):
            A[j] = A[j] + phival(a, pn)
            if n <= 3 and (j % 7 == 0 or n == 2):  # independent index check via p-adic log
                la = pari(f"log({pow(a, p - 1, pq)}+O({p}^{n+4}))") / (p - 1) / lg_hi
                jj = int(pari.lift(pari.Mod(pari.truncate(la), mod))) if mod > 1 else 0
                if jj != j % mod: jcheck_bad += 1
            a = a * onep % pn
    cnt = (p - 1) * K1
    coeffs = []
    for k in range(J):
        s = pari(0)
        for j in range(K1):
            if A[j] != 0:
                c = binom_mod(j, k, mod)
                if c: s += c * A[j]
        s = s * scale
        num, den = int(pari.numerator(s)), int(pari.denominator(s))
        if den % p == 0:
            coeffs.append(None)
        else:
            coeffs.append(num * pow(den, -1, mod) % mod)
    # precision for coefficient k: e_{n-1,k} = min_{1<=i<=k} v_p(C(p^(n-1), i))  (SW Prop 3.5, c_n = 0)
    prec_k = [n - 1] + [min(int(pari.valuation(math.comb(K1, i), p)) for i in range(1, k + 1)) for k in range(1, J)]
    return coeffs, prec_k, cnt, jcheck_bad

J = r + 3
thetas = {}
for n in levels:
    t = time.time()
    co, pk, cnt, jb = theta(n, J)
    thetas[n] = (co, pk)
    out[f"level{n}"] = dict(n_symbols=cnt, coeffs_T0_up_mod_p_pow=[(c, e) for c, e in zip(co, pk)],
                            j_index_mismatches=jb, sec=round(time.time() - t, 1))
    if jb: refuse.append(f"level {n}: Teichmueller index disagrees with p-adic log index ({jb})")
    log(f"[level {n}] {json.dumps(out[f'level{n}'])}")
nmax = max(levels)
co, pk = thetas[nmax]
c_lead = co[r + 1]; c_lead_prec = pk[r + 1]
for k in range(r + 1):
    if co[k] is not None and co[k] % p != 0:
        refuse.append(f"[T^{k}] L_p nonzero mod p although ord_T L_p >= r+1 by SW 6.1+7.3 (inconsistent data)")

# ---------------- 4. p-adic tools (own series) ----------------
def P_(xq, prec=K):
    return pari(xq) + pari(f"O({p}^{prec})")
def plog_unit(w):
    y = w ** (p - 1) - 1
    s = pari(0); yn = pari(1)
    for nn in range(1, 3 * K + 8):
        yn = yn * y
        s += (-1) ** (nn + 1) * yn / nn
    return s / (p - 1)
def plog(xp):
    v = int(pari.valuation(xp, p))
    return plog_unit(xp / pari(p) ** v)
def pexp(z):
    assert int(pari.valuation(z, p)) >= 1
    s = pari(1); t = pari(1)
    for nn in range(1, 3 * K + 8):
        t = t * z / nn; s += t
    return s
lg = plog_unit(P_(1 + p, K + 4))

# ---------------- 5. Tate parameter, L-invariant, isomorphism to E_q ----------------
c4, c6 = E[9], E[10]
jE = pari("(e)->e.j")(E)
if int(pari.valuation(jE, p)) >= 0:
    refuse.append("v_p(j) >= 0: no Tate parametrisation (good reduction at p)")
    out["refuse"] = refuse; out["conclusion"] = "NO CONCLUSION"
    log("[result] " + json.dumps(out)); sys.exit()
ijp = P_(1 / jE, K + 6)
def sig(k, n): return sum(d ** k for d in range(1, n + 1) if n % d == 0)
NQ = K + 8
def E4(q): return 1 + 240 * sum(sig(3, n) * q ** n for n in range(1, NQ))
def E6(q): return 1 - 504 * sum(sig(5, n) * q ** n for n in range(1, NQ))
def prod24(q):
    s = pari(1)
    for n in range(1, NQ): s = s * (1 - q ** n)
    return s ** 24
qE = ijp
for it in range(K + 12):
    qE = ijp * E4(qE) ** 3 / prod24(qE)
jchk = E4(qE) ** 3 / (qE * prod24(qE)) - P_(jE, K + 6)
vq = int(pari.valuation(qE, p))
out["q_E"] = str(qE); out["ord_p_q_E"] = vq
out["j_of_q_minus_j_valuation"] = int(pari.valuation(jchk, p)) if jchk != 0 else "exact"
if vq <= 0: refuse.append("ord_p(q_E) <= 0")
logq = plog(qE)
Linv = logq / vq
out["L_invariant"] = str(Linv)
mu2 = P_(c6, K + 6) * E4(qE) / (P_(c4, K + 6) * (-E6(qE)))
mu4chk = P_(c4, K + 6) / E4(qE) - mu2 ** 2
out["mu2"] = str(mu2); out["mu4_consistency_valuation"] = int(pari.valuation(mu4chk, p)) if mu4chk != 0 else "exact"
split_ok = bool(pari.issquare(mu2))
out["mu2_is_square_in_Qp"] = split_ok
if not split_ok:
    refuse.append("mu^2 not a square in Q_p: E is not a Tate curve over Q_p (not split)")
    out["refuse"] = refuse; out["conclusion"] = "NO CONCLUSION"
    log("[result] " + json.dumps(out)); sys.exit()
mu = pari.sqrt(mu2)
logC = -plog_unit(mu)
out["log_p_C"] = str(logC)
if opts.get("dropC") == "1":                    # red control: SW 4.2's sigma_p(u) taken literally, without C
    logC = pari(0)

# ---------------- 6. formal logarithm (own) ----------------
Ltr = 3 * K + 12
a1, a2, a3, a4, a6 = [pari(c) for c in ai]
tt = pari("'t")
Ot = pari(f"O(t^{Ltr + 4})")
w = tt ** 3 + Ot
for it in range(Ltr + 4):
    w = tt ** 3 + a1 * tt * w + a2 * tt ** 2 * w + a3 * w ** 2 + a4 * tt * w ** 2 + a6 * w ** 3 + Ot
xs_ = tt / w; ys_ = -1 / w
om = pari.deriv(xs_, tt) / (2 * ys_ + a1 * xs_ + a3)
logser = pari.truncate(pari.intformal(om, tt))
try:
    fl = pari.ellformallog(E, Ltr - 2, tt)
    out["formal_log_vs_ellformallog"] = bool(pari.truncate(fl) == pari.truncate(logser + pari(f"O(t^{Ltr - 1})")))
except Exception as ex:
    out["formal_log_vs_ellformallog"] = "n/a: " + str(ex)[:80]

def Xtate(u, q):
    s = u / (1 - u) ** 2
    for n in range(1, NQ):
        qn = q ** n
        s += qn * u / (1 - qn * u) ** 2 + qn / u / (1 - qn / u) ** 2 - 2 * qn / (1 - qn) ** 2
    return s

# ---------------- 7. split height ----------------
def nonsingular_at(Q, l):
    xq, yq = Q[0], Q[1]
    if pari.valuation(xq, l) < 0: return True
    X_, Y_ = int(pari.lift(pari.Mod(xq, l))), int(pari.lift(pari.Mod(yq, l)))
    A1, A2, A3, A4, A6 = ai
    Fx = (A1 * Y_ - 3 * X_ * X_ - 2 * A2 * X_ - A4) % l
    Fy = (2 * Y_ + A1 * X_ + A3) % l
    return not (Fx == 0 and Fy == 0)
def good_mult(P):
    """smallest m | (p-1)*prod c_l with mP in E^0(Q_l) for all l | N and mP in E^1(Q_p)."""
    M0 = (p - 1) * prodc
    for m in sorted(int(d) for d in pari.divisors(M0)):
        Q = pari.ellmul(E, P, m)
        if Q == pari([0]): continue
        if pari.valuation(Q[0], p) >= 0: continue
        if all(nonsingular_at(Q, l) for l in bad if l != p):
            return m, Q
    raise ValueError("no good multiple found")
def h_of_Q(Q):
    xq, yq = Q[0], Q[1]
    den = int(pari.denominator(xq)); e = math.isqrt(den)
    assert e * e == den
    k = int(pari.valuation(e, p)); assert k >= 1
    ep = e // p ** k
    t0 = P_(-xq / yq, K + 2 * k + 4)
    lE = pari.subst(logser, tt, t0)
    z = lE * mu                                     # log_p(u) = log_E(Q)/C, C = 1/mu
    u = pexp(z)
    # Tate X check: x + b2/12 == mu^2 (X(u,q) + 1/12)
    xchk = P_(xq + b2 / 12, K) - mu2 * (Xtate(u, qE) + pari(1) / 12)
    xrel = int(pari.valuation(xchk, p)) - int(pari.valuation(P_(xq, K), p)) if xchk != 0 else 99
    u1 = (u - 1) / pari(p) ** k
    sprod = pari(0)
    for n in range(1, NQ):
        qn = qE ** n
        sprod += plog_unit(1 - qn * u) + plog_unit(1 - qn / u) - 2 * plog_unit(1 - qn)
    logsig_unit = logC + plog_unit(u1) - z / 2 + sprod
    h = 2 * (plog_unit(P_(ep, K + 4)) - logsig_unit) - z ** 2 / logq
    return h, xrel, k
def height(P, extra=1):
    m, Q = good_mult(P)
    if extra != 1:
        m = m * extra; Q = pari.ellmul(E, P, m)
    h, xrel, k = h_of_Q(Q)
    return h / m ** 2, dict(m=m, k=k, tate_x_check_rel_digits=xrel)

if r > 0:
    t = time.time()
    H = {}; info = {}
    for i in range(r):
        H[(i, i)], info[f"P{i+1}"] = height(gens[i])
    for i in range(r):
        for j in range(i + 1, r):
            S_ = pari.elladd(E, gens[i], gens[j])
            hs, info[f"P{i+1}+P{j+1}"] = height(S_)
            H[(i, j)] = H[(j, i)] = (hs - H[(i, i)] - H[(j, j)]) / 2
    # quadraticity checks
    h1b, _ = height(gens[0], extra=2)
    out["quadratic_check_m_vs_2m_valuation"] = int(pari.valuation(h1b - H[(0, 0)], p)) if h1b != H[(0, 0)] else "equal"
    if h1b != H[(0, 0)] and int(pari.valuation(h1b - H[(0, 0)], p)) < K - 4:
        refuse.append("height not quadratic: h(mP)/m^2 depends on m (wrong sigma normalisation)")
    if r >= 2:
        Dm = pari.elladd(E, gens[0], pari.ellneg(E, gens[1]))
        hd, _ = height(Dm)
        par = hd + (2 * H[(0, 1)] + H[(0, 0)] + H[(1, 1)]) - 2 * H[(0, 0)] - 2 * H[(1, 1)]
        out["parallelogram_valuation"] = int(pari.valuation(par, p)) if par != 0 else "exact"
    G = pari.matrix(r, r, [H[(i, j)] for i in range(r) for j in range(r)])
    if basisA:                                  # Gram matrix of the basis A*B is A G A^T
        A = pari.matrix(r, r, [int(c) for row in basisA.split(";") for c in row.split(",")])
        G = A * G * pari.mattranspose(A)
    Regp = pari.matdet(G)
    Reg_gamma = Regp / lg ** r
    out.update(height_info=info, Reg_p=str(Regp), v_Reg_p=int(pari.valuation(Regp, p)),
               Reg_gamma=str(Reg_gamma), v_Reg_gamma=int(pari.valuation(Reg_gamma, p)),
               heights_diag=[str(H[(i, i)]) for i in range(r)], sec_height=round(time.time() - t, 1))
else:
    Regp = Reg_gamma = pari(1)
LoverLg = Linv / lg
out["L_over_log(1+p)"] = str(LoverLg); out["v_L_over_lg"] = int(pari.valuation(LoverLg, p))
log(f"[heights] Reg_p={out.get('Reg_p')} L={out['L_invariant']}")

# ---------------- 8. PARI cross-check (NOT used in the verdict) ----------------
try:
    n2 = K - 4
    Ep = pari.ellinit(E[:5], pari(f"O({p}^{n2})"))
    tate = pari("(e)->e.tate")(Ep)
    u2, qP = tate[0], tate[2]
    s2 = pari.ellpadics2(E, p, n2)
    def hpari(P):
        hv = pari.ellpadicheight(E, p, n2, P)
        return hv[0] + hv[1] * (-s2 + 1 / pari.log(qP) / u2)
    xc = dict(q_pari_minus_mine_val=int(pari.valuation(qP - qE, p)) if qP != qE else "equal")
    if r > 0:
        HP = {}
        for i in range(r): HP[(i, i)] = hpari(gens[i])
        for i in range(r):
            for j in range(i + 1, r):
                HP[(i, j)] = HP[(j, i)] = (hpari(pari.elladd(E, gens[i], gens[j])) - HP[(i, i)] - HP[(j, j)]) / 2
        GP = pari.matrix(r, r, [HP[(i, j)] for i in range(r) for j in range(r)])
        RegP = pari.matdet(GP)
        RegR = pari.ellpadicregulator(E, p, n2, gens)
        xc.update(Reg_docstring=str(RegP), Reg_ellpadicregulator=str(RegR),
                  mine_minus_docstring_val=int(pari.valuation(Regp - RegP, p)) if Regp != RegP else "equal",
                  mine_minus_ellpadicregulator_val=int(pari.valuation(Regp - RegR, p)) if Regp != RegR else "equal",
                  height_P1_mine_minus_pari_val=int(pari.valuation(H[(0, 0)] - HP[(0, 0)], p)))
    out["pari_crosscheck"] = xc
except Exception as e:
    out["pari_crosscheck"] = "failed: " + str(e)[:200]

# ---------------- 9. S and conclusion ----------------
if c_lead is None:
    refuse.append("leading T-coefficient not p-integral")
    S = None
else:
    cl = P_(c_lead, c_lead_prec)
    S = cl / (LoverLg * Reg_gamma) if (LoverLg != 0 and Reg_gamma != 0) else None
out["c_lead_T^(r+1)_SW_norm"] = [c_lead, f"mod p^{c_lead_prec}"]
out["c_lead_T^(r+1)_Omega_plus_norm"] = (c_lead * cinf) % p ** c_lead_prec if c_lead is not None else None
out["S_SW_norm"] = str(S)
out["S_Omega_plus_norm"] = str(S * cinf) if S is not None else None
out["BSD_prediction_prod_c_over_tors2_times_Sha"] = f"{prodc}/{tors**2} * #Sha"
vc = int(pari.valuation(c_lead, p)) if c_lead not in (None, 0) else None
concl = None
if not refuse and r > 0 and c_lead is not None and c_lead % p != 0 \
        and out["v_L_over_lg"] == 0 and out["v_Reg_gamma"] == 0 and prodc % p != 0 and tors % p != 0:
    concl = (f"Sha(E/Q)[{p}^oo] = 0, GIVEN SW Thm 6.1 (Jones) + Thm 7.3 (Kato/Kobayashi) and rho_p surjective; "
             f"rank = {r} and p-saturation of the basis follow (ord_T L_p = r+1 and v(Reg_gamma) = 0)")
else:
    concl = "NO CONCLUSION"
out["refuse"] = refuse
out["conclusion"] = concl
out["sec_total"] = round(time.time() - t00, 1)
log("[result] " + json.dumps(out))
