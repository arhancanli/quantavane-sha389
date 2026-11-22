"""Second-algorithm check of the main number: PARI's ellpadicL / ellpadicbsd (mspadic moments, not level-p^n
Mazur-Tate elements) for 389a1 at p = 389.  Hard time limit via SIGALRM.
L_p^{(3)}(chi^0) = int log_p^3 dmu = 3! * c3 * lg^3 when c0 = c1 = c2 = 0 (T = (1+p)^s - 1).
Usage: python -I ctrl_ellpadicL_389.py n [limit_seconds]"""
import sys, json, time, signal
import cypari2
pari = cypari2.Pari(size=8 * 10**8, sizemax=1_800_000_000)
n = int(sys.argv[1]); lim = int(sys.argv[2]) if len(sys.argv) > 2 else 600
signal.signal(signal.SIGALRM, lambda *a: (print(json.dumps({"timeout_after_s": lim}), flush=True), sys.exit(3)))
signal.alarm(lim)
p = 389
E = pari.ellinit([0, 1, 1, -2, 0])
lg = pari(f"log(1+{p}+O({p}^{n+6}))")
out = dict(p=p, n=n)
t = time.time()
L3 = pari.ellpadicL(E, p, n, 0, 3)
out["ellpadicL_r3"] = str(L3); out["sec_L3"] = round(time.time() - t, 1)
c3_mine_SW = pari(355) + pari(f"O({p})")          # ref_split.py level-2 result, SW Omega_E normalisation
out["ratio_ellpadicL3_over_6_c3SW_lg3"] = str(L3 / (6 * c3_mine_SW * lg ** 3))
for k in (0, 1, 2):
    out[f"ellpadicL_r{k}"] = str(pari.ellpadicL(E, p, n, 0, k))
print(json.dumps(out), flush=True)
t = time.time()
rk, Lp = pari.ellpadicbsd(E, p, n)
out.update(ellpadicbsd_r=int(rk), ellpadicbsd_Lp=str(Lp), sec_bsd=round(time.time() - t, 1))
R = pari.ellpadicregulator(E, p, n + 2, [[-1, 1], [0, 0]])
out["ellpadicbsd_Lp_over_ellpadicregulator"] = str(Lp / R)
print(json.dumps(out), flush=True)
