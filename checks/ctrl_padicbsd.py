"""Independent-algorithm control: compare ref_split.py's leading coefficient with PARI's ellpadicbsd
(overconvergent / mspadic route, a different algorithm from Mazur-Tate elements at level p^n).
Usage: python -I ctrl_padicbsd.py <ref_split log file> n
Reads [result] JSON from the log (curve, p, c_lead mod p^e, L, prod_c, torsion, c_inf), then:
  Lp_pari = ellpadicbsd(E,p,n)[2]   (doc: conjecturally R_p * #Sha, R_p = ellpadicregulator)
  mine    = c_lead * lg^(r+1) * tors^2 / (L * prod_c)  [SW Omega_E normalisation], lg = log_p(1+p)
and prints Lp_pari / mine and both Sha estimates.
"""
import sys, json
import cypari2
pari = cypari2.Pari(size=4 * 10**8, sizemax=1_500_000_000)
d = None
for line in open(sys.argv[1]):
    if line.startswith("[result]"): d = json.loads(line[9:])
n = int(sys.argv[2])
p = d["p"]; r = d["r_points"]
E = pari.ellinit(d["curve"])
c, e = d["c_lead_T^(r+1)_SW_norm"][0], int(d["c_lead_T^(r+1)_SW_norm"][1].split("^")[1])
Lp_mine_c = pari(c) + pari(f"O({p}^{e})")
L = pari(d["L_invariant"])
lg = pari(f"log(1+{p}+O({p}^{n+4}))")
mine = Lp_mine_c * lg ** (r + 1) * d["torsion"] ** 2 / (L * d["prod_c"])
rk, Lp = pari.ellpadicbsd(E, p, n)
gens = None
res = dict(curve=d["curve"], p=p, n=n, ellpadicbsd_r=int(rk), ellpadicbsd_Lp=str(Lp), mine_R_times_Sha=str(mine),
           ratio_pari_over_mine=str(Lp / mine))
regs = pari(d["Reg_p"])
res["Sha_estimate_pari"] = str(Lp / regs)
res["Sha_estimate_mine"] = str(mine / regs)
print(json.dumps(res))
