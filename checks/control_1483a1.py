"""M2-ss pipeline control (22 Nov 2026): reproduce Stein-Wuthrich, Math. Comp. 2013, section 12.2 (curve 1483a1, rank 2,
generators (-1,0), (0,1), supersingular at p = 5) with PARI ellpadicbsd / ellpadicregulator, and fix the basis convention.
SW print, basis (omega, phi(omega)) with phi = F/p:
  L_p T^2 coefficient = (1 + O(5)) omega + (4*5 + O(5^2)) phi(omega)
  Reg_gamma = (1 + 2*5 + 3*5^2 + 5^3 + O(5^5)) omega + (4*5 + 4*5^2 + 4*5^3 + 5^4 + 2*5^5 + O(5^6)) phi(omega)
  and (1-phi)^-2 L_p^* = #Sha * Reg with Sha(5) trivial.
PARI returns D-valued quantities in the basis (omega, F omega), F = p*phi, so its second coordinate should be the SW one
divided by p (if the normalisations otherwise agree). SW's Reg_gamma is normalised by log_p(gamma)^r; PARI's is not
(checked below by printing both). Red control: a non-saturated basis [5*P1, P2] multiplies the regulator by 25."""
import sys, json
import cypari2
pari = cypari2.Pari(size=10**8, sizemax=1_500_000_000)
p, n = 5, int(sys.argv[1]) if len(sys.argv) > 1 else 6
E = pari.ellinit([0, 1, 1, 2, 2])          # 1483a1: y^2 + y = x^3 + x^2 + 2x + 2 (SW print reads "y2 + x", a text-extraction slip; this model has conductor 1483 and both printed generators)
G = [[-1, 0], [0, 1]]
print("conductor", pari.ellglobalred(E)[0], "on curve", [int(pari.ellisoncurve(E, P)) for P in G], "a_5", pari.ellap(E, 5))
r, Lp = pari.ellpadicbsd(E, p, n)
R = pari.ellpadicregulator(E, p, n, G)
Rbad = pari.ellpadicregulator(E, p, n, [pari.ellmul(E, G[0], p), G[1]])
print("r =", r); print("Lp (omega, F omega) =", Lp); print("Reg (omega, F omega) =", R)
v = lambda x: int(pari.valuation(x, p))
cross = Lp[0] * R[1] - Lp[1] * R[0]
S0 = Lp[0] / R[0]
print("proportionality Lp[0]*R[1]-Lp[1]*R[0] =", cross, " (valuation", v(cross) if cross != 0 else "exact 0", ")")
print("S = Lp[0]/R[0] =", S0, " v_p(S) =", v(S0))
print("red control: Reg with [5P1,P2] =", Rbad, " v_p(Lp[0]/Rbad[0]) =", v(Lp[0] / Rbad[0]))
# SW coordinates: phi(omega) = F(omega)/p, so x*F(omega) = (p*x)*phi(omega)
print("PARI second coordinates times p (= SW phi(omega) coordinates): Lp", Lp[1] * p, " Reg", R[1] * p)
log_gamma = pari.log(pari('1 + 5 + O(5^12)'))
print("log_p(1+p) =", log_gamma, "; Reg/log^2 (SW Reg_gamma if PARI unnormalised):", R[0] / log_gamma**2, R[1] * p / log_gamma**2)
