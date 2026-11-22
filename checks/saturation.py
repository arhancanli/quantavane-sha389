"""Saturation of the basis P1 = (-1,1), P2 = (0,0) of 389a1 at every prime up to a bound (PARI ellsaturation), plus
rank bounds (ellrank), torsion, Tamagawa numbers and the reduction type at 389. A regulator computed on a basis that is
not saturated at p would carry extra powers of p (see the red controls in the programs).
Usage: python -I saturation.py [bound=15000]"""
import sys, json, time
import cypari2
pari = cypari2.Pari(size=2 * 10**8, sizemax=1_200_000_000)
B = int(sys.argv[1]) if len(sys.argv) > 1 else 15000
E = pari.ellinit([0, 1, 1, -2, 0])
P = [[-1, 1], [0, 0]]
t0 = time.time()
S = pari.ellsaturation(E, P, B)
idx2 = pari.matdet(pari.ellheightmatrix(E, P)) / pari.matdet(pari.ellheightmatrix(E, S))
gr = pari.ellglobalred(E)
print(json.dumps(dict(bound=B, saturated_basis=str(S), index_squared=str(pari.round(idx2)), rank_bounds=str(pari.ellrank(E)[:2]),
                      torsion=int(pari.elltors(E)[0]), conductor=int(gr[0]), tamagawa_product=int(gr[2]),
                      a_389=int(pari.ellap(E, 389)), localred_389=str(pari.elllocalred(E, 389)),
                      disc=int(E.disc()), sec=round(time.time() - t0, 1))))
