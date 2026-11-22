"""Checks the outputs written by run_checks.sh against the expected values. Prints one PASS/FAIL line per check and
exits non-zero if any check fails. Usage: python -I checks/verify.py OUTDIR"""
import sys, json, os
D = sys.argv[1]
fails = 0


def load(name):
    path = os.path.join(D, name + ".log")
    if not os.path.exists(path):
        return None
    rows = [json.loads(l) for l in open(path) if l.lstrip().startswith(("{", "["))]
    return rows[-1] if rows else open(path).read()


def check(label, ok):
    global fails
    print(("PASS  " if ok else "FAIL  ") + label)
    fails += 0 if ok else 1


txt = load("primes")
check("supersingular primes of 389a1 below 9600 are 107, 599, 1049, 2957, 6661, 8263, 9397, 9551",
      isinstance(txt, str) and "[107, 599, 1049, 2957, 6661, 8263, 9397, 9551]" in txt)
s = load("saturation")
check("basis (-1,1), (0,0) saturated at every prime <= 12000; rank 2; torsion 1; Tamagawa 1; a_389 = 1",
      s["index_squared"] == "1" and s["rank_bounds"] == "[2, 2]" and s["torsion"] == 1 and s["tamagawa_product"] == 1
      and s["a_389"] == 1)
t = load("sw122_pari")
check("SW 2013 section 12.2 (1483a1, p = 5): PARI's Reg_gamma = SW's printed (1 + 2*5 + 3*5^2 + 5^3 + ...) omega + "
      "(4*5 + 4*5^2 + 4*5^3 + 5^4 + 2*5^5 + ...) phi(omega), and S = 1 + O(5^5)",
      isinstance(t, str) and "1 + 2*5 + 3*5^2 + 5^3 + " in t and "4*5 + 4*5^2 + 4*5^3 + 5^4 + 2*5^5 + " in t
      and "S = Lp[0]/R[0] = 1 + O(5^5)" in t)
a = load("A_1483a1_p5")
check("program A, 1483a1 p = 5: S = 1, c+ = [0,0,4,4], c- = [0,0,1,4]",
      a["S_mod_p"] == 1 and a["cplus_mod_p"] == [0, 0, 4, 4] and a["cminus_mod_p"] == [0, 0, 1, 4])
a = load("A_11a1_p19")
check("program A, 11a1 p = 19 (rank 0 interpolation, L/Omega+ = 1/5): c+_0 = 15, c-_0 = 8",
      a["cplus_mod_p"][0] == 15 and a["cminus_mod_p"][0] == 8)
a = load("A_389a1_p107")
check("program A, 389a1 p = 107: S = 2 (= c_oo #Sha with #Sha = 1), second coordinate consistent",
      a["S_mod_p"] == 2 and a["consistency_second_coord"] is True)
check("program A red control R1: basis [107 P1, P2] gives v(Reg_gamma) = (2, 2), no conclusion",
      a["red_R1_v_Reg_nonsaturated"] == [2, 2])
check("program A red control R2: sign-flipped c- is inconsistent", a["red_R2_flipped_sign_consistent"] is False)
b = load("B_1483a1_p5")
check("program B, 1483a1 p = 5: PASS with S = 1", b["S_modp"] == 1 and b["verdict"].startswith("PASS"))
b = load("B_11a1_p19")
check("program B, 11a1 p = 19: PASS (S = L/Omega_E = 1/5 = 4 mod 19)", b["S_modp"] == 4 and b["verdict"].startswith("PASS"))
b = load("B_389a1_p107")
check("program B, 389a1 p = 107: PASS with S = 1 (Omega_E normalisation)", b["S_modp"] == 1 and b["verdict"].startswith("PASS"))
b = load("B_red_p107_scaled_by_p")
check("program B red control: symbols scaled by 107 are refused", b["verdict"].startswith("REFUSE"))
for p in (599, 1049, 2957):
    a, b = load(f"A_389a1_p{p}"), load(f"B_389a1_p{p}")
    if a is None and b is None:
        continue
    check(f"p = {p}: program A S = 2 and program B PASS with S = 1",
          a["S_mod_p"] == 2 and b["S_modp"] == 1 and b["verdict"].startswith("PASS"))
print(f"{fails} check(s) failed" if fails else "all checks passed")
sys.exit(1 if fails else 0)
