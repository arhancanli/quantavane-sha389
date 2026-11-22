# Sha(389a1)[p^∞] = 0 for every prime p < 14,633

**Status: a Tier 5 verification (computational; the lowest step of our five-step scale). It does not show that Sha is
finite, and it is not yet peer reviewed. Verified by two independently written programs and an internal referee.**

Author: Arhan Canli (Quantavane, arhan@quantavane.com). On the use of AI, see
[Authorship and the role of AI](#authorship-and-the-role-of-ai).

## The statement

Let E be the elliptic curve 389a1, y² + y = x³ + x² − 2x, of conductor 389 and rank 2, with E(Q) = Z(−1, 1) ⊕ Z(0, 0).
The primes p < 14,633 at which E has supersingular reduction are 107, 599, 1049, 2957, 6661, 8263, 9397 and 9551; E has split multiplicative
reduction at its only bad prime, 389.

**Statement.** Sha(E/Q)[p^∞] = 0 for p = 107, 599, 1049, 2957, 6661, 8263, 9397, 9551 and 389.

Together with Stein and Wuthrich's Theorem 12.3 (Math. Comp. 82 (2013), 1757–1792: Sha(E/Q)[p] = 0 for p = 2 and every
good ordinary prime p < 48,859 except 16,231), this gives:

**Sha(389a1)[p^∞] = 0 for every prime p < 14,633.**

Stein and Wuthrich left out exactly the supersingular primes (their Remark 12.5: "not because our algorithms do not
apply (they do apply), but because our implementations are significantly slower in this case"), the bad prime 389 and
the ordinary prime 16,231. Stein's 2013 JMM talk lists "generalize elliptic curve computations to supersingular
primes" and "verify that Sha(389a)[16,231] = 0" as future work.

## What this is and what it is not

- It is a finite computation that, combined with published theorems (listed below), proves the statement for eight
  primes. The method is Stein and Wuthrich's Algorithm 11.1 in the supersingular case, fed by Pollack's ± p-adic
  L-functions computed from Mazur–Tate elements, and at p = 389 their split-multiplicative case (exceptional zero,
  L-invariant, Schneider's height). Nothing in the method is new.
- It is not a proof that Sha(389a1) is finite. That is open: no elliptic curve over Q of rank at least 2 is known to
  have a finite Tate–Shafarevich group. Each prime is handled separately.
- It assumes nothing conjectural, but it relies on theorems we cite without re-reading every proof (marked below).
- It has not been reviewed by anyone outside Quantavane.

## Results

S is the ratio of the leading coefficient of (1 − φ)⁻² L_p(E, T) to the normalised p-adic regulator Reg_γ, reduced
mod p (first coordinate in the basis (ω, Fω)). If S is a p-adic unit, Sha(E/Q)(p) = 0. Program A normalises the
modular symbols by Ω⁺ (so it expects S = c_∞ · #Sha = 2); program B normalises by Ω_E = 2Ω⁺ (so it expects S = 1).

| p | c⁺ mod p (T⁰..T³), program A | Reg_γ mod p | S, program A | S, program B | second coordinate |
|---|---|---|---|---|---|
| 107 | 0, 0, 21, 60 | 43, 59 | 2 | 1 | consistent (c⁻ from level p³) |
| 599 | 0, 0, 56, 383 | 571, 416 | 2 | 1 | consistent (c⁻ from level p³) |
| 1049 | 0, 0, 847, 972 | 101, 761 | 2 | 1 | not computed |
| 2957 | 0, 0, 2547, 2802 | 205, 698 | 2 | 1 | not computed |
| 6661 | 0, 0, 3666, 3665 | 4828, 705 | 2 | 1 | not computed |
| 8263 | 0, 0, 3529, 6834 | 2367, 640 | 2 | 1 | not computed |
| 9397 | 0, 0, 6604, 2356 | 6095, 204 | 2 | 1 | not computed |
| 9551 | 0, 0, 2258, 968 | 8422, 3323 | 2 | 1 | not computed |

In every case S equals what the p-adic Birch–Swinnerton-Dyer formula predicts for #Sha = 1, which is a stronger check
than "S is a unit". Logs: `logs/supersingular/`.

**p = 389 (split multiplicative).** L_p has an exceptional zero, so the leading coefficient is the T³ coefficient:
[T⁰..T⁴] L_p ≡ 0, 0, 0, 321, 102 (mod 389) in the Ω⁺ normalisation (program A; program B gives 355, 51 = 321/2, 102/2
with Ω_E). The L-invariant is L = log_p(q_E)/ord_p(q_E) = 191·389 + 370·389² + …, Reg_γ = 232 + 3·389 + …, and
S = [T³]L_p / ((L/log_p(1 + p)) · Reg_γ) = 2 (program A, Ω⁺) = 1 (program B, Ω_E): again the BSD value for #Sha = 1.
Logs: `logs/split389/`.

## How it works (details in [ROUTE.md](ROUTE.md))

1. Kato's divisibility (in the supersingular case: Kato 2004, Thm 12.5; Kobayashi 2003, Thm 6.3), as stated in Stein–
   Wuthrich Thm 7.3: L_p(E, T) = f_E(T) d(T) with d in Λ, for E semistable at p and ρ_p surjective.
2. Perrin-Riou (Stein–Wuthrich Thm 6.2): ord f_E ≥ r, with equality exactly when the p-adic regulator is non-zero and
   Sha(E/Q)(p) is finite; the leading term then has the valuation of ∏c_v · #Sha(p) · Reg_p.
3. Pollack (Duke Math. J. 118, 2003): for a_p = 0, L_α = log⁺ L⁺ + α log⁻ L⁻ with L± in Λ. At level p² and p³ the
   Mazur–Tate elements give c±_j mod p (exact identities, ROUTE.md section 3), so the T² coefficient of L_p mod p
   needs only p² (or p³) modular symbol values and under 200 MB of memory. (PARI's overconvergent `ellpadicbsd` needs
   about 3.1 MB · p², i.e. about 35 GB at p = 107: `logs/pari_ellpadicbsd_memory_389a1_20261122_0446.log`.)
4. The regulator: PARI's `ellpadicregulator` on the basis (−1, 1), (0, 0), which is saturated at every prime up to
   50,000 (`logs/saturation_50000_*.log`), divided by log_p(1 + p)².
5. At p = 389 (split multiplicative): Kato's divisibility with the extra zero of the Coleman map (Kobayashi 2006), the
   leading-term theorem of Schneider, Perrin-Riou and Jones with the L-invariant, and the T³ coefficient of L_p from the
   level-389² Mazur–Tate element (ROUTE.md section 6).

## The programs

Supersingular primes:

- `programs/pm_route.py` (program A, run with PARI 2.15.4): symbols from `msfromell`, the Mazur–Tate projection by a
  table of (1 + p)^s and Teichmüller representatives, binomial coefficients by Lucas' theorem, the closed form for
  (1 − φ)⁻², and `ellpadicregulator`.
- `programs/ref_m2ss.py` (program B, written separately by the internal referee, run with PARI 2.17.2): its own plus
  eigensymbol (`msinit` + `mssplit`) scaled by twisted L-values, the projection through a primitive root and the p-adic
  logarithm, Horner's rule in F_p[T]/(T^D), (1 − φ)⁻² as a p-adic matrix in the basis (ω, η), and the regulator
  assembled from `ellpadicheight` (checked against `ellpadicregulator`). It computes S for three functionals, which
  Perrin-Riou's proportionality forces to agree.

p = 389:

- `programs/split_mult.py` (program A): Mazur–Tate element of level p² from `msfromell`, the L-invariant from PARI's
  Tate parametrisation, and the split-multiplicative height in PARI's form ellpadicheight · [1, −s₂ + 1/(log q · u²)]
  (checked against `ellpadicregulator`).
- `programs/ref_split.py` (program B, written separately by the internal referee): its own eigensymbol scaled by
  twisted L-values (the scale checked on nine discriminants, integrality checked on all generators), q_E by inverting
  j, and the height through Tate's parametrisation with its own formal logarithm and σ-product; PARI's heights only as a
  cross-check.

## Controls

- Stein–Wuthrich section 12.2 (1483a1, p = 5, supersingular, rank 2): PARI reproduces their printed Reg_γ digit for
  digit (`checks/control_1483a1.py`); both programs give S = 1 there.
- Rank-0 interpolation (11a1, p = 19): the programs' c±₀ match L(E, 1)/Ω⁺ = 1/5 exactly mod 19.
- c_∞ = 2 (the normalisation that 389a1 needs): two rank-2 curves with Δ > 0 and a supersingular prime 5
  (conductors 1613 and 1627) where program B's leading terms equal PARI's overconvergent `ellpadicbsd` mod 5
  (`checks/ctrl_cinf2_search.py`, `logs/controls/`).
- Red controls (each must be refused): a basis that is not saturated at p ([p·P₁, P₂]) loses two powers of p in the
  regulator; flipping the sign of c⁻ breaks the second-coordinate check; symbols scaled by p are refused.
- Split multiplicative: Stein–Wuthrich section 12.1 (446d1, p = 223): PARI's split height reproduces their printed
  Reg_p = 153·223² + 125·223³ + 124·223⁴ + … and L = 179·223 + 85·223² + 30·223³ + … digit for digit, and both programs
  give their printed T³ coefficient 139 (program A: 55 = 2·139 with Ω⁺) and S = ∏c_v. 91b1 at p = 7 (rank 1): program A's
  [T²] L_p agrees with PARI's overconvergent `ellpadicbsd` to O(7⁴), and S = 1/9 = 1/#E(Q)²_tors mod 7³. The referee's
  program adds conductors 185 (p = 37) and 92 (p = 23). Red controls: mis-scaled symbols, a non-saturated basis, a
  nonsplit prime, a good prime and a height without Stein–Wuthrich's constant C all give no conclusion.

## Running the checks

Needs Python 3 with [cypari2](https://github.com/sagemath/cypari2) (PARI/GP 2.15 or later).

    PY=python3 zsh run_checks.sh            # p = 107 and p = 389 with both programs, every control (~1 minute)
    FULL=1 PY=python3 zsh run_checks.sh     # adds p = 599, 1049, 2957 (~10 more minutes)
    FULL=2 PY=python3 zsh run_checks.sh     # as FULL=1, plus 6661, 8263, 9397, 9551 (about 2 more hours)

Each run writes its logs to `out/` and ends with `checks/verify.py`, which prints one PASS/FAIL line per check.

## Credits

The method and Theorem 12.3 are W. Stein and C. Wuthrich's (*Algorithms for the arithmetic of elliptic curves using
Iwasawa theory*, Math. Comp. 82 (2013)). The ± p-adic L-functions are R. Pollack's (*On the p-adic L-function of a
modular form at a supersingular prime*, Duke Math. J. 118 (2003)); the integrality statement we use is as restated by
M. Kurihara and R. Pollack. The leading-term theorem and the Dieudonné-module-valued p-adic BSD formula are
B. Perrin-Riou's (*Arithmétique des courbes elliptiques à réduction supersingulière en p*, Experiment. Math. 12
(2003); with D. Bernardi, C. R. Acad. Sci. Paris 317 (1993), for the conjecture). Mazur–Tate elements and the p-adic L-function come from B. Mazur and J. Tate
(*Refined conjectures of the "Birch and Swinnerton-Dyer type"*, Duke Math. J. 54 (1987)) and Mazur–Tate–Teitelbaum
(Invent. Math. 84 (1986)). The divisibility is K. Kato's (Astérisque 295 (2004)) and S. Kobayashi's (Invent. Math. 152
(2003)); in the split multiplicative case the extra zero of the Coleman map is in S. Kobayashi's *An elementary proof
of the Mazur–Tate–Teitelbaum conjecture for elliptic curves* (Doc. Math. Extra Vol. Coates, 2006). The height there is
P. Schneider's (Invent. Math. 69 (1982)) with A. Werner's correction (Doc. Math. 3 (1998)), and the leading-term theorem
is due to Schneider, Perrin-Riou and J. W. Jones (Duke Math. J. 59 (1989)). The curve data follow J. Cremona's tables. Every computation uses PARI/GP (`msfromell`, `mseval`, `msinit`,
`ellpadicregulator`, `ellpadicheight`, `ellpadicbsd`, `ellsaturation`, `ellrank`) by the PARI Group.

## Not yet done

- Outside review: none yet. W. Stein and C. Wuthrich were told of this check before this repository was made public.
- More primes: the nine supersingular primes between 14,633 and 45,161 are being run with both programs; the ordinary
  prime 16,231 (where ord_p Reg_p = 3) needs another route (overconvergent modular
  symbols, as Stein and Wuthrich suggest in their Remark 12.4).
- Lean: this is out of Lean's reach for now (the inputs are deep theorems in Iwasawa theory, not finite statements).

## Authorship and the role of AI

Author and maintainer: Arhan Canli (Quantavane, https://quantavane.com), who takes responsibility for the content.

The programs, the computations, the literature search and the text were produced by AI agents (Claude, by Anthropic)
working for Quantavane under Arhan Canli's direction; an internal referee, also an AI agent, wrote program B
separately and checked the argument and the citations. No mathematician outside Quantavane has reviewed this
repository. Every number here can be re-checked with `run_checks.sh`.

## Licence and citation

The text is under CC BY 4.0 (`LICENSE`); the programs are under the MIT licence (`LICENSE-CODE`). Please cite as:
Arhan Canli, *Sha(389a1)[p^∞] = 0 for every prime p < 14,633*,
https://github.com/arhancanli/quantavane-sha389 (2026), together with Stein and Wuthrich (Math. Comp. 2013).
