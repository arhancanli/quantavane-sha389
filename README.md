# Sha(389a1)[p^∞] = 0 at the supersingular primes below 6661

**Status: a Tier 5 verification (computational; the lowest step of our five-step scale). It does not show that Sha is
finite, and it is not yet peer reviewed. Verified by two independently written programs and an internal referee.**

Author: Arhan Canli (Quantavane, arhan@quantavane.com). On the use of AI, see
[Authorship and the role of AI](#authorship-and-the-role-of-ai).

## The statement

Let E be the elliptic curve 389a1, y² + y = x³ + x² − 2x, of conductor 389 and rank 2, with E(Q) = Z(−1, 1) ⊕ Z(0, 0).
The primes p < 6661 at which E has supersingular reduction are 107, 599, 1049 and 2957.

**Statement.** Sha(E/Q)[p^∞] = 0 for p = 107, 599, 1049 and 2957.

Together with Stein and Wuthrich's Theorem 12.3 (Math. Comp. 82 (2013), 1757–1792: Sha(E/Q)[p] = 0 for p = 2 and every
good ordinary prime p < 48,859 except 16,231), this gives:

**Sha(389a1)[p^∞] = 0 for every prime p < 6661 except p = 389.**

Stein and Wuthrich left out exactly the supersingular primes (their Remark 12.5: "not because our algorithms do not
apply (they do apply), but because our implementations are significantly slower in this case"), the bad prime 389 and
the ordinary prime 16,231.

## What this is and what it is not

- It is a finite computation that, combined with published theorems (listed below), proves the statement for four
  primes. The method is Stein and Wuthrich's Algorithm 11.1 in the supersingular case, fed by Pollack's ± p-adic
  L-functions computed from Mazur–Tate elements. Nothing in the method is new.
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

In every case S equals what the p-adic Birch–Swinnerton-Dyer formula predicts for #Sha = 1, which is a stronger check
than "S is a unit". Logs: `logs/supersingular/`.

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

## The two programs

- `programs/pm_route.py` (program A, run with PARI 2.15.4): symbols from `msfromell`, the Mazur–Tate projection by a
  table of (1 + p)^s and Teichmüller representatives, binomial coefficients by Lucas' theorem, the closed form for
  (1 − φ)⁻², and `ellpadicregulator`.
- `programs/ref_m2ss.py` (program B, written separately by the internal referee, run with PARI 2.17.2): its own plus
  eigensymbol (`msinit` + `mssplit`) scaled by twisted L-values, the projection through a primitive root and the p-adic
  logarithm, Horner's rule in F_p[T]/(T^D), (1 − φ)⁻² as a p-adic matrix in the basis (ω, η), and the regulator
  assembled from `ellpadicheight` (checked against `ellpadicregulator`). It computes S for three functionals, which
  Perrin-Riou's proportionality forces to agree.

## Controls

- Stein–Wuthrich section 12.2 (1483a1, p = 5, supersingular, rank 2): PARI reproduces their printed Reg_γ digit for
  digit (`checks/control_1483a1.py`); both programs give S = 1 there.
- Rank-0 interpolation (11a1, p = 19): the programs' c±₀ match L(E, 1)/Ω⁺ = 1/5 exactly mod 19.
- c_∞ = 2 (the normalisation that 389a1 needs): two rank-2 curves with Δ > 0 and a supersingular prime 5
  (conductors 1613 and 1627) where program B's leading terms equal PARI's overconvergent `ellpadicbsd` mod 5
  (`checks/ctrl_cinf2_search.py`, `logs/controls/`).
- Red controls (each must be refused): a basis that is not saturated at p ([p·P₁, P₂]) loses two powers of p in the
  regulator; flipping the sign of c⁻ breaks the second-coordinate check; symbols scaled by p are refused.

## Running the checks

Needs Python 3 with [cypari2](https://github.com/sagemath/cypari2) (PARI/GP 2.15 or later).

    PY=python3 zsh run_checks.sh            # p = 107 with both programs, every control and red control (~1 minute)
    FULL=1 PY=python3 zsh run_checks.sh     # adds p = 599, 1049, 2957 (~10 more minutes)

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
(2003)). The curve data follow J. Cremona's tables. Every computation uses PARI/GP (`msfromell`, `mseval`, `msinit`,
`ellpadicregulator`, `ellpadicheight`, `ellpadicbsd`, `ellsaturation`, `ellrank`) by the PARI Group.

## Not yet done

- Outside review: none yet. W. Stein and C. Wuthrich were told of this check before this repository was made public.
- More primes: the supersingular primes 6661, 8263, 9397, 9551 and the nine between 14,633 and 45,161 are being run with
  both programs; the bad prime 389 (split multiplicative, an exceptional zero) and 16,231 need other routes.
- Lean: this is out of Lean's reach for now (the inputs are deep theorems in Iwasawa theory, not finite statements).

## Authorship and the role of AI

Author and maintainer: Arhan Canli (Quantavane, https://quantavane.com), who takes responsibility for the content.

The programs, the computations, the literature search and the text were produced by AI agents (Claude, by Anthropic)
working for Quantavane under Arhan Canli's direction; an internal referee, also an AI agent, wrote program B
separately and checked the argument and the citations. No mathematician outside Quantavane has reviewed this
repository. Every number here can be re-checked with `run_checks.sh`.

## Licence and citation

The text is under CC BY 4.0 (`LICENSE`); the programs are under the MIT licence (`LICENSE-CODE`). Please cite as:
Arhan Canli, *Sha(389a1)[p^∞] = 0 at the supersingular primes below 6661*,
https://github.com/arhancanli/quantavane-sha389 (2026), together with Stein and Wuthrich (Math. Comp. 2013).
