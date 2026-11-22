# The route, written out

E = 389a1: y² + y = x³ + x² − 2x, N = 389 (prime, so E is semistable), E(Q) = Z(−1, 1) ⊕ Z(0, 0), Tamagawa product 1,
trivial torsion, c_∞ = 2 (Δ = 389 > 0). Below, p is one of 107, 599, 1049, 2957, where a_p = 0.

## 1. Inputs (each checked against print unless marked)

I1. Kato's divisibility, as stated in Stein–Wuthrich (SW) Thm 7.3: "Suppose that E has semistable reduction at p and
    that ρ_p is surjective. Then there exists a series d(T) in Λ such that L_p(E, T) = f_E(T) d(T)." In the
    supersingular case L_p and f_E are D_p(E)-valued (SW 3.5 and 6.2). For supersingular p the divisibility is Kato
    2004, Thm 12.5 with the explicit reciprocity law, or Kobayashi 2003, Thm 6.3 for L± (as cited by Kurihara–Pollack,
    Thm 1.3); neither was re-read. For p ≥ 5, ρ_p is surjective iff ρ̄_p is (SW Prop. 7.2).

I2. Perrin-Riou (SW Thm 6.2): ord_{T=0} f_E ≥ r = rank E(Q), with equality iff Reg_p(E/Q) ≠ 0 and Sha(E/Q)(p) is finite;
    then the leading term of (1 − φ)⁻² f_E has the valuation of ∏c_v · #Sha(p) · Reg_p. Proportionality in D_p
    (Perrin-Riou 2003: "quotient de deux vecteurs liés") makes the first-coordinate ratio legitimate.

I3. SW Algorithm 11.1, supersingular case: b_p = ord_p((1 − φ)⁻² L_p*) − ord_p(Reg_γ) − Σ ord_p(c_v), and
    #Sha(E/Q)(p) ≤ p^(b_p), for E semistable at p, p > 2, ρ̄_p surjective and E(Q) known. Reg_γ = Reg_p / log_p(κ(γ))^r
    is SW's normalised regulator (SW 4.4, (5.3), 12.2).

I4. The D_p-valued series (SW 3.5): with L_α = G + αH (G, H in Q_p[[T]]), L_p = Gω + a_p Hω − pHφ(ω); φ = F/p has
    characteristic polynomial X² − (a_p/p)X + 1/p. PARI returns D_p-valued data in the basis (ω, Fω).

I5. Pollack (Duke 2003): for a_p = 0, L_α = log⁺ L⁺ + α log⁻ L⁻ with L± in Λ and
    log⁺ = (1/p) ∏_{n≥1} Φ_{p^(2n)}(1 + T)/p, log⁻ = (1/p) ∏_{n≥1} Φ_{p^(2n−1)}(1 + T)/p.
    Integrality as stated in Kurihara–Pollack (eq. (3), citing Pollack's Thm 5.1) and in Burungale–Skinner–Tian–Wan
    (arXiv:2409.01350, section 6.1.1). Pollack's paper itself was not re-read.

I6. Modular symbols: PARI `msfromell(E, 1)` gives x⁺ normalised by Ω⁺ = E.omega[1], x⁺([0] − [∞]) = L(E, 1)/Ω⁺. We use
    [r]⁺ := x⁺([r] − [∞]). For 389a1, Ω_E = 2Ω⁺.

I7. ρ̄_p is surjective: E is semistable with ord_389(Δ) = 1, so inertia at 389 acts on E[p] through a transvection;
    E has no rational p-isogeny (`ellisomat`: the isogeny class is trivial; at supersingular p > 2, ρ̄_p is irreducible
    anyway), and an irreducible subgroup of GL₂(F_p) containing a transvection contains SL₂(F_p) (Serre 1972, section 2,
    Prop. 15; not re-read); the determinant is the cyclotomic character, which is surjective.

## 2. The measure and the Mazur–Tate elements

μ(a + pⁿZ_p) = α⁻ⁿ[a/pⁿ]⁺ − α^−(n+1)[a/p^(n−1)]⁺, α² = −p, and L_α(T) = ∫(1 + T)^s(x) dμ with ⟨x⟩ = (1 + p)^s(x).
For ζ of exact order p^m (m ≥ 1), the character x ↦ ζ^s(x) has conductor p^(m+1), so the α^−(n+1) term integrates to 0
and L_α(ζ − 1) = α^−(m+1) Θ_{m+1}(ζ − 1), where Θ_n(T) = Σ_{a ∈ (Z/pⁿ)ˣ} [a/pⁿ]⁺ (1 + T)^s(a).

## 3. The coefficients mod p

- m = 1: log⁻(ζ − 1) = 0, log⁺(ζ − 1) = 1/p and α⁻² = −1/p, so L⁺(ζ − 1) = −Θ₂(ζ − 1). Both sides are in Λ, so
  L⁺ + Θ₂ = Φ_p(1 + T)·g with g in Λ; Φ_p(1 + T) is Eisenstein, so c⁺_j ≡ −[T^j]Θ₂ (mod p) for j ≤ p − 2.
- m = 2: log⁺(ζ − 1) = 0, log⁻(ζ − 1) = Φ_p(ζ)/p² and α⁻³ = α/p², so Φ_p(ζ)L⁻(ζ − 1) = Θ₃(ζ − 1); hence
  c⁻_j ≡ [T^(p−1+j)]Θ₃ (mod p) for j ≤ p² − 2p.
- ord_T L_p ≥ 2 (I1 and I2), so G, H = O(T²) and L± = O(T²) (log±(0) ≠ 0). The T² coefficient of L_p is
  (c⁺₂/p)ω − (c⁻₂/p)Fω. In the basis (ω, Fω), φ has matrix [[0, −1], [1/p, 0]], and exactly
  (1 − φ)⁻²[(c⁺/p), −(c⁻/p)] = [−c⁺ + p(c⁺ + 2c⁻), 2c⁺ + c⁻ − pc⁻] / (1 + p)²,
  so its leading digits are [−c⁺₂, 2c⁺₂ + c⁻₂] mod p.
- S := (first coordinate of (1 − φ)⁻² L_p*) / (first coordinate of Reg_γ) mod p. By I1–I3 the ratio is
  c_∞ · d(0) · (unit) · #Sha(p) · ∏c_v (the factor c_∞ = 2 comes from Ω⁺ in place of Ω_E and is a unit), so S a
  p-adic unit gives #Sha(p) = 1. When c⁻₂ is computed, the second coordinate must satisfy the same ratio.

## 4. Memory

PARI's `ellpadicbsd` (overconvergent modular symbols) needs about 3.1 MB · p² for this curve at n = 2 (measured at
p = 11, 13, 17, 19, 23: 369 to 1625 MB), so about 35 GB at p = 107. The route above needs p² symbol values (p³ for the
second-coordinate check) and under 200 MB.

## 5. Literature check (22 Nov 2026)

Read: SW 2013 (Thm 12.3, Remarks 12.4 and 12.5); Perrin-Riou 2003 (arXiv:math/0109227, section 6: systematic
computations for conductor < 150 and p ≤ 19; 389a1 is out of that range); the Sage documentation (p = 5, 7 for 389a);
Burungale–Skinner–Tian–Wan section 1. Not read: Pollack 2003 itself; Sprung (Algebra Number Theory 2017);
Balakrishnan–Müller–Stein (Math. Comp. 2016); later supersingular algorithm papers. A web search found no computation
of Sha(389a)[p^∞] at a supersingular p. If you know of one, please tell us and we will credit it here.
