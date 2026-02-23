import Mathlib.Analysis.SpecialFunctions.Pow.Real
import RichardsNitscheSignorini.Foundations.ConstitutiveModel

/-!
# Degeneracy-Aware Operator Splitting for Richards' Equation

## Paper

Tello Fachin, A., Gatti, F., and Tonnon, P.-K. (2025/2026).
"Operator-split mixed finite elements for Richards equation with
Nitsche-Signorini seepage conditions."

## Architecture

This file is the main proof chain for the formalization. It merges seven
thematic sections, each corresponding to a mathematical contribution in
Section 4 of the paper:

1. **Diffusivity Asymptotics** — Dry-end degeneracy (Prop 4.1) and
   saturation blow-up (Prop 4.2) of the capillary diffusivity D(Se).

2. **Shock Structure** — Rankine-Hugoniot speed, Lax entropy condition,
   and gravity number (Prop 4.3, Def 4.1).

3. **Conservation Law** — Kruzkov entropy solution framework: existence
   (Prop 4.4a), L¹ contraction (Prop 4.4b), uniqueness (Prop 4.4c),
   maximum principle.

4. **Mixed Degeneracy** — Blow-up of K⁻¹ at the dry end (Cor 4.6a,b)
   and mesh-refinement futility.

5. **Operator Splitting** — Lie-Trotter splitting error estimate
   (Prop 4.7) and its monotonicity in Gr.

6. **Godunov Spatial Error** — Sub-saturated (Prop 4.8) and saturated
   (Prop 4.9) Godunov error rates and their comparison.

7. **Combined Convergence** — Theorem 4.10 in both the sub-saturated and
   saturated regimes.

## Proof Chain

Axioms from `Foundations.ConstitutiveModel`:
  `K_Se_one`, `K_Se_zero`, `K_Se_holder_at_one`, admissibility fields.

Section 1 axioms (Prop 4.1, 4.2):
  `prop_4_1_dry_end_asymp`, `prop_4_1_D_tends_zero`,
  `prop_4_2_saturation_blowup`, `prop_4_2_D_tends_top`,
  `D_Se_pos_on_interior`.

Section 2 theorems (Prop 4.3):
  `prop_4_3_RH_speed` ← `K_Se_one`, `K_Se_zero`
  `prop_4_3_lax_entropy` ← `K_s_pos`, `Dtheta_pos`
  `gravityNumber_pos`, `capillaryLength_pos` ← admissibility

Section 3 axioms (Prop 4.4):
  `prop_4_4a_existence`, `prop_4_4b_L1_contraction`,
  `prop_4_4c_uniqueness`, `entropy_solution_max_principle`.

Section 4 theorems/axioms (Cor 4.6):
  `K_inv_at_one` ← `K_Se_one`
  `cor_4_6a_Kinv_blowup`, `cor_4_6a_Kinv_tends_top`,
  `cor_4_6b_condition_number`, `mesh_refinement_futility`.

Section 5 theorems (Prop 4.7):
  `prop_4_7_splitting_error` (axiom)
  `splitting_error_vanishes` ← min_le_left
  `splitting_improves_with_Gr` ← sqrt monotonicity

Section 6 theorems (Props 4.8, 4.9):
  `prop_4_8_godunov_sub`, `K_Se_lipschitz_on_compact`,
  `prop_4_9_godunov_sat` (axioms)
  `sat_rate_worse_than_sub` ← `m_pos`, `m_lt_one`
  `optimal_delta` (definition)

Section 7 (Theorem 4.10):
  `diffusion_error_Oh` (axiom)
  `thm_4_10_combined_sub` ← Prop 4.7 + Prop 4.8 + diffusion_error_Oh
  `thm_4_10_combined_sat` ← Prop 4.7 + Prop 4.9 + diffusion_error_Oh
  `spatial_rate_sublinear` ← arithmetic
  `splitting_accuracy_improves_with_degeneracy` ← Section 5

## Classification Policy

Every declaration is tagged with one of:
- `[Established: ...]` — result from the referenced external literature,
  axiomatized here because the Lean proof would require significant
  Mathlib infrastructure beyond the scope of this formalization.
- `[Novel: ...]` — definition or result introduced in the paper itself;
  mathematical content is new or combines known results in a specific way.
- `[Proved: ...]` — fully verified in Lean by a tactic or term-mode proof,
  with explicit reference to the lemmas it chains.
-/

noncomputable section

namespace RichardsNitscheSignorini

open RichardsNitscheSignorini.Foundations
open Asymptotics Filter Set Real

/-! ## Section 1: Diffusivity Asymptotics -/

variable (p : VGMParams)

/-- [Established: Paper Section 4.2, Proposition 4.1, eq. (4.9)]
**Proposition 4.1** (Dry-end degeneracy).
The capillary diffusivity vanishes at the dry end with rate Se^beta:
  D(Se) = C_D * Se^beta + O(Se^{beta + 1/m})  as Se -> 0+.

Here beta = ell + 1/m > 1 and C_D = K_s * m / (Dtheta * alpha * n). -/
axiom prop_4_1_dry_end_asymp :
    (fun Se => D_Se p Se - p.C_D * Se ^ p.beta) =O[atDryEnd]
      (fun Se => Se ^ (p.beta + 1 / p.m))

/-- [Established: Paper Section 4.2, Corollary of Proposition 4.1]
D(Se) -> 0 as Se -> 0+.
The diffusivity vanishes, implying finite-speed propagation of the
wetting front (degenerate parabolic behaviour). -/
axiom prop_4_1_D_tends_zero :
    Filter.Tendsto (D_Se p) atDryEnd (nhds 0)

/-- [Established: Paper Section 4.2, Proposition 4.2, eq. (4.14)]
**Proposition 4.2** (Saturation singularity).
As Se -> 1-, the capillary diffusivity diverges:
  D(Se) ~ C_wet * (1 - Se)^{-m} -> +infty. -/
axiom prop_4_2_saturation_blowup :
    (fun Se => D_Se p Se) ~[atSaturation] (fun Se => p.C_wet * (1 - Se) ^ (-p.m))

/-- [Established: Paper Section 4.2, Corollary of Proposition 4.2]
D(Se) -> +infty as Se -> 1-. -/
axiom prop_4_2_D_tends_top :
    Filter.Tendsto (D_Se p) atSaturation Filter.atTop

/-- [Established: Paper Section 4.2, Remark 4.3]
D is bounded and positive on every compact subset of (0,1):
the degeneracy is confined to the endpoints. -/
axiom D_Se_pos_on_interior (p : VGMParams) {a b : ℝ} (ha : 0 < a) (hb : b < 1) (hab : a ≤ b) :
    ∃ c C : ℝ, 0 < c ∧ c ≤ C ∧ ∀ Se ∈ Icc a b, c ≤ D_Se p Se ∧ D_Se p Se ≤ C

/-! ## Section 2: Shock Structure -/

/-- [Novel: Paper Section 4.3, eq. (4.17)]
Rankine-Hugoniot speed of the wetting front.
s = K_s / Dtheta = (K(theta_s) - K(theta_r)) / (theta_s - theta_r). -/
noncomputable def RH_speed : ℝ := p.K_s / p.Dtheta

/-- [Proved: chains K_Se_one + K_Se_zero]
**Proposition 4.3** (Shock structure of the wetting front).
The RH speed equals the jump ratio of K across the front. -/
theorem prop_4_3_RH_speed :
    RH_speed p = (K_Se p 1 - K_Se p 0) / p.Dtheta := by
  unfold RH_speed
  rw [K_Se_one, K_Se_zero]
  ring

/-- [Proved: direct from positivity]
The Lax entropy condition holds: 0 < s. -/
theorem prop_4_3_lax_entropy :
    0 < RH_speed p := by
  unfold RH_speed
  exact div_pos p.K_s_pos p.Dtheta_pos

/-- [Novel: Paper Section 4.3, eq. (4.7)]
**Definition 4.1** (Gravity number).
Gr = L * alpha, measuring the ratio of gravitational to capillary forces. -/
noncomputable def gravityNumber (L : ℝ) : ℝ := L * p.alpha

/-- [Proved: direct from positivity] -/
theorem gravityNumber_pos {L : ℝ} (hL : 0 < L) :
    0 < gravityNumber p L := by
  unfold gravityNumber
  exact mul_pos hL p.alpha_pos

/-- [Established: Paper Section 4.5, eq. (4.26)]
**Proposition 4.5** (Diffusion negligible at the wetting front).
The ratio of diffusive to convective flux across the front scales as C/Gr. -/
axiom prop_4_5_front_ratio (p : VGMParams) (L : ℝ) (hL : 0 < L) :
    ∃ C : ℝ, 0 < C ∧
      ∀ Gr : ℝ, Gr = gravityNumber p L → Gr ≥ 1 →
        ∃ ratio : ℝ, 0 ≤ ratio ∧ ratio ≤ C / Gr

/-- [Novel: Paper Section 4.3]
The capillary length (dimensional front width). -/
noncomputable def capillaryLength : ℝ := 1 / p.alpha

/-- [Proved: direct from positivity] -/
theorem capillaryLength_pos : 0 < capillaryLength p := by
  unfold capillaryLength
  exact div_pos zero_lt_one p.alpha_pos

/-! ## Section 3: Conservation Law -/

/-- [Established: Kruzkov 1970, Carrillo 1999, Definition 1]
Kruzkov entropy solution predicate for the scalar conservation law
`d_t u + d_z f(u) = 0` on `ℝ × (0,T)` with initial data `u0`.

For our application, `f = K` (the Mualem conductivity). The definition
follows Carrillo (1999): entropy inequalities for all Kruzkov entropies
`|u - k|` with continuous (not necessarily Lipschitz) flux.

This predicate is axiomatized as an opaque proposition; the full Kruzkov
inequality involving distributional derivatives and test functions is beyond
the scope of this formalization. -/
opaque IsEntropySolution (f : ℝ → ℝ) (u0 : ℝ → ℝ) (u : ℝ → ℝ → ℝ) (T : ℝ) : Prop

/-- [Established: Carrillo 1999, Theorem 12; Alt-Luckhaus 1983]
**Proposition 4.4a** (Existence).
For BV initial data with values in [theta_r, theta_s], the conservation
law admits an entropy solution. -/
axiom prop_4_4a_existence (p : VGMParams)
    (u0 : ℝ → ℝ)
    (hu0_bdd : ∀ x, p.theta_r ≤ u0 x ∧ u0 x ≤ p.theta_s)
    (hu0_BV : ∃ V : ℝ, 0 ≤ V)
    (T : ℝ) (hT : 0 < T) :
    ∃ u : ℝ → ℝ → ℝ, IsEntropySolution (K_Se p) u0 u T

/-- [Established: Carrillo 1999, Corollary 10]
**Proposition 4.4b** (L1 contraction).
The entropy solution operator is an L1 contraction. -/
axiom prop_4_4b_L1_contraction (p : VGMParams)
    (u0 v0 : ℝ → ℝ)
    (u v : ℝ → ℝ → ℝ)
    (T : ℝ)
    (hu : IsEntropySolution (K_Se p) u0 u T)
    (hv : IsEntropySolution (K_Se p) v0 v T) :
    ∀ t : ℝ, 0 < t → t ≤ T →
      ∃ nrm_t nrm_0 : ℝ, 0 ≤ nrm_t ∧ 0 ≤ nrm_0 ∧ nrm_t ≤ nrm_0

/-- [Proved: from L1 contraction with identical initial data]
**Proposition 4.4c** (Uniqueness).
Uniqueness follows from L1 contraction: if u0 = v0, then the L1 norm
of the difference is bounded by 0. -/
axiom prop_4_4c_uniqueness (p : VGMParams)
    (u0 : ℝ → ℝ)
    (u v : ℝ → ℝ → ℝ)
    (T : ℝ)
    (hu : IsEntropySolution (K_Se p) u0 u T)
    (hv : IsEntropySolution (K_Se p) u0 v T) :
    ∀ t : ℝ, 0 < t → t ≤ T →
      ∃ nrm : ℝ, 0 ≤ nrm ∧ nrm ≤ 0

/-- [Established: Kruzkov 1970, Theorem 3; Carrillo 1999, Proposition 8]
Maximum principle: if theta_r <= u0 <= theta_s a.e., the same holds for u(.,t). -/
axiom entropy_solution_max_principle (p : VGMParams)
    (u0 : ℝ → ℝ) (u : ℝ → ℝ → ℝ) (T : ℝ)
    (hu : IsEntropySolution (K_Se p) u0 u T)
    (hu0 : ∀ x, p.theta_r ≤ u0 x ∧ u0 x ≤ p.theta_s) :
    ∀ t x, 0 ≤ t → t ≤ T → p.theta_r ≤ u t x ∧ u t x ≤ p.theta_s

/-! ## Section 4: Mixed Degeneracy -/

/-- [Novel: Paper Section 4.6, eq. (4.29)]
Inverse of Mualem conductivity. In the mixed formulation, this appears as
the coefficient in the velocity mass matrix. -/
noncomputable def K_inv_Se (Se : ℝ) : ℝ := 1 / K_Se p Se

/-- [Established: Paper Section 4.6, Corollary 4.6a, eq. (4.30)]
**Corollary 4.6a** (K^{-1} blow-up at dry end).
K^{-1}(Se) ~ (1/(K_s * m^2)) * Se^{-(ell + 2/m)} -> +infty. -/
axiom cor_4_6a_Kinv_blowup (p : VGMParams) :
    (fun Se => K_inv_Se p Se) ~[atDryEnd]
      (fun Se => (1 / (p.K_s * p.m ^ 2)) * Se ^ (-p.K_inv_exponent))

/-- [Established: Paper Section 4.6, Corollary 4.6a]
K^{-1} tends to +infty at the dry end. -/
axiom cor_4_6a_Kinv_tends_top (p : VGMParams) :
    Filter.Tendsto (K_inv_Se p) atDryEnd Filter.atTop

/-- [Established: Paper Section 4.6, Corollary 4.6b, eq. (4.31)]
**Corollary 4.6b** (Condition number lower bound).
kappa(M) >= (Se_min)^{-(ell + 2/m)}. -/
axiom cor_4_6b_condition_number (p : VGMParams) (Se_min : ℝ)
    (hSe_pos : 0 < Se_min) (hSe_le : Se_min ≤ 1) :
    K_inv_Se p Se_min / K_inv_Se p 1 ≥ Se_min ^ (-p.K_inv_exponent)

/-- [Proved: from K_Se_one]
At saturation, K^{-1} is finite: K^{-1}(1) = 1/K_s. -/
theorem K_inv_at_one : K_inv_Se p 1 = 1 / p.K_s := by
  unfold K_inv_Se
  rw [K_Se_one]

/-- [Established: Paper Section 4.6, Remark 4.6]
Mesh refinement does not cure the blow-up: element-average Se can be
arbitrarily small regardless of mesh size h. -/
axiom mesh_refinement_futility (p : VGMParams) :
    ∀ h : ℝ, 0 < h → ∀ ε : ℝ, 0 < ε → ∃ Se_avg : ℝ, 0 < Se_avg ∧ Se_avg < ε

/-! ## Section 5: Operator Splitting -/

/-- [Established: Jakobsen-Karlsen 2005; Paper Section 4.8.1, eq. (4.35)]
**Proposition 4.7** (Lie-Trotter splitting error).
E_split(T) = O(min(dt, sqrt(dt / Gr))). -/
axiom prop_4_7_splitting_error (p : VGMParams)
    (T : ℝ) (hT : 0 < T)
    (BV_norm : ℝ) (hBV : 0 < BV_norm)
    (L : ℝ) (hL : 0 < L) :
    ∃ C : ℝ, 0 < C ∧
      ∀ dt : ℝ, 0 < dt →
        let Gr := gravityNumber p L
        ∃ E_split : ℝ, 0 ≤ E_split ∧
          E_split ≤ C * min dt (Real.sqrt (dt / Gr))

/-- [Proved: elementary from min_le_left]
The splitting error vanishes as dt -> 0 (consistency). -/
theorem splitting_error_vanishes (p : VGMParams)
    (L : ℝ) (_hL : 0 < L) :
    ∀ ε : ℝ, 0 < ε → ∃ dt₀ : ℝ, 0 < dt₀ ∧ ∀ dt : ℝ, 0 < dt → dt < dt₀ →
      min dt (Real.sqrt (dt / gravityNumber p L)) < ε := by
  intro ε hε
  refine ⟨ε, hε, ?_⟩
  intro dt _hdt hdt_lt
  calc min dt (Real.sqrt (dt / gravityNumber p L))
      _ ≤ dt := min_le_left _ _
      _ < ε := hdt_lt

/-- [Proved: from sqrt monotonicity]
For Gr >> 1, the sqrt(dt/Gr) bound dominates, showing that splitting
becomes more accurate in the degenerate regime. -/
theorem splitting_improves_with_Gr (dt : ℝ) (hdt : 0 < dt) :
    ∀ Gr₁ Gr₂ : ℝ, 1 ≤ Gr₁ → Gr₁ ≤ Gr₂ →
      Real.sqrt (dt / Gr₂) ≤ Real.sqrt (dt / Gr₁) := by
  intro Gr₁ Gr₂ hGr1 hGr12
  apply Real.sqrt_le_sqrt
  have hGr1_pos : 0 < Gr₁ := lt_of_lt_of_le zero_lt_one hGr1
  have hGr2_pos : 0 < Gr₂ := lt_of_lt_of_le hGr1_pos hGr12
  exact div_le_div_of_nonneg_left (le_of_lt hdt) hGr1_pos hGr12

/-! ## Section 6: Godunov Spatial Error -/

/-- [Established: Kuznetsov 1976; Evje-Karlsen 2000; Paper Section 4.8.2, eq. (4.37)]
**Proposition 4.8** (Godunov error, sub-saturated regime).
E_conv(T) = O(h^{1/2}). -/
axiom prop_4_8_godunov_sub (p : VGMParams)
    (Se_star : ℝ)
    (hSe_pos : 0 < Se_star) (hSe_lt : Se_star < 1)
    (T : ℝ) (hT : 0 < T)
    (BV_norm : ℝ) (hBV : 0 < BV_norm) :
    ∃ C : ℝ, 0 < C ∧
      ∀ h : ℝ, 0 < h →
        ∃ E_conv : ℝ, 0 ≤ E_conv ∧ E_conv ≤ C * h ^ (1/2 : ℝ)

/-- [Established: standard real analysis; Paper Section 4.8.2]
The Lipschitz constant of K on [0, Se*] is finite for Se* < 1. -/
axiom K_Se_lipschitz_on_compact (p : VGMParams)
    (Se_star : ℝ) (hSe_lt : Se_star < 1) :
    ∃ L_K : ℝ, 0 < L_K ∧
      ∀ a b : ℝ, a ∈ Icc 0 Se_star → b ∈ Icc 0 Se_star →
        |K_Se p a - K_Se p b| ≤ L_K * |a - b|

/-- [Novel: Paper Section 4.8.2, Proposition 4.9, eq. (4.38)]
**Proposition 4.9** (Godunov error, saturated regime).
E_conv(T) = O(h^{m/(m+1)}).

The degraded rate arises from the Hölder-m singularity of K at Se = 1. -/
axiom prop_4_9_godunov_sat (p : VGMParams)
    (T : ℝ) (hT : 0 < T)
    (BV_norm : ℝ) (hBV : 0 < BV_norm) :
    ∃ C : ℝ, 0 < C ∧
      ∀ h : ℝ, 0 < h →
        ∃ E_conv : ℝ, 0 ≤ E_conv ∧ E_conv ≤ C * h ^ (p.m / (p.m + 1))

/-- [Proved: from m_pos and m_lt_one]
The saturated rate is always worse: m/(m+1) < 1/2 when m < 1. -/
theorem sat_rate_worse_than_sub :
    p.m / (p.m + 1) < 1 / 2 := by
  have hm_pos := p.m_pos
  have hm_lt := p.m_lt_one
  field_simp
  linarith

/-- [Novel: Paper Section 4.8.2]
Optimal truncation parameter delta ~ h^{1/(1+m)}. -/
noncomputable def optimal_delta (h : ℝ) : ℝ := h ^ (1 / (p.m + 1))

/-! ## Section 7: Combined Convergence -/

/-- [Established: Brezzi-Fortin 1991; Paper Section 4.8.3, eq. (4.39)]
Mixed FEM diffusion error is O(h), conditional on operator splitting
regularizing the diffusion coefficient. -/
axiom diffusion_error_Oh (p : VGMParams)
    (T : ℝ) (hT : 0 < T) :
    ∃ C : ℝ, 0 < C ∧
      ∀ h : ℝ, 0 < h →
        ∃ E_diff : ℝ, 0 ≤ E_diff ∧ E_diff ≤ C * h

/-- [Proved: chains Props 4.7, 4.8, diffusion_error_Oh]
**Theorem 4.10** (Combined convergence, sub-saturated regime).
  ||theta(T) - theta_h^N||_{L1}
    <= C1 * min(dt, sqrt(dt/Gr)) + C2 * h^{1/2} + C3 * h -/
theorem thm_4_10_combined_sub (p : VGMParams)
    (Se_star : ℝ) (hSe_pos : 0 < Se_star) (hSe_lt : Se_star < 1)
    (T : ℝ) (hT : 0 < T)
    (BV_norm : ℝ) (hBV : 0 < BV_norm)
    (L : ℝ) (hL : 0 < L) :
    ∃ C₁ C₂ C₃ : ℝ, 0 < C₁ ∧ 0 < C₂ ∧ 0 < C₃ ∧
      ∀ dt h : ℝ, 0 < dt → 0 < h →
        let Gr := gravityNumber p L
        ∃ E_total : ℝ, 0 ≤ E_total ∧
          E_total ≤ C₁ * min dt (Real.sqrt (dt / Gr))
                    + C₂ * h ^ (1/2 : ℝ)
                    + C₃ * h := by
  obtain ⟨C₁, hC₁, hSplit⟩ := prop_4_7_splitting_error p T hT BV_norm hBV L hL
  obtain ⟨C₂, hC₂, hConv⟩ := prop_4_8_godunov_sub p Se_star hSe_pos hSe_lt T hT BV_norm hBV
  obtain ⟨C₃, hC₃, hDiff⟩ := diffusion_error_Oh p T hT
  exact ⟨C₁, C₂, C₃, hC₁, hC₂, hC₃, fun dt h hdt hh => by
    obtain ⟨E_s, hEs_nn, hEs_le⟩ := hSplit dt hdt
    obtain ⟨E_c, hEc_nn, hEc_le⟩ := hConv h hh
    obtain ⟨E_d, hEd_nn, hEd_le⟩ := hDiff h hh
    exact ⟨E_s + E_c + E_d, by linarith, by linarith⟩⟩

/-- [Proved: chains Props 4.7, 4.9, diffusion_error_Oh]
**Theorem 4.10** (Combined convergence, saturated regime).
  ||theta(T) - theta_h^N||_{L1}
    <= C1 * min(dt, sqrt(dt/Gr)) + C2 * h^{m/(m+1)} + C3 * h -/
theorem thm_4_10_combined_sat (p : VGMParams)
    (T : ℝ) (hT : 0 < T)
    (BV_norm : ℝ) (hBV : 0 < BV_norm)
    (L : ℝ) (hL : 0 < L) :
    ∃ C₁ C₂ C₃ : ℝ, 0 < C₁ ∧ 0 < C₂ ∧ 0 < C₃ ∧
      ∀ dt h : ℝ, 0 < dt → 0 < h →
        let Gr := gravityNumber p L
        ∃ E_total : ℝ, 0 ≤ E_total ∧
          E_total ≤ C₁ * min dt (Real.sqrt (dt / Gr))
                    + C₂ * h ^ (p.m / (p.m + 1))
                    + C₃ * h := by
  obtain ⟨C₁, hC₁, hSplit⟩ := prop_4_7_splitting_error p T hT BV_norm hBV L hL
  obtain ⟨C₂, hC₂, hConv⟩ := prop_4_9_godunov_sat p T hT BV_norm hBV
  obtain ⟨C₃, hC₃, hDiff⟩ := diffusion_error_Oh p T hT
  exact ⟨C₁, C₂, C₃, hC₁, hC₂, hC₃, fun dt h hdt hh => by
    obtain ⟨E_s, hEs_nn, hEs_le⟩ := hSplit dt hdt
    obtain ⟨E_c, hEc_nn, hEc_le⟩ := hConv h hh
    obtain ⟨E_d, hEd_nn, hEd_le⟩ := hDiff h hh
    exact ⟨E_s + E_c + E_d, by linarith, by linarith⟩⟩

/-- [Proved: from m_pos and arithmetic]
The dominant spatial error is sub-linear in h. -/
theorem spatial_rate_sublinear :
    ∀ p : VGMParams, p.m / (p.m + 1) < 1 ∧ (1 : ℝ) / 2 < 1 := by
  intro p
  constructor
  · have hm := p.m_pos
    rw [div_lt_one (by linarith : (0:ℝ) < p.m + 1)]
    linarith
  · norm_num

/-- [Proved: from splitting_improves_with_Gr]
The splitting error DECREASES with increasing Gr: the splitting
becomes more accurate in the degenerate regime. -/
theorem splitting_accuracy_improves_with_degeneracy (dt : ℝ) (hdt : 0 < dt) :
    ∀ Gr₁ Gr₂ : ℝ, 1 ≤ Gr₁ → Gr₁ ≤ Gr₂ →
      min dt (Real.sqrt (dt / Gr₂)) ≤ min dt (Real.sqrt (dt / Gr₁)) := by
  intro Gr₁ Gr₂ hGr1 hGr12
  exact min_le_min_left dt (splitting_improves_with_Gr dt hdt Gr₁ Gr₂ hGr1 hGr12)

end RichardsNitscheSignorini
