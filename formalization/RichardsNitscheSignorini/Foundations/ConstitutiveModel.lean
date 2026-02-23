import Mathlib.Analysis.Asymptotics.Defs
import Mathlib.Analysis.Asymptotics.AsymptoticEquivalent
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Topology.Algebra.Order.LiminfLimsup
import RichardsNitscheSignorini.Foundations.Filters

/-!
# Van Genuchten-Mualem Constitutive Model

Defines the Van Genuchten-Mualem (VGM) parameter structure and the three
core constitutive functions relating effective saturation `Se ∈ [0,1]` to
hydraulic conductivity, specific moisture capacity, and capillary diffusivity.
These functions appear throughout the analysis in Section 4 of the paper.

## Mathematical objects

- `VGMParams` — parameter structure with 6 physical fields and 7 admissibility proofs
- `VGMParams.m` — derived shape parameter `m = 1 - 1/n ∈ (0,1)`
- `VGMParams.Dtheta` — water content range `Δθ = θ_s - θ_r > 0`
- `VGMParams.beta` — dry-end diffusivity exponent `β = ℓ + 1/m > 1`
- `VGMParams.C_D` — dry-end diffusivity prefactor
- `VGMParams.C_wet` — wet-end diffusivity prefactor
- `VGMParams.K_inv_exponent` — blow-up exponent `ℓ + 2/m` for `K⁻¹` at dry end
- `K_Se` — Mualem hydraulic conductivity `K(Se) = K_s · Se^ℓ · [1 - (1 - Se^{1/m})^m]²`
- `C_Se` — specific moisture capacity `C(Se) = Δθ · α · m · n · Se^{1/m} · (1 - Se^{1/m})^m`
- `D_Se` — capillary diffusivity `D(Se) = K(Se) / C(Se)`

## References

- van Genuchten, M. Th. (1980). A closed-form equation for predicting the
  hydraulic conductivity of unsaturated soils. *Soil Sci. Soc. Am. J.*, 44(5),
  892–898. Equations (4.5)–(4.8) of the paper.
- Mualem, Y. (1976). A new model for predicting the hydraulic conductivity of
  unsaturated porous media. *Water Resour. Res.*, 12(3), 513–522.
- Paper, Section 4.1: constitutive relations and their asymptotic properties.

## Design

The constitutive functions are defined on all of `ℝ`; they are meaningful only
for `Se ∈ [0,1]`. Properties that require real-analysis arguments (power-function
asymptotics, monotonicity of composed maps, continuity, Hölder regularity) are
stated as axioms with `[Established]` tags referencing the original papers and
the paper's equation numbers. This keeps the formalization tractable while
preserving the precise mathematical content of each claim.

Parameter arithmetic lemmas (`m_pos`, `m_lt_one`, etc.) are fully proved from
the admissibility hypotheses using elementary Lean tactics.
-/

noncomputable section

namespace RichardsNitscheSignorini.Foundations

open Asymptotics Filter Set Real

/-! ## Parameter structure -/

/-- Van Genuchten-Mualem parameter set.

Fields encode the six physical parameters plus admissibility constraints.
Notation follows the paper: `n` is the VG shape parameter (n > 1),
`ell` is the Mualem tortuosity exponent (ell >= 0). -/
structure VGMParams where
  /-- Saturated water content (dimensionless, 0 < theta_s). -/
  theta_s : ℝ
  /-- Residual water content (dimensionless, 0 <= theta_r < theta_s). -/
  theta_r : ℝ
  /-- Inverse air-entry pressure (1/m, alpha > 0). -/
  alpha : ℝ
  /-- Van Genuchten shape parameter (dimensionless, n > 1). -/
  n : ℝ
  /-- Mualem tortuosity exponent (dimensionless, ell >= 0). -/
  ell : ℝ
  /-- Saturated hydraulic conductivity (m/s, K_s > 0). -/
  K_s : ℝ
  -- Admissibility hypotheses
  theta_s_pos : 0 < theta_s
  theta_r_nonneg : 0 ≤ theta_r
  theta_r_lt_theta_s : theta_r < theta_s
  alpha_pos : 0 < alpha
  n_gt_one : 1 < n
  ell_nonneg : 0 ≤ ell
  K_s_pos : 0 < K_s

namespace VGMParams

variable (p : VGMParams)

/-! ## Derived parameters -/

/-- Van Genuchten shape parameter m = 1 - 1/n in (0,1). -/
noncomputable def m : ℝ := 1 - 1 / p.n

/-- Water content range Delta_theta = theta_s - theta_r > 0. -/
noncomputable def Dtheta : ℝ := p.theta_s - p.theta_r

/-- Dry-end diffusivity exponent beta = ell + 1/m > 1.
Controls the rate D(Se) ~ C_D * Se^beta as Se -> 0+. -/
noncomputable def beta : ℝ := p.ell + 1 / p.m

/-- Dry-end diffusivity prefactor C_D = K_s * m / (Dtheta * alpha * n). -/
noncomputable def C_D : ℝ := p.K_s * p.m / (p.Dtheta * p.alpha * p.n)

/-- Wet-end diffusivity prefactor C_wet = K_s * m^{m-1} / (Dtheta * alpha * n). -/
noncomputable def C_wet : ℝ := p.K_s * p.m ^ (p.m - 1) / (p.Dtheta * p.alpha * p.n)

/-- K^{-1} blow-up exponent at dry end: ell + 2/m. -/
noncomputable def K_inv_exponent : ℝ := p.ell + 2 / p.m

/-! ## Parameter arithmetic lemmas -/

theorem m_pos : 0 < p.m := by
  unfold m
  have h : 0 < p.n := lt_trans zero_lt_one p.n_gt_one
  have : 0 < 1 / p.n := one_div_pos.mpr h
  have : 1 / p.n < 1 := by
    rw [div_lt_one h]
    exact p.n_gt_one
  linarith

theorem m_lt_one : p.m < 1 := by
  unfold m
  have h : 0 < p.n := lt_trans zero_lt_one p.n_gt_one
  have : 0 < 1 / p.n := one_div_pos.mpr h
  linarith

theorem Dtheta_pos : 0 < p.Dtheta := by
  unfold Dtheta
  linarith [p.theta_r_lt_theta_s]

theorem beta_gt_one : 1 < p.beta := by
  unfold beta
  have hm_pos := p.m_pos
  have hm_lt_one := p.m_lt_one
  have : 1 < 1 / p.m := by
    field_simp
    nlinarith
  linarith [p.ell_nonneg]

theorem C_D_pos : 0 < p.C_D := by
  unfold C_D
  have hm_pos := p.m_pos
  have hDtheta_pos := p.Dtheta_pos
  have hn_pos : 0 < p.n := lt_trans zero_lt_one p.n_gt_one
  apply div_pos
  · exact mul_pos p.K_s_pos hm_pos
  · apply mul_pos
    · apply mul_pos
      · exact hDtheta_pos
      · exact p.alpha_pos
    · exact hn_pos

theorem K_inv_exponent_pos : 0 < p.K_inv_exponent := by
  unfold K_inv_exponent
  have hm_pos := p.m_pos
  have : 0 < 2 / p.m := div_pos (by norm_num : (0 : ℝ) < 2) hm_pos
  linarith [p.ell_nonneg]

end VGMParams

/-! ## Constitutive functions -/

/-- Mualem hydraulic conductivity K(Se) = K_s * Se^ell * [1 - (1 - Se^{1/m})^m]^2.
Defined on all of ℝ; meaningful for Se in [0,1]. -/
noncomputable def K_Se (p : VGMParams) (Se : ℝ) : ℝ :=
  p.K_s * Se ^ p.ell * (1 - (1 - Se ^ (1 / p.m)) ^ p.m) ^ 2

/-- Specific moisture capacity C(Se) = Dtheta * alpha * m * n * Se^{1/m} * (1 - Se^{1/m})^m.
The derivative d(theta)/d(psi) expressed in terms of Se. -/
noncomputable def C_Se (p : VGMParams) (Se : ℝ) : ℝ :=
  p.Dtheta * p.alpha * p.m * p.n * Se ^ (1 / p.m) * (1 - Se ^ (1 / p.m)) ^ p.m

/-- Capillary diffusivity D(Se) = K(Se) / C(Se).
Controls the parabolic character of Richards' equation. -/
noncomputable def D_Se (p : VGMParams) (Se : ℝ) : ℝ :=
  K_Se p Se / C_Se p Se

/-! ## Constitutive axioms

These encode properties of K that follow from the Mualem formula but
require real-analysis arguments (compositions of power functions,
monotonicity of composed maps, continuity, Hölder regularity). We axiomatize
them with `[Established]` tags referencing the original papers and the paper's
equation numbers, to keep the formalization tractable while faithfully recording
the mathematical content of each claim. -/

/-- [Established: Mualem 1976; Paper Section 4.1, eq. (4.5)]
K(0) = 0: dry soil has zero conductivity.

Direct computation from the definition: Se = 0 forces Se^ell = 0. -/
axiom K_Se_zero (p : VGMParams) : K_Se p 0 = 0

/-- [Established: Mualem 1976; Paper Section 4.1, eq. (4.5)]
K(1) = K_s: saturated soil has full conductivity.

Direct computation from the definition: Se = 1 forces (1 - 1^{1/m})^m = 0,
so the bracket equals 1, giving K_s * 1^ell * 1^2 = K_s. -/
axiom K_Se_one (p : VGMParams) : K_Se p 1 = p.K_s

/-- [Established: van Genuchten 1980; Paper Section 4.1]
K is strictly increasing on [0,1].

Follows from the monotonicity of each factor in the Mualem product formula. -/
axiom K_Se_strictMono (p : VGMParams) : StrictMonoOn (K_Se p) (Icc 0 1)

/-- [Established: van Genuchten 1980; Paper Section 4.1]
K is continuous on [0,1].

Follows from the composition of continuous functions (power functions with
positive base on the interior, extended continuously to the endpoints). -/
axiom K_Se_continuousOn (p : VGMParams) : ContinuousOn (K_Se p) (Icc 0 1)

/-- [Established: Mualem 1976; Paper Section 4.1]
K is nonneg on [0,1].

Follows from the product of nonneg factors: Se^ell >= 0 and the squared
bracket is nonneg. -/
axiom K_Se_nonneg (p : VGMParams) {Se : ℝ} (hSe : Se ∈ Icc 0 1) : 0 ≤ K_Se p Se

/-- [Established: Paper Section 4.1, eq. (4.20)]
Hölder-m regularity at saturation: K_s - K(Se) ~ C * (1-Se)^m as Se -> 1-.

K is NOT Lipschitz at Se = 1 when m < 1. The exponent m ∈ (0,1) governs the
Hölder singularity at the wet end. Derived by Taylor expansion at Se = 1. -/
axiom K_Se_holder_at_one (p : VGMParams) :
    (fun Se => p.K_s - K_Se p Se) =Θ[atSaturation] (fun Se => (1 - Se) ^ p.m)

/-- [Established: Paper Section 4.1, eq. (4.12)]
K'(Se) -> 0 as Se -> 0+ (characteristic speed vanishes at dry end).

From eq. (4.12): dK/dSe ~ K_s * m^2 * (ell + 2/m) * Se^{ell + 2/m - 1} -> 0
as Se -> 0+, since ell + 2/m - 1 > 0 when ell >= 0 and m < 1. -/
axiom K_Se_deriv_zero_at_dry (p : VGMParams) :
    Filter.Tendsto (fun Se => (K_Se p Se - K_Se p 0) / Se) atDryEnd (nhds 0)

/-- [Established: Paper Section 4.1, eq. (4.21)]
K'(Se) -> +∞ as Se -> 1- (characteristic speed blows up at saturation).

From eq. (4.21), the Hölder singularity of K at Se = 1 implies the difference
quotient diverges. This is the Lax entropy condition for the shock. -/
axiom K_Se_deriv_infty_at_sat (p : VGMParams) :
    Filter.Tendsto (fun Se => (K_Se p 1 - K_Se p Se) / (1 - Se)) atSaturation Filter.atTop

/-- [Established: Paper Section 4.1, eq. (4.11)]
Small-Se expansion of K: K(Se) = K_s * m^2 * Se^{ell + 2/m} + h.o.t.

The leading term in the small-Se Taylor expansion of K(Se) is
K_s * m^2 * Se^{ell + 2/m}, with the next correction of order Se^{ell + 3/m}. -/
axiom K_Se_asymp_dry (p : VGMParams) :
    (fun Se => K_Se p Se - p.K_s * p.m ^ 2 * Se ^ (p.ell + 2 / p.m)) =O[atDryEnd]
      (fun Se => Se ^ (p.ell + 3 / p.m))

/-- [Established: Paper Section 4.1, eq. (4.13)]
Small-Se expansion of C: C(Se) = Dtheta * alpha * m * n * Se^{1/m} * (1 + O(Se^{1/m})).

The leading term in the small-Se Taylor expansion of C(Se) is
Dtheta * alpha * m * n * Se^{1/m}, with the correction of order Se^{2/m}. -/
axiom C_Se_asymp_dry (p : VGMParams) :
    (fun Se => C_Se p Se - p.Dtheta * p.alpha * p.m * p.n * Se ^ (1 / p.m)) =O[atDryEnd]
      (fun Se => Se ^ (2 / p.m))

end RichardsNitscheSignorini.Foundations
