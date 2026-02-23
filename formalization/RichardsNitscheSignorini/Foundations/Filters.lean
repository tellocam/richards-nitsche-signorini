import Mathlib.Order.Filter.Basic
import Mathlib.Topology.Order.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# Endpoint Filters

Filters for asymptotic analysis at saturation endpoints of the effective saturation
domain `(0, 1)`. These filters encode the limiting behaviour as `Se → 0⁺` (dry end)
and `Se → 1⁻` (saturation end) and are used throughout the formalization to state
and prove asymptotic propositions from Sections 4.1–4.2 of the paper.

## Mathematical objects

- `atDryEnd` — filter on `ℝ` representing the limit `Se → 0⁺` (approach from above)
- `atSaturation` — filter on `ℝ` representing the limit `Se → 1⁻` (approach from below)

## References

- Paper, Section 4.1–4.2 (dry-end and wet-end asymptotic behaviour of `kr`, `pc`, etc.)

## Design

Both filters are defined via Mathlib's `nhdsWithin`, which gives the neighbourhood
filter of a point restricted to a set. This avoids introducing ad-hoc filter
infrastructure and connects directly to Mathlib's `Filter.Tendsto` API, so that
asymptotic statements can be expressed as plain `Tendsto f atDryEnd (nhds L)`.
-/

noncomputable section

namespace RichardsNitscheSignorini.Foundations

/-- Filter for `Se → 0⁺` (approaching the dry end from above).
Used for all dry-end asymptotic statements (Props 4.1, 4.5, Cor 4.6a). -/
noncomputable def atDryEnd : Filter ℝ := nhdsWithin 0 (Set.Ioi 0)

/-- Filter for `Se → 1⁻` (approaching saturation from below).
Used for wet-end asymptotic statements (Prop 4.2). -/
noncomputable def atSaturation : Filter ℝ := nhdsWithin 1 (Set.Iio 1)

end RichardsNitscheSignorini.Foundations
