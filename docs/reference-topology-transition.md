# Reference topology and four-vector transition

This is a coordinated breaking change in ThreeBodyDecays and ThreeBodyDecaysIO.
ThreeBodyDecaysIO 0.9 remains available for existing angle-based documents. The
new IO API targets ThreeBodyDecays 0.16 and requires Julia 1.10 or later.

## Ownership

- ThreeBodyDecays owns the amplitude convention and stores `reference_k` on each
  `ThreeBodyDecay`. Its value is part of the model definition, not an evaluation
  option. All chains contributing coherently use that model reference.
- ThreeBodyDecaysIO owns phase-space domain parsing, indexed four-vector input,
  extraction of three-body angles and invariants, and document validation.
- FourVectors supplies boosts and rotations. InstructionalDecayTrees is an
  independent reference for tests, not a runtime dependency.
- CascadeDecaysIO remains the separate implementation for general cascades.

## Core prerequisites

The model-state migration is reviewed in
[ThreeBodyDecays PR #128](https://github.com/mmikhasenko/ThreeBodyDecays.jl/pull/128).
While that API is unreleased, the IO test workflow pins its exact development
revision instead of treating any package labeled 0.16 as sufficient.

[ThreeBodyDecays PR #127](https://github.com/mmikhasenko/ThreeBodyDecays.jl/pull/127)
corrects the full-turn phase of parent rotations between chains 2 and 3. For a
half-integer parent spin, omitting this phase can change coherent intensities
when changing reference frames. The new reference-topology API must include
this correction and preserve its rotation-composition regression tests.

The incoming [0.16 API cleanup, PR #119](https://github.com/mmikhasenko/ThreeBodyDecays.jl/pull/119)
removes the `VertexFunction` alias. Use `Vertex` in new code. The transition must
also be tested together with that branch before the two core changes are merged.

## Model migration

Construct models with an explicit reference, for example
`ThreeBodyDecay(description; reference_k = 1)`. The old public `refζs` keyword is
removed. Slicing preserves the reference; concatenating models with different
references is an error. Copying a coefficient list into a differently defined
model is not an automatic convention conversion.

The normal polarized-event path also takes the model itself:

```julia
angles, invariants = angles_invariants(vectors, model)
A = amplitude(model, angles, invariants)
```

Neither call accepts a reference override. The separate low-level conversion
with masses and an explicit `k` is for geometry work, not for changing the
reference of an existing fitted model.

An old common reference `(k, k, k, k)` maps to `reference_k = k`. A mixed tuple
such as the old default `(1, 2, 3, 1)` has no literal representation in the new
API. Choose and document a common reference and revalidate the model, especially
its helicity components and polarization-dependent observables. An invariant
unpolarized checksum alone does not establish identical amplitude components.

## Document migration

The implementation follows the draft
[amplitude-serialization PR #111](https://github.com/RUB-EP1/amplitude-serialization/pull/111)
at commit `2d75724dde80ae94c97d2ce0ab2673a8da489d9b`. This is a pinned draft
contract, not a claim that the standard has already been finalized.

- Store masses only in `domains` entries whose type is `phase_space`.
- Store particle identity and spin in `decay_description.kinematics`.
- Identify the parent by index 0 and the three final particles by indices 1–3.
  Match all records by index, never by their position in a JSON array.
- Represent a point as `{name, domain, final_state: [{index, p: {E, px, py, pz}}]}`.
- Validate full distributions using `misc.amplitude_model_checksums`. Old scalar
  lineshape checks and angle-parameter points are not a fallback input format.
- Resolve the domain from each checksum point. Without a point or explicit
  domain selector, multiple phase-space domains are ambiguous and must not be
  selected by array position.

For serialized validation events, use the parent rest frame. The conversion API
must document any support for moving-parent input separately, including its
boost and spin-axis convention. Domain validation checks energies, on-shell
masses, and total four-momentum before evaluation.

## Integration and release order

1. Review the parent full-turn correction and the core API cleanup together with
   the reference-topology migration. Keep amplitude-component tests as well as
   coherent intensity and overlap tests.
2. Test the IO branch against the exact proposed core revision. Check the pinned
   standard's three reference models, shuffled index records, all three reference
   topologies, malformed domains/events, and nontrivial spatial orientations.
3. Merge and release the core version first. Confirm the registered/released
   artifact contains the required API rather than relying on its version number.
4. Replace the IO CI prerelease-core pin with ordinary package resolution and
   rerun the complete IO suite against the released core.
5. Release IO 0.10 with this migration guide. Update consumers of the old schema
   deliberately; do not silently reinterpret legacy points or fitted conventions.

No change to either default branch or release is implied by preparing the
coordinated development branches.
