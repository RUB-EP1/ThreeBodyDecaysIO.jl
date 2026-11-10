# ThreeBodyDecaysIO

Read and write three-body amplitude models and evaluate indexed four-vector
phase-space points using ThreeBodyDecays.jl and FourVectors.jl.

Version 0.10 is a breaking transition requiring Julia 1.10 and ThreeBodyDecays 0.16.
It implements the draft phase-space layout in
[amplitude-serialization PR111](https://github.com/RUB-EP1/amplitude-serialization/pull/111),
pinned in tests to `2d75724dde80ae94c97d2ce0ab2673a8da489d9b`.
See the [coordinated migration guide](docs/reference-topology-transition.md).

## Read, evaluate, write, validate

```julia
using ThreeBodyDecaysIO, ThreeBodyDecays, JSON

document = JSON.parsefile("model.json")
distribution = read_model(document; distribution="default_model", domain="default")
point = first(document["parameter_points"])
intensity = distribution(point)

# Polarized amplitudes use the same stored reference topology.
p = four_vectors(point, distribution.domain)
angles, invariants = angles_invariants(p, distribution.model)
A = amplitude(distribution.model, angles, invariants)

# The writer obtains reference_topology from the model, including when the
# first chain has a different topology. Masses are written only in domains.
output = add_hs3_fields(distribution.model, "default_model"; domain_name="default")
merge!(output, validation_section(distribution.model, [p];
    point_names=["check"], model_name="default_model", domain_name="default"))
open("roundtrip.json", "w") do io
    JSON.print(io, output, 4)
end
results = validate(output)
@assert all(r -> r.passed, results)
```

Omit selectors only if exactly one distribution and one phase-space domain exist.
Validation selects each checksum's distribution and its point's named domain.
Each result contains `point`, `distribution`, `value`, `expected`, `passed`, and
`error`. A malformed checksum produces a failed result and later checks continue.
Missing or empty checksums raise an error. Default checksum tolerances are
`rtol=1e-10`, `atol=1e-12`.

`read_model` returns `HadronicUnpolarizedIntensity` with `.model` and `.domain`.
The reference topology is chosen when constructing the core model using
`ThreeBodyDecay(pairs; reference_k=k)` and is inseparable from its couplings.
Evaluation has no reference-topology override.

## Phase-space conventions

A `phase_space` domain holds parent index 0 and final-state indices 1, 2, 3 with
masses. Distribution kinematics hold matching identities and spins, without
masses. A point holds `name`, `domain`, and `final_state`, whose entries have
`index` and `p = {E, px, py, pz}`. Arrays may be reordered; indices determine
particle identity. Duplicate, missing, unexpected, or noninteger indices fail.
Legacy scalar angle points and flat `four_vectors` objects are not supported.

Standard points must be in the parent rest frame. Positive finite energies,
on-shell masses, total mass, and momentum conservation are checked. Absolute
tolerances are `1e-6` for squared masses and `1e-8` for momentum/energy closure;
units are those of the input domain. These are validation tolerances, not a
projection onto physical phase space.

`angles_invariants(vectors, masses; k)` accepts three FourVectors in index order.
For `[[i,j],k]`, the aligned isobar points along +z, spectator k along -z, and
daughter i lies in the positive-x half-plane. The active Euler rotation is
`Rz(α) Ry(β) Rz(γ)`. At either polar axis α is chosen as zero; a collinear daughter
uses γ=0. Principal `atan`/`acos` values define the Euler chart. Four-vectors
fix an SO(3) orientation, not a unique SU(2) lift; polarized components use this
declared chart and the core amplitude convention, without extra spectator phases.
An isobar at rest has undefined orientation and is rejected.

For explicit low-level conversion of a moving parent, pass `frame=:canonical`.
This applies a rotation-free boost into the parent rest frame. Polarization axes
are those of this canonical frame; it does not infer a parent helicity frame.
The standard distribution reader continues to require rest-frame points.

## Functions and serialization

Supported function types are BreitWigner, MultichannelBreitWigner, Polynomial,
BlattWeisskopf, MomentumPower, and generic_function. Generic expressions accept
numeric analytic operations, one independent variable, and finite numeric
parameters. They do not execute arbitrary Julia code. Supply `custom_functions`
as a name-to-`NamedArgFunc` dictionary for an explicit custom implementation.
Arbitrary closures cannot be serialized.

`serializeToDict(model)` returns a decay description and function appendix.
`add_hs3_fields(description, functions, name; domain=PhaseSpaceDomain(...))`
assembles a document when using that lower-level interface. Chain names are
preserved independently of shared propagator names.

Cascade models remain in the separate
[CascadeDecaysIO.jl](https://github.com/RUB-EP1/CascadeDecaysIO.jl) package.
InstructionalDecayTrees.jl is an independent convention cross-check, not a dependency.
