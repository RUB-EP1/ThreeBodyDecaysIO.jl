# Optional InstructionalDecayTrees oracle

`check_instructional_decay_trees.jl` compares IO orientation extraction with
InstructionalDecayTrees execution of the two instruction paths for `((i,j),k)`:
measure the isobar in the parent frame, then daughter `i` in the isobar helicity
frame. It checks 20 deterministic physical events for each reference topology:
180 angle assertions. It compares root azimuth/polar cosine and inner daughter
azimuth. It does not compare amplitude tensors or fitted coupling conventions.

This standalone cross-check is excluded from `runtests.jl`.
InstructionalDecayTrees is not an IO runtime or test dependency. CascadeDecays
is not needed by this oracle.

Tested with Julia 1.11.5, FourVectors 0.3.4 and InstructionalDecayTrees 0.4.1
(tree `3f64d77c542dadce88241d72bf355cb07ea2ed64`).

Run in an isolated environment, replacing the two local checkout placeholders:

```julia
using Pkg
Pkg.activate(mktempdir())
Pkg.develop([
    PackageSpec(path = "/path/to/ThreeBodyDecays.jl"),
    PackageSpec(path = "/path/to/ThreeBodyDecaysIO.jl"),
])
Pkg.add([
    PackageSpec(name = "InstructionalDecayTrees", version = "0.4.1"),
    PackageSpec(name = "FourVectors", version = "0.3.4"),
])
include("/path/to/ThreeBodyDecaysIO.jl/test/oracles/check_instructional_decay_trees.jl")
```

Expected result: 180 passing assertions. No package compatibility changes are
needed. Keep the temporary environment's Manifest for exact reproduction.
