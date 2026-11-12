# ThreeBodyDecaysIO

## Overview

`ThreeBodyDecaysIO` is a Julia module designed for the serialization and deserialization of models related to three-body decays in particle physics. It leverages the `ThreeBodyDecays` framework to provide comprehensive support for reading and writing model descriptions, facilitating the exchange of complex decay models and their kinematic configurations.

## Features

-   Serialization and deserialization of three-body decay models to and from JSON format.
-   Support for detailed kinematics, lineshapes, and interaction chains descriptions.
-   Utilities for validation and manipulation of model components.

## Installation

To install `ThreeBodyDecaysIO`, use the Julia package manager. From the Julia REPL, type the following:

```julia
Pkg.add(url="https://github.com/RUB-EP1/ThreeBodyDecaysIO.jl")
```

## Usage

### Writing and reading a model

`serializeToDict(model)` returns a decay description and a function appendix.
Use `add_hs3_fields` to build the full document:

```julia
using ThreeBodyDecaysIO, ThreeBodyDecays, HadronicLineshapes, JSON

description, functions = serializeToDict(model)
document = add_hs3_fields(description, functions, "my_model")
open("model.json", "w") do io
    JSON.print(io, document, 4)
end

# In-memory round trip. The function appendix is the reader's workspace.
restored = dict2instance(ThreeBodyDecay, description; workspace = functions)
```

When reading JSON, deserialize each entry in `document["functions"]` with
`dict2instance` and a supported type such as `BreitWigner`, then pass the resulting
name-to-function dictionary as `workspace`. Chain names are preserved independently
of the propagator names used to reference that workspace.

### Functions with named arguments

`dict2instance` returns `NamedArgFunc` wrappers for Breit-Wigner, multichannel
Breit-Wigner, polynomial, and generic functions. Each wrapper keeps one argument
name and accepts a dictionary containing its value:

```julia
f = dict2instance(Polynomial, Dict("coefficients" => [1.0, 2.0], "x" => "s"))
f(Dict("s" => 3.0)) # 7.0
serialized, appendix = serializeToDict(f)
```

Breit-Wigner and polynomial wrappers support serialization back to dictionaries.
Generic expressions support reading and evaluation; arbitrary Julia closures
cannot be serialized.

### Cascade models

For CascadeDecays models, use the separate
[CascadeDecaysIO.jl](https://github.com/RUB-EP1/CascadeDecaysIO.jl) package,
which owns cascade serialization, reading, and validation. ThreeBodyDecaysIO
has no CascadeDecays dependency or extension.

Core functionality supports Julia 1.7, ThreeBodyDecays 0.13/0.15, and JSON 0.21/1.
Parameters 0.12 remains allowed because the supported ThreeBodyDecays releases
require it.

## Running Tests

To ensure `ThreeBodyDecaysIO` is working correctly, you can run its test suite:

```julia
using Pkg
Pkg.test("ThreeBodyDecaysIO")
```

This will execute a series of tests, verifying the functionality of model serialization/deserialization, kinematics parsing, and more. See details in [`test/runtests.jl`](test/runtests.jl)

## Running pre-commit hooks

To check code formatting, spelling, and other style issues before committing, you can run the pre-commit hooks manually:

```sh
pre-commit run --all-files
```

This will run all configured checks on the entire codebase and automatically fix some issues. See `.pre-commit-config.yaml` for details.

To skip pre-commit hooks for a specific commit (e.g., for quick commits), use the `--no-verify` flag:

```sh
git commit -m "your message" --no-verify
```

This allows you to commit without running the pre-commit checks when needed.

## Contributing

Contributions to `ThreeBodyDecaysIO` are welcome. To contribute, please fork the repository, make your changes, and submit a pull request.

## ThreeBodyDecays 0.16

Requires ThreeBodyDecays 0.16. Vertex serialization uses `Vertex`. Model serialization writes the model's stored `reference_k` as `reference_topology`, and deserialization restores that reference. A description must include `reference_topology`; it is never inferred from the first chain. For legacy fitted models, verify that the serialized topology matches the coupling convention before importing: the former mixed reference default has no general single-topology replacement.
