# Optional standalone oracle; see README.md in this directory.
using ThreeBodyDecaysIO
using ThreeBodyDecays
import InstructionalDecayTrees
import FourVectors
using Random, Test

const FV = FourVectors
const IDT = InstructionalDecayTrees

@testset "InstructionalDecayTrees three-body orientation oracle" begin
    Random.seed!(3441)
    ms = ThreeBodyMasses(1, 2, 3; m0 = 7.9)
    for sample = 1:20, k = 1:3
        sigmas = randomPoint(ms)
        i, j = ij_from_k(k)
        # Instruction paths for ((i, j), k): measure the isobar in the
        # current parent frame, then daughter i in the isobar helicity frame.
        programs = (
            (IDT.MeasureCosThetaPhi(:root, (i, j)),),
            (IDT.ToHelicityFrame((i, j)), IDT.MeasureCosThetaPhi(:inner, (i,))),
        )
        aligned = map(FV.FourVector, aligned_four_vectors(sigmas, ms; k))
        event = map(aligned) do p
            p |> FV.Rz(-0.93) |> FV.Ry(0.62) |> FV.Rz(0.71)
        end
        angles, _ = angles_invariants(event, ms; k)
        measured = map(programs) do program
            _, result = IDT.apply_decay_instruction(program, event)
            only(values(result))
        end
        @test angles.α ≈ measured[1].ϕ atol = 1e-12
        @test angles.cosβ ≈ measured[1].cosθ atol = 1e-12
        @test angles.γ ≈ measured[2].ϕ atol = 1e-12
    end
end
