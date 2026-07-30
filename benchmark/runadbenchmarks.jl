import Pkg

const BENCHMARK_DIR = @__DIR__
const REPO_ROOT = normpath(joinpath(BENCHMARK_DIR, ".."))

Pkg.activate(BENCHMARK_DIR; io=devnull)
if VERSION < v"1.12"
    Pkg.develop(Pkg.PackageSpec(path=REPO_ROOT); io=devnull)
    Pkg.instantiate(; io=devnull)
else
    Pkg.resolve(; io=devnull, workspace=true)
    Pkg.instantiate(; io=devnull, workspace=true)
end

using BenchmarkTools
using ChainRulesCore: rrule
using DifferentiationInterface:
    AutoFiniteDiff,
    AutoForwardDiff,
    AutoZygote,
    derivative
import FiniteDiff
import ForwardDiff
using FrameTransformations
using Printf: @printf
using StaticArrays: SVector
import Zygote

include("scenarios.jl")

function build_ad_suite()
    synthetic = build_synthetic_benchmark_scenario(4)
    ephemeris = build_de440_vector_benchmark_scenario(4)

    synthetic_position =
        time -> vector3(synthetic.fr, :Origin, :P2, :A, time)
    synthetic_position_compiled =
        compile_translation(synthetic.fr, :Origin, :P2, :A, Val(1))
    ephemeris_position =
        time -> vector3(ephemeris.fr, :Earth, :Moon, :ME421, time)
    ephemeris_position_compiled =
        compile_translation(
            ephemeris.fr, :Earth, :Moon, :ME421, Val(1)
        )
    synthetic_scalar = time -> sum(synthetic_position(time))
    ephemeris_scalar = time -> sum(ephemeris_position(time))

    forward_backend = AutoForwardDiff()
    finite_backend = AutoFiniteDiff(fdtype=Val(:central))
    reverse_backend = AutoZygote()

    _, synthetic_pullback = rrule(
        vector3,
        synthetic.fr,
        :Origin,
        :P2,
        :A,
        synthetic.t0,
    )
    output_tangent = Ref(SVector(1.0, 1.0, 1.0))

    suite = BenchmarkGroup()
    suite["synthetic"] = BenchmarkGroup()
    suite["ephemeris"] = BenchmarkGroup()
    suite["chainrules"] = BenchmarkGroup()

    suite["synthetic"]["analytic velocity"] = @benchmarkable vector6(
        $(synthetic.fr), :Origin, :P2, :A, $(synthetic.t0)
    )[4:6]
    suite["synthetic"]["forward direct"] = @benchmarkable derivative(
        $synthetic_position, $forward_backend, $(synthetic.t0)
    )
    suite["synthetic"]["forward compiled"] = @benchmarkable derivative(
        $synthetic_position_compiled, $forward_backend, $(synthetic.t0)
    )
    suite["synthetic"]["finite direct"] = @benchmarkable derivative(
        $synthetic_position, $finite_backend, $(synthetic.t0)
    )
    suite["synthetic"]["reverse scalar direct"] = @benchmarkable derivative(
        $synthetic_scalar, $reverse_backend, $(synthetic.t0)
    )
    suite["ephemeris"]["analytic velocity"] = @benchmarkable vector6(
        $(ephemeris.fr), :Earth, :Moon, :ME421, $(ephemeris.t0)
    )[4:6]
    suite["ephemeris"]["forward direct"] = @benchmarkable derivative(
        $ephemeris_position, $forward_backend, $(ephemeris.t0)
    )
    suite["ephemeris"]["forward compiled"] = @benchmarkable derivative(
        $ephemeris_position_compiled, $forward_backend, $(ephemeris.t0)
    )
    suite["ephemeris"]["reverse scalar direct"] = @benchmarkable derivative(
        $ephemeris_scalar, $reverse_backend, $(ephemeris.t0)
    )
    suite["chainrules"]["rrule construction"] = @benchmarkable rrule(
        vector3,
        $(synthetic.fr),
        :Origin,
        :P2,
        :A,
        $(synthetic.t0),
    )
    suite["chainrules"]["reused pullback"] =
        @benchmarkable $synthetic_pullback($output_tangent[])

    return suite
end

function print_ad_summary(group::BenchmarkGroup, prefix::String="")
    for key in sort!(collect(keys(group)); by=string)
        label = isempty(prefix) ? string(key) : "$prefix / $key"
        value = group[key]
        if value isa BenchmarkGroup
            print_ad_summary(value, label)
        else
            estimate = median(value)
            @printf(
                "%-42s %12.1f ns %8d B %4d alloc\n",
                label,
                estimate.time,
                estimate.memory,
                estimate.allocs,
            )
        end
    end
end

function main()
    BenchmarkTools.DEFAULT_PARAMETERS.seconds =
        parse(Float64, get(ENV, "BENCHMARK_SECONDS", "1.0"))
    BenchmarkTools.DEFAULT_PARAMETERS.samples =
        parse(Int, get(ENV, "BENCHMARK_SAMPLES", "10"))
    results = run(build_ad_suite(); verbose=true)
    print_ad_summary(results)
    return results
end

main()
