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
using DifferentiationInterface: AutoForwardDiff, derivative
import ForwardDiff
using Printf: @printf
using StaticArrays: SVector

include("scenarios.jl")

const backend = AutoForwardDiff()

function benchmark_int(name, default)
    value = tryparse(Int, get(ENV, name, string(default)))
    isnothing(value) && throw(ArgumentError("$name must parse as Int."))
    return value
end

function benchmark_float(name, default)
    value = tryparse(Float64, get(ENV, name, string(default)))
    isnothing(value) && throw(ArgumentError("$name must parse as Float64."))
    return value
end

function configure_benchmarktools!()
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = benchmark_float("BENCHMARK_SECONDS", 1.0)
    BenchmarkTools.DEFAULT_PARAMETERS.samples = benchmark_int("BENCHMARK_SAMPLES", 10)
    return nothing
end

function add_case!(cases, name, direct, prepared, compiled)
    push!(cases, (; name, direct, prepared, compiled))
    return cases
end

function benchmark_cases()
    cases = NamedTuple[]

    for order in (2, 3, 4)
        synthetic = build_synthetic_benchmark_scenario(order)

        synthetic_rotation_position(t) = synthetic.rotation_direct(t)[1]
        synthetic_rotation_position_compiled(t) = synthetic.rotation_compiled(t)[1]
        synthetic_rotation_position_prepared(t) = synthetic.rotation_prepared(t)[1]
        synthetic_translation_position(t) = synthetic.translation_direct(t)[SVector(1, 2, 3)]
        synthetic_translation_position_compiled(t) = synthetic.translation_compiled(t)[SVector(1, 2, 3)]
        synthetic_translation_position_prepared(t) = synthetic.translation_prepared(t)[SVector(1, 2, 3)]

        one_hop_rotation_direct = rotation_direct_callable(
            synthetic.fr, :ICRF, :A, order)
        one_hop_rotation_prepared = prepare_rotation(
            synthetic.fr, :ICRF, :A, Val(order))
        one_hop_rotation_compiled = compile_rotation(
            synthetic.fr, :ICRF, :A, Val(order))
        one_hop_translation_direct = translation_direct_callable(
            synthetic.fr, :Origin, :P1, :ICRF, order)
        one_hop_translation_prepared = prepare_translation(
            synthetic.fr, :Origin, :P1, :ICRF, Val(order))
        one_hop_translation_compiled = compile_translation(
            synthetic.fr, :Origin, :P1, :ICRF, Val(order))

        add_case!(
            cases, "synthetic O$order rot one",
            () -> one_hop_rotation_direct(synthetic.t0),
            () -> one_hop_rotation_prepared(synthetic.t0),
            () -> one_hop_rotation_compiled(synthetic.t0),
        )
        add_case!(
            cases, "synthetic O$order vec one",
            () -> one_hop_translation_direct(synthetic.t0),
            () -> one_hop_translation_prepared(synthetic.t0),
            () -> one_hop_translation_compiled(synthetic.t0),
        )
        add_case!(
            cases, "synthetic O$order rot$(synthetic.suffix)",
            () -> synthetic.rotation_direct(synthetic.t0),
            () -> synthetic.rotation_prepared(synthetic.t0),
            () -> synthetic.rotation_compiled(synthetic.t0),
        )
        add_case!(
            cases, "synthetic O$order vec$(synthetic.suffix)",
            () -> synthetic.translation_direct(synthetic.t0),
            () -> synthetic.translation_prepared(synthetic.t0),
            () -> synthetic.translation_compiled(synthetic.t0),
        )
        add_case!(
            cases, "synthetic O$order dir$(synthetic.suffix)",
            () -> synthetic.direction_direct(synthetic.t0),
            () -> synthetic.direction_prepared(synthetic.t0),
            () -> synthetic.direction_compiled(synthetic.t0),
        )
        add_case!(
            cases, "synthetic O$order d/dt rot",
            () -> derivative(synthetic_rotation_position, backend, synthetic.t0),
            () -> derivative(synthetic_rotation_position_prepared, backend, synthetic.t0),
            () -> derivative(synthetic_rotation_position_compiled, backend, synthetic.t0),
        )
        add_case!(
            cases, "synthetic O$order d/dt vec",
            () -> derivative(synthetic_translation_position, backend, synthetic.t0),
            () -> derivative(synthetic_translation_position_prepared, backend, synthetic.t0),
            () -> derivative(synthetic_translation_position_compiled, backend, synthetic.t0),
        )
    end

    for order in (2, 3, 4)
        de440 = build_de440_benchmark_scenario(order)
        de440_vector = build_de440_vector_benchmark_scenario(order)

        de440_rotation_position(t) = de440.rotation_direct(t)[1]
        de440_rotation_position_compiled(t) = de440.rotation_compiled(t)[1]
        de440_rotation_position_prepared(t) = de440.rotation_prepared(t)[1]
        de440_translation_position(t) = de440_vector.translation_direct(t)[SVector(1, 2, 3)]
        de440_translation_position_compiled(t) = de440_vector.translation_compiled(t)[SVector(1, 2, 3)]
        de440_translation_position_prepared(t) = de440_vector.translation_prepared(t)[SVector(1, 2, 3)]

        add_case!(
            cases, "de440 O$order rot$(de440.suffix)",
            () -> de440.rotation_direct(de440.t0),
            () -> de440.rotation_prepared(de440.t0),
            () -> de440.rotation_compiled(de440.t0),
        )
        add_case!(
            cases, "de440 O$order vec$(de440_vector.suffix)",
            () -> de440_vector.translation_direct(de440_vector.t0),
            () -> de440_vector.translation_prepared(de440_vector.t0),
            () -> de440_vector.translation_compiled(de440_vector.t0),
        )
        add_case!(
            cases, "de440 O$order d/dt rot",
            () -> derivative(de440_rotation_position, backend, de440.t0),
            () -> derivative(de440_rotation_position_prepared, backend, de440.t0),
            () -> derivative(de440_rotation_position_compiled, backend, de440.t0),
        )
        add_case!(
            cases, "de440 O$order d/dt vec",
            () -> derivative(de440_translation_position, backend, de440_vector.t0),
            () -> derivative(de440_translation_position_prepared, backend, de440_vector.t0),
            () -> derivative(de440_translation_position_compiled, backend, de440_vector.t0),
        )
    end

    return Tuple(cases)
end

function add_construction_suite!(suite)
    scenario = build_synthetic_benchmark_scenario(2)
    frames = scenario.fr
    construction = BenchmarkGroup()
    construction["prepare rotation"] = @benchmarkable prepare_rotation(
        $frames, :ICRF, :C, Val(2))
    construction["compile rotation"] = @benchmarkable compile_rotation(
        $frames, :ICRF, :C, Val(2))
    construction["prepare translation"] = @benchmarkable prepare_translation(
        $frames, :Origin, :P2, :A, Val(2))
    construction["compile translation"] = @benchmarkable compile_translation(
        $frames, :Origin, :P2, :A, Val(2))
    suite["construction"] = construction
    return suite
end

function build_suite(cases)
    suite = BenchmarkGroup()

    for case in cases
        group = BenchmarkGroup()
        group["direct"] = @benchmarkable $(case.direct)()
        group["prepared"] = @benchmarkable $(case.prepared)()
        group["compiled"] = @benchmarkable $(case.compiled)()
        suite[case.name] = group
    end

    return suite
end

trial_summary(trial) = BenchmarkTools.median(trial)

function print_summary(results, cases)
    println()
    println("Median steady-state timing (ns)")
    @printf("%-24s %14s %14s %14s\n", "benchmark", "direct", "prepared", "compiled")
    println(repeat("-", 72))

    for case in cases
        direct = trial_summary(results[case.name]["direct"])
        prepared = trial_summary(results[case.name]["prepared"])
        compiled = trial_summary(results[case.name]["compiled"])

        @printf(
            "%-24s %14.1f %14.1f %14.1f\n",
            case.name,
            direct.time,
            prepared.time,
            compiled.time,
        )
    end


    println()
    println("Median steady-state allocations (count, bytes)")
    @printf(
        "%-24s %10s %10s %10s %10s %10s %10s\n",
        "benchmark", "dir count", "dir bytes", "prep count", "prep bytes",
        "comp count", "comp bytes",
    )
    println(repeat("-", 92))
    for case in cases
        direct = trial_summary(results[case.name]["direct"])
        prepared = trial_summary(results[case.name]["prepared"])
        compiled = trial_summary(results[case.name]["compiled"])
        @printf(
            "%-24s %10d %10d %10d %10d %10d %10d\n",
            case.name, direct.allocs, direct.memory, prepared.allocs,
            prepared.memory, compiled.allocs, compiled.memory,
        )
    end
end

function print_construction_summary(results)
    println()
    println("Median route construction (ns, count, bytes)")
    @printf("%-22s %14s %12s %12s\n", "operation", "time", "allocations", "bytes")
    println(repeat("-", 64))
    for name in (
        "prepare rotation", "compile rotation",
        "prepare translation", "compile translation",
    )
        estimate = trial_summary(results["construction"][name])
        @printf(
            "%-22s %14.1f %12d %12d\n",
            name, estimate.time, estimate.allocs, estimate.memory,
        )
    end
end

function main()
    configure_benchmarktools!()

    cases = benchmark_cases()
    suite = build_suite(cases)
    add_construction_suite!(suite)

    println("Running benchmark suite...")
    results = run(suite; verbose = true)

    print_summary(results, cases)
    print_construction_summary(results)
end

main()
