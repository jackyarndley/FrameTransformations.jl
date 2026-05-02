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
using ForwardDiff
using Printf: @printf
using StaticArrays: SVector

include("scenarios.jl")

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

function add_case!(cases, name, direct, compiled)
    push!(cases, (name = name, direct = direct, compiled = compiled))
    return cases
end

function benchmark_cases()
    cases = NamedTuple[]

    for order in (2, 3, 4)
        synthetic = build_synthetic_benchmark_scenario(order)

        synthetic_rotation_position(t) = synthetic.rotation_direct(t)[1]
        synthetic_rotation_position_compiled(t) = synthetic.rotation_compiled(t)[1]
        synthetic_translation_position(t) = synthetic.translation_direct(t)[SVector(1, 2, 3)]
        synthetic_translation_position_compiled(t) = synthetic.translation_compiled(t)[SVector(1, 2, 3)]

        add_case!(cases, "synthetic O$order rot$(synthetic.suffix)", () -> synthetic.rotation_direct(synthetic.t0), () -> synthetic.rotation_compiled(synthetic.t0))
        add_case!(cases, "synthetic O$order vec$(synthetic.suffix)", () -> synthetic.translation_direct(synthetic.t0), () -> synthetic.translation_compiled(synthetic.t0))
        add_case!(cases, "synthetic O$order dir$(synthetic.suffix)", () -> synthetic.direction_direct(synthetic.t0), () -> synthetic.direction_compiled(synthetic.t0))
        add_case!(cases, "synthetic O$order d/dt rot", () -> ForwardDiff.derivative(synthetic_rotation_position, synthetic.t0), () -> ForwardDiff.derivative(synthetic_rotation_position_compiled, synthetic.t0))
        add_case!(cases, "synthetic O$order d/dt vec", () -> ForwardDiff.derivative(synthetic_translation_position, synthetic.t0), () -> ForwardDiff.derivative(synthetic_translation_position_compiled, synthetic.t0))
    end

    for order in (2, 3, 4)
        de440 = build_de440_benchmark_scenario(order)
        de440_vector = build_de440_vector_benchmark_scenario(order)

        de440_rotation_position(t) = de440.rotation_direct(t)[1]
        de440_rotation_position_compiled(t) = de440.rotation_compiled(t)[1]
        de440_translation_position(t) = de440_vector.translation_direct(t)[SVector(1, 2, 3)]
        de440_translation_position_compiled(t) = de440_vector.translation_compiled(t)[SVector(1, 2, 3)]

        add_case!(cases, "de440 O$order rot$(de440.suffix)", () -> de440.rotation_direct(de440.t0), () -> de440.rotation_compiled(de440.t0))
        add_case!(cases, "de440 O$order vec$(de440_vector.suffix)", () -> de440_vector.translation_direct(de440_vector.t0), () -> de440_vector.translation_compiled(de440_vector.t0))
        add_case!(cases, "de440 O$order d/dt rot", () -> ForwardDiff.derivative(de440_rotation_position, de440.t0), () -> ForwardDiff.derivative(de440_rotation_position_compiled, de440.t0))
        add_case!(cases, "de440 O$order d/dt vec", () -> ForwardDiff.derivative(de440_translation_position, de440_vector.t0), () -> ForwardDiff.derivative(de440_translation_position_compiled, de440_vector.t0))
    end

    return Tuple(cases)
end

function build_suite(cases)
    suite = BenchmarkGroup()

    for case in cases
        group = BenchmarkGroup()
        group["direct"] = @benchmarkable $(case.direct)()
        group["compiled"] = @benchmarkable $(case.compiled)()
        suite[case.name] = group
    end

    return suite
end

trial_summary(trial) = BenchmarkTools.median(trial)

function print_summary(results, cases)
    println()
    println("Median timing summary (ns, bytes)")
    @printf("%-24s %14s %14s %10s %14s %14s\n", "benchmark", "direct", "compiled", "speedup", "direct alloc", "compiled alloc")
    println(repeat("-", 98))

    for case in cases
        direct = trial_summary(results[case.name]["direct"])
        compiled = trial_summary(results[case.name]["compiled"])
        speedup = direct.time / compiled.time

        @printf(
            "%-24s %14.1f %14.1f %9.2fx %14d %14d\n",
            case.name,
            direct.time,
            compiled.time,
            speedup,
            direct.memory,
            compiled.memory,
        )
    end
end

function main()
    configure_benchmarktools!()

    cases = benchmark_cases()
    suite = build_suite(cases)

    println("Running benchmark suite...")
    results = run(suite; verbose = true)

    print_summary(results, cases)
end

main()
