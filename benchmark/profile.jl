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

using ForwardDiff
using Profile
using StaticArrays: SVector

include("scenarios.jl")

function benchmark_int(name, default)
    value = tryparse(Int, get(ENV, name, string(default)))
    isnothing(value) && throw(ArgumentError("$name must parse as Int."))
    return value
end

profile_iterations() = benchmark_int("PROFILE_ITERATIONS", 50_000)
profile_ad_iterations() = benchmark_int("PROFILE_AD_ITERATIONS", 20_000)

function profile_block(name, f)
    Profile.clear()
    @profile f()
    println()
    println("== ", name, " ==")
    Profile.print(format = :flat, sortedby = :count, maxdepth = 16)
end

function main()
    synthetic = build_synthetic_benchmark_scenario()
    de440 = build_de440_benchmark_scenario()
    n = profile_iterations()
    n_ad = profile_ad_iterations()

    synthetic.translation_direct(synthetic.t0)
    synthetic.translation_compiled(synthetic.t0)
    de440.rotation_direct(de440.t0)
    de440.rotation_compiled(de440.t0)
    ForwardDiff.derivative(t -> synthetic.translation_compiled(t)[SVector(1, 2, 3)], synthetic.t0)

    Profile.init(delay = 0.0001)

    profile_block("synthetic direct vector12", () -> begin
        for _ in 1:n
            synthetic.translation_direct(synthetic.t0)
        end
    end)

    profile_block("synthetic compiled vector12", () -> begin
        for _ in 1:n
            synthetic.translation_compiled(synthetic.t0)
        end
    end)

    profile_block("de440 direct rotation12", () -> begin
        for _ in 1:n
            de440.rotation_direct(de440.t0)
        end
    end)

    profile_block("de440 compiled rotation12", () -> begin
        for _ in 1:n
            de440.rotation_compiled(de440.t0)
        end
    end)

    profile_block("synthetic compiled ForwardDiff", () -> begin
        pos(t) = synthetic.translation_compiled(t)[SVector(1, 2, 3)]
        for _ in 1:n_ad
            ForwardDiff.derivative(pos, synthetic.t0)
        end
    end)
end

main()
