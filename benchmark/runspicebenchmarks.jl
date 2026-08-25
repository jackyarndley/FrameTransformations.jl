import Pkg

const BENCHMARK_DIR = @__DIR__
const REPO_ROOT = normpath(joinpath(BENCHMARK_DIR, ".."))
const TEST_ASSETS = joinpath(REPO_ROOT, "test", "assets")

Pkg.activate(BENCHMARK_DIR; io=devnull)
if VERSION < v"1.12"
    Pkg.develop(Pkg.PackageSpec(path=REPO_ROOT); io=devnull)
    Pkg.instantiate(; io=devnull)
else
    Pkg.resolve(; io=devnull, workspace=true)
    Pkg.instantiate(; io=devnull, workspace=true)
end

using BenchmarkTools
using Ephemerides: EphemerisProvider, ephem_vector3, ephem_vector6, prepare_ephemeris
using FrameTransformations
using Printf: @printf
using StaticArrays: SVector

# Mirrors the SPICE workloads added in duncaneddy/brahe#376:
# See benchmarks/spice_benchmarks.rs at commit 1204781560b2 in duncaneddy/brahe.

const ICRF = 1
const SSB = 0
const EARTH_BARYCENTER = 3
const MARS_BARYCENTER = 4
const SUN = 10
const MOON = 301
const EARTH = 399

# TDB seconds past J2000, matching the dates used by the Brahe benchmark. Epoch
# construction and UTC-to-TDB conversion are deliberately outside the timed workloads.
const SINGLE_EPOCH = 802_051_269.184 # 2025-06-01T00:00:00 UTC
const SEQUENTIAL_EPOCH = 788_961_669.184 # 2025-01-01T00:00:00 UTC
const SEQUENTIAL_STEP = 60.0

const GM_SUN = 1.32712440018e20 # m^3 / s^2
const GM_MOON = 4.902800118e12 # m^3 / s^2
const SATELLITE_POSITION = SVector(7.0e6, 0.0, 0.0) # m

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

function spice_kernel_path()
    if haskey(ENV, "SPICE_KERNEL")
        path = abspath(ENV["SPICE_KERNEL"])
        isfile(path) || throw(ArgumentError("SPICE_KERNEL does not exist: $path"))
        return path
    end

    de440s = joinpath(TEST_ASSETS, "de440s.bsp")
    isfile(de440s) && return de440s

    fallback = joinpath(TEST_ASSETS, "de432s.bsp")
    isfile(fallback) || throw(ArgumentError(
        "No SPK kernel found. Set SPICE_KERNEL to a de440s.bsp path."
    ))
    @warn "Using bundled DE432s. Set SPICE_KERNEL to de440s.bsp for an exact " *
          "Brahe PR #376 comparison."
    return fallback
end

function build_frame_system(ephemeris)
    frames = FrameSystem{2,Float64}()
    add_axes_icrf!(frames)
    add_point!(frames, :SSB, SSB, ICRF)
    add_point_ephemeris!(frames, ephemeris, :Sun, SUN)
    add_point_ephemeris!(frames, ephemeris, :EarthBarycenter, EARTH_BARYCENTER)
    add_point_ephemeris!(frames, ephemeris, :Earth, EARTH)
    add_point_ephemeris!(frames, ephemeris, :Moon, MOON)
    add_point_ephemeris!(frames, ephemeris, :MarsBarycenter, MARS_BARYCENTER)
    return frames
end

function prepared_chain(ephemeris, target)
    earth_from_earth_barycenter = prepare_ephemeris(
        ephemeris, EARTH_BARYCENTER, EARTH
    )

    if target == MOON
        target_from_earth_barycenter = prepare_ephemeris(
            ephemeris, EARTH_BARYCENTER, target
        )
        position = let target_route = target_from_earth_barycenter,
            earth_route = earth_from_earth_barycenter
            t -> ephem_vector3(target_route, t) - ephem_vector3(earth_route, t)
        end
        state = let target_route = target_from_earth_barycenter,
            earth_route = earth_from_earth_barycenter
            t -> ephem_vector6(target_route, t) - ephem_vector6(earth_route, t)
        end
    else
        earth_barycenter_from_ssb = prepare_ephemeris(
            ephemeris, SSB, EARTH_BARYCENTER
        )
        target_from_ssb = prepare_ephemeris(ephemeris, SSB, target)
        position = let target_route = target_from_ssb,
            earth_barycenter_route = earth_barycenter_from_ssb,
            earth_route = earth_from_earth_barycenter
            t -> ephem_vector3(target_route, t) -
                 ephem_vector3(earth_barycenter_route, t) -
                 ephem_vector3(earth_route, t)
        end
        state = let target_route = target_from_ssb,
            earth_barycenter_route = earth_barycenter_from_ssb,
            earth_route = earth_from_earth_barycenter
            t -> ephem_vector6(target_route, t) -
                 ephem_vector6(earth_barycenter_route, t) -
                 ephem_vector6(earth_route, t)
        end
    end

    return (; position, state)
end

function body_case(frames, ephemeris, name, target)
    prepared = prepared_chain(ephemeris, target)
    direct_position = t -> vector3(frames, EARTH, target, ICRF, t)
    direct_state = t -> vector6(frames, EARTH, target, ICRF, t)
    compiled_position = compile_translation(frames, EARTH, target, ICRF, Val(1))
    compiled_state = compile_translation(frames, EARTH, target, ICRF, Val(2))

    reference_position = prepared.position(SINGLE_EPOCH)
    reference_state = prepared.state(SINGLE_EPOCH)
    @assert direct_position(SINGLE_EPOCH) ≈ reference_position
    @assert compiled_position(SINGLE_EPOCH) ≈ reference_position
    @assert direct_state(SINGLE_EPOCH) ≈ reference_state
    @assert compiled_state(SINGLE_EPOCH) ≈ reference_state

    return (;
        name,
        prepared_position=prepared.position,
        prepared_state=prepared.state,
        direct_position,
        direct_state,
        compiled_position,
        compiled_state,
    )
end

function build_cases()
    kernel = spice_kernel_path()
    ephemeris = EphemerisProvider(kernel)
    frames = build_frame_system(ephemeris)
    cases = (
        body_case(frames, ephemeris, "sun", SUN),
        body_case(frames, ephemeris, "moon", MOON),
        body_case(frames, ephemeris, "mars", MARS_BARYCENTER),
    )
    return (; kernel, frames, ephemeris, cases)
end

function add_single_query!(suite, case)
    body = BenchmarkGroup()
    body["position"] = BenchmarkGroup()
    body["state"] = BenchmarkGroup()

    for path in ("prepared", "direct", "compiled")
        position = getproperty(case, Symbol(path, "_position"))
        state = getproperty(case, Symbol(path, "_state"))
        body["position"][path] = @benchmarkable $position($SINGLE_EPOCH)
        body["state"][path] = @benchmarkable $state($SINGLE_EPOCH)
    end

    suite[case.name] = body
    return suite
end

function sequential_position(function_object, initial_epoch, count, step)
    accumulator = SVector(0.0, 0.0, 0.0)
    for index in 0:(count - 1)
        accumulator += function_object(initial_epoch + index * step)
    end
    return accumulator
end

function sequential_state(function_object, initial_epoch, count, step)
    accumulator = 0.0
    for index in 0:(count - 1)
        state = function_object(initial_epoch + index * step)
        accumulator += sqrt(sum(abs2, state))
    end
    return accumulator
end

@inline function third_body_acceleration(body_position_km, satellite_position_m, gm)
    body_position_m = 1_000.0 * body_position_km
    relative_position = body_position_m - satellite_position_m
    relative_scale = sum(abs2, relative_position)^(-1.5)
    body_scale = sum(abs2, body_position_m)^(-1.5)
    return gm * (relative_scale * relative_position - body_scale * body_position_m)
end

function sun_moon_acceleration(sun_position, moon_position, epoch)
    sun_acceleration = third_body_acceleration(
        sun_position(epoch), SATELLITE_POSITION, GM_SUN
    )
    moon_acceleration = third_body_acceleration(
        moon_position(epoch), SATELLITE_POSITION, GM_MOON
    )
    return sun_acceleration, moon_acceleration
end

function build_suite(setup, sequential_count)
    suite = BenchmarkGroup()
    suite["single_query"] = BenchmarkGroup()
    for case in setup.cases
        add_single_query!(suite["single_query"], case)
    end

    sun = setup.cases[1]
    moon = setup.cases[2]
    suite["sequential"] = BenchmarkGroup()
    suite["third_body"] = BenchmarkGroup()

    for path in ("prepared", "direct", "compiled")
        sun_position = getproperty(sun, Symbol(path, "_position"))
        sun_state = getproperty(sun, Symbol(path, "_state"))
        moon_position = getproperty(moon, Symbol(path, "_position"))

        group = BenchmarkGroup()
        group["position"] = @benchmarkable sequential_position(
            $sun_position, $SEQUENTIAL_EPOCH, $sequential_count, $SEQUENTIAL_STEP
        )
        group["state"] = @benchmarkable sequential_state(
            $sun_state, $SEQUENTIAL_EPOCH, $sequential_count, $SEQUENTIAL_STEP
        )
        suite["sequential"][path] = group

        suite["third_body"][path] = @benchmarkable sun_moon_acceleration(
            $sun_position, $moon_position, $SINGLE_EPOCH
        )
    end

    return suite
end

trial_summary(trial) = BenchmarkTools.median(trial)

function print_single_query_summary(results)
    println()
    println("Single-query median timing (ns, bytes)")
    @printf(
        "%-6s %-8s %12s %12s %12s %10s %10s %10s\n",
        "body", "query", "prepared", "direct", "compiled", "prep alloc", "dir alloc",
        "comp alloc",
    )
    println(repeat("-", 94))

    for body in ("sun", "moon", "mars"), query in ("position", "state")
        group = results["single_query"][body][query]
        prepared = trial_summary(group["prepared"])
        direct = trial_summary(group["direct"])
        compiled = trial_summary(group["compiled"])
        @printf(
            "%-6s %-8s %12.1f %12.1f %12.1f %10d %10d %10d\n",
            body, query, prepared.time, direct.time, compiled.time,
            prepared.memory, direct.memory, compiled.memory,
        )
    end
end

function print_workload_summary(results, sequential_count)
    println()
    println("Sequential Sun-relative-to-Earth median timing ($sequential_count epochs, ms)")
    @printf("%-10s %14s %14s\n", "path", "position", "state")
    println(repeat("-", 40))
    for path in ("prepared", "direct", "compiled")
        position = trial_summary(results["sequential"][path]["position"])
        state = trial_summary(results["sequential"][path]["state"])
        @printf("%-10s %14.3f %14.3f\n", path, position.time / 1e6, state.time / 1e6)
    end

    println()
    println("Sun + Moon third-body adapter median timing (ns, bytes)")
    @printf("%-10s %14s %14s\n", "path", "time", "allocation")
    println(repeat("-", 40))
    for path in ("prepared", "direct", "compiled")
        trial = trial_summary(results["third_body"][path])
        @printf("%-10s %14.1f %14d\n", path, trial.time, trial.memory)
    end
end

function main()
    configure_benchmarktools!()
    sequential_count = benchmark_int("SPICE_SEQUENTIAL_EPOCHS", 10_000)
    sequential_count > 0 ||
        throw(ArgumentError("SPICE_SEQUENTIAL_EPOCHS must be positive."))

    setup = build_cases()
    suite = build_suite(setup, sequential_count)

    println("Brahe PR #376-compatible SPICE benchmark")
    println("Kernel: ", setup.kernel)
    println("Sequential epochs: ", sequential_count, " at ", SEQUENTIAL_STEP, " s steps")
    println("Running benchmark suite...")
    results = run(suite; verbose=true)

    print_single_query_summary(results)
    print_workload_summary(results, sequential_count)
end

main()
