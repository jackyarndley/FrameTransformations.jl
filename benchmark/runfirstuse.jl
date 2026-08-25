import Pkg

const BENCHMARK_DIR = @__DIR__
const REPO_ROOT = normpath(joinpath(BENCHMARK_DIR, ".."))

function activate_benchmark()
    Pkg.activate(BENCHMARK_DIR; io=devnull)
    if VERSION < v"1.12"
        Pkg.develop(Pkg.PackageSpec(path=REPO_ROOT); io=devnull)
        Pkg.instantiate(; io=devnull)
    else
        Pkg.resolve(; io=devnull, workspace=true)
        Pkg.instantiate(; io=devnull, workspace=true)
    end
end

function measure(function_object, label)
    GC.gc()
    result = @timed function_object()
    println(
        rpad(label, 36),
        lpad(round(result.time * 1e3; digits=3), 12), " ms  ",
        lpad(result.bytes, 12), " bytes",
    )
    return result.value
end

function build_frames()
    frame_system = FrameSystem{2,Float64}()
    add_axes!(frame_system, :I, 1)
    add_axes_rotating!(frame_system, :A, 2, 1, t -> angle_to_dcm(t, :Z))
    add_axes_rotating!(frame_system, :B, 3, 2, t -> angle_to_dcm(0.5t, :X))
    add_axes_rotating!(frame_system, :C, 4, 3, t -> angle_to_dcm(0.3t, :Y))
    add_point!(frame_system, :O, 1, 1)
    add_point_dynamical!(
        frame_system, :P1, 2, 1, 1, t -> SVector(cos(t), sin(t), one(t)))
    add_point_dynamical!(
        frame_system, :P2, 3, 2, 2, t -> SVector(t^2, 2t, 3one(t)))
    add_direction!(
        frame_system, :D, 1, t -> SVector(cos(t), sin(t), zero(t)))
    return frame_system
end

function run_worker(phase)
    activate_benchmark()

    if phase == "package_load"
        measure("package load") do
            Base.eval(Main, :(using FrameTransformations))
        end
        return
    end

    Base.eval(Main, :(using FrameTransformations))
    Base.eval(Main, :(using ForwardDiff, ReferenceFrameRotations, StaticArrays))
    return Base.invokelatest(run_loaded_phase, phase)
end

function run_loaded_phase(phase)
    if phase == "frame_construction"
        measure(build_frames, "frame-system construction")
        return
    end

    frames = build_frames()
    if phase == "direct_rotation"
        measure("first direct rotation") do
            rotation6(frames, :I, :C, 0.2)
        end
    elseif phase == "prepare_rotation"
        measure("prepare rotation") do
            prepare_rotation(frames, :I, :C, Val(2))
        end
    elseif phase == "prepared_rotation_eval"
        operation = prepare_rotation(frames, :I, :C, Val(2))
        measure("first prepared rotation eval") do
            operation(0.2)
        end
    elseif phase == "compile_rotation"
        measure("compile rotation") do
            compile_rotation(frames, :I, :C, Val(2))
        end
    elseif phase == "compiled_rotation_eval"
        operation = compile_rotation(frames, :I, :C, Val(2))
        measure("first compiled rotation eval") do
            operation(0.2)
        end
    elseif phase == "direct_translation"
        measure("first direct translation") do
            vector6(frames, :O, :P2, :A, 0.2)
        end
    elseif phase == "prepare_translation"
        measure("prepare translation") do
            prepare_translation(frames, :O, :P2, :A, Val(2))
        end
    elseif phase == "prepared_translation_eval"
        operation = prepare_translation(frames, :O, :P2, :A, Val(2))
        measure("first prepared translation eval") do
            operation(0.2)
        end
    elseif phase == "compile_translation"
        measure("compile translation") do
            compile_translation(frames, :O, :P2, :A, Val(2))
        end
    elseif phase == "compiled_translation_eval"
        operation = compile_translation(frames, :O, :P2, :A, Val(2))
        measure("first compiled translation eval") do
            operation(0.2)
        end
    elseif phase == "composite"
        rotation = prepare_rotation(frames, :I, :C, Val(2))
        translation = prepare_translation(frames, :O, :P2, :A, Val(2))
        direction = prepare_direction(frames, :D, :C, Val(2))
        composite = measure("composite model construction") do
            (; rotation, translation, direction)
        end
        measure("first composite dynamics call") do
            time = 0.2
            rotated = composite.rotation(time) * Translation{2}(
                SVector(1.0, 2.0, 3.0, 0.1, 0.2, 0.3))
            SVector(rotated) + composite.translation(time) + composite.direction(time)
        end
        measure("first composite ForwardDiff") do
            ForwardDiff.derivative(
                time -> sum(composite.translation(time)[SVector(1, 2, 3)]), 0.2)
        end
    else
        throw(ArgumentError("unknown first-use phase $phase"))
    end
end

function main()
    if length(ARGS) == 2 && ARGS[1] == "--worker"
        run_worker(ARGS[2])
        return
    end

    println("Fresh-process first-use latency")
    println(rpad("phase", 36), lpad("time", 12), "     bytes")
    println(repeat("-", 68))
    for phase in (
        "package_load",
        "frame_construction",
        "direct_rotation",
        "prepare_rotation",
        "prepared_rotation_eval",
        "compile_rotation",
        "compiled_rotation_eval",
        "direct_translation",
        "prepare_translation",
        "prepared_translation_eval",
        "compile_translation",
        "compiled_translation_eval",
        "composite",
    )
        command = `$(Base.julia_cmd()) --startup-file=no $PROGRAM_FILE --worker $phase`
        run(command)
    end
end

main()
