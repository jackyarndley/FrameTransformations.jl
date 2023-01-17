using Basic
using Test

using RemoteFiles 

@RemoteFileSet KERNELS "Generic Spice Kernels" begin
    IAU = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc" dir=joinpath(@__DIR__, "assets")
    LEAP = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls" dir=joinpath(@__DIR__, "assets")
    GMS = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc" dir=joinpath(@__DIR__, "assets")
end

download(KERNELS, verbose=true, force=true)


@testset "Basic" verbose=true begin
    @eval begin
        modules = [:Tempo, :Orient, :Frames, :Utils]
        for m in modules
            @testset "$m" verbose=true begin 
                include("$m/$m.jl")
            end         
        end
    end
end;