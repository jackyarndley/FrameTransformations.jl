using Basic
using Test

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