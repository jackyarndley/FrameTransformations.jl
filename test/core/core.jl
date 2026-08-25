using Test 
using SafeTestsets
 
@testset "Core" verbose=true begin
    @safetestset "Translation" begin include("translation.jl") end
    @safetestset "Rotation" begin include("rotation.jl") end
    @safetestset "Graph" begin include("graph.jl") end
    @safetestset "Transform (generic)" begin include("transform.jl") end
    @safetestset "Registration and exact order" begin include("registration.jl") end
    @safetestset "Compiled fast-path" begin include("compiled.jl") end
end;
