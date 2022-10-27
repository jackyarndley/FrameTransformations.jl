module Orient

    using StaticArrays
    using LinearAlgebra
    using ReferenceFrameRotations

    using Basic.Tempo 
    using Basic.Utils: arcsec2rad

    include("abstract.jl")
    include("iau.jl")
    include("types.jl")

    include("Earth/long_term.jl")
    include("Earth/obliquity.jl")
    include("Earth/precession.jl")

    include("ecliptic.jl")
end