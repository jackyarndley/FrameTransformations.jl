using DelimitedFiles, 
      RemoteFiles, 
      Interpolations

"""
    get_iers_eop(; force_download = false)

Download and parse the IERS EOP C04 data. 

The files are downloaded using the `RemoteFile` package with weekly updates. Hence, if one 
desires to force a download before the scheduled time, then set the keyword `force_download`  
to `true`.

!!! note
    The files will be downloaded from the default URL. If the user want to use
    another one, then use the specialized functions
    [`get_iers_eop_iau_1980`](@ref), [`get_iers_eop_iau_2000A`](@ref)

    See also: [`get_iers_eop_IAU2000A`](@ref)

### Returns
A structure ([`EOPData_IAU1980`](@ref) or [`EOPData_IAU2000A`](@ref), depending
on `data_type`) with the interpolations of the EOP parameters. Notice that the
interpolation indexing is set to the Julian Day.

"""
function get_iers_eop(; force_download = false)
    return get_iers_eop_iau_2000A(force_download = force_download)
end

"""
    get_iers_eop_iau_2000A(url::String="https://datacenter.iers.org/data/csv/finals2000A.all.csv"; 
        force_download = false)

Get the IERS EOP C04 IAU2000A data from the URL `url`.

If `url` is omitted, then it defaults to https://datacenter.iers.org/data/csv/finals2000A.all.csv.
The file is downloaded using the `RemoteFile` package with weekly updates. Hence, if one desires 
to force a download before the scheduled time, then set the keyword `force_download` to `true`.

!!! note
    The interpolation of every field in [`EOPData_IAU2000A`](@ref) between two
    points in the grid is linear. If extrapolation is needed, then if will use
    the nearest value (flat extrapolation).

See also: [`get_iers_eop`](@ref)

### Returns
The structure `EOPData` with the interpolations of the EOP parameters. Notice that the 
interpolation indexing is set to the Julian Day.
"""
function get_iers_eop_IAU2000A(
    url::String = "https://datacenter.iers.org/data/csv/finals2000A.all.csv";
    force_download = false
)
    @RemoteFile(
        _eop_iau2000A,
        url,
        file="EOP_IAU2000A.txt",
        dir="ext/eop",
        updates=:fridays
    )

    # Download the data
    download(_eop_iau2000A; force = force_download, force_update = true)

    # Parse the data removing the header.
    eop, ~ = readdlm(path(_eop_iau2000A), ';'; header = true)

    # Create the EOP Data structure by creating the interpolations:
    # - The interpolation will be linear between two points in the grid.
    # - The extrapolation will be flat, considering the nearest point.
    
    knots::Vector{Float64} = Vector{Float64}(eop[:, 1] .+ 2400000.5)

    return EOPData(
        _create_iers_eop_interpolation(knots, eop[:, 6]),
        _create_iers_eop_interpolation(knots, eop[:, 8]),
        _create_iers_eop_interpolation(knots, eop[:, 11]),
        _create_iers_eop_interpolation(knots, eop[:, 13]),
        _create_iers_eop_interpolation(knots, eop[:, 20]),
        _create_iers_eop_interpolation(knots, eop[:, 22]),
    )

end

"""
    EOPData{T}
EOP Data for IAU 2000A.

!!! note
    Each field will be an `AbstractInterpolation` indexed by the Julian Day.

# Fields
- `x, y`: Polar motion with respect to the crust [arcsec].
- `UT1_UTC`: Irregularities of the rotation angle [s].
- `LOD`: Length of day offset [ms].
- `dX, dY`: Celestial pole offsets referred to the model IAU2000A [milliarcsec].
"""
struct EOPData{T}
    x::T
    y::T
    UT1_UTC::T
    LOD::T
    dX::T
    dY::T
end

function Base.show(io::IO, eop::EOPData{T}) where T
    # Check if IO has support for colors.
    println(io, " ")
    println(io, "  EOPData ", "│ ", "Timespan")
    println(io, " ─────────┼──────────────────────────────────────────────")
    println(io, "        x ", "│ ", _get_iers_eop_timespan(eop.x))
    println(io, "        y ", "│ ", _get_iers_eop_timespan(eop.y))
    println(io, "  UT1-UTC ", "│ ", _get_iers_eop_timespan(eop.UT1_UTC))
    println(io, "      LOD ", "│ ", _get_iers_eop_timespan(eop.LOD))
    println(io, "       dX ", "│ ", _get_iers_eop_timespan(eop.dX))
    println(io, "       dY ", "│ ", _get_iers_eop_timespan(eop.dY))

    return nothing
end

# Get timespan
function _get_iers_eop_timespan(itp::AbstractInterpolation)
    str = string(first(first(itp.itp.knots))) * " - " *
          string(last(first(itp.itp.knots)))
    return str
end

# Create the interpolation object for the `knots` and `field` from IERS.
function _create_iers_eop_interpolation(
    knots::AbstractVector,
    field::AbstractVector
)
    # Obtain the last available index of the field.
    last_id = findlast(!isempty, field)
    last_id === nothing && (last_id = length(field))

    # Convert the field to a `Vector{Float64}`.
    field_float::Vector{Float64} = Vector{Float64}(field[1:last_id])

    # Create the interpolation object.
    interp = extrapolate(interpolate(
        (knots[1:last_id],),
        field_float,
        Gridded(Linear())),
        Flat()
    )

    return interp
end