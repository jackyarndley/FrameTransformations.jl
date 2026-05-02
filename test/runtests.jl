const GROUP = replace(lowercase(get(ENV, "GROUP", "all")), "/" => "", "\\" => "")
const GROUPS = ("all", "core", "definitions")

GROUP in GROUPS || throw(ArgumentError("GROUP must be one of $(join(GROUPS, ", ")); got $(repr(GROUP))"))

@time begin
    if GROUP == "all" || GROUP == "core"
        include("core/core.jl")
    end
    if GROUP == "all" || GROUP == "definitions"
        include("definitions/definitions.jl")
    end
end;
