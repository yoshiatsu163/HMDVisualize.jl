Base.@kwdef mutable struct BondData
    origin::Vector{Point3f0} = Point3f[]
    direction::Vector{Point3f0} = Point3f[]
    radius::Vector{Float32} = Float32[]
end

#function add_bdata!(b::BondData, origin::AbstractVector{F}, direction::AbstractVector{F}, radius::F) where {F<:AbstractFloat}

#end

#function add_bdata!(b::BondData, origin::Point3f0, direction::Point3f0, radius::Float32)
function add_bdata!(b::BondData, origin::AbstractVector{F}, direction::AbstractVector{F}, radius::F) where {F<:AbstractFloat}
    push!(b.origin, origin)
    push!(b.direction, direction)
    push!(b.radius, radius)
end

function bondscatter!(axis::LScene, bonds::BondData; color::NTuple{3, UInt8}, bond_radius::Number, quality::Integer)
    # bond mesh template
    # origin, extremity, radius1, radius2, segments
    m = Makie._mantle(Point3f(zeros(3)), Point3f((0,0,1)), bond_radius, bond_radius, quality)

    # scaling factor for `m`
    scales = Vec3f[]
    for v in bonds.direction
        l = norm(v)
        push!(scales, Vec3f(1/l, 1/l, l))
    end

    # rotation for `m`
    rots = normalize(bonds.direction)

    meshscatter!(axis, bonds.origin; rotation = rots, markersize = scales, color = hex(color), marker = m)
end

function bondscatter!(axis::LScene, bonds::Observable{BondData}; color::NTuple{3, UInt8}, bond_radius::Number, quality::Integer)
    # bond mesh template
    # origin, extremity, radius1, radius2, segments
    m = Makie._mantle(Point3f(zeros(3)), Point3f((0,0,1)), bond_radius, bond_radius, quality)
    println(1)

    # scaling factor for `m`
    scales = lift(bonds) do stub
        map(bonds[].direction) do v
            l = norm(v)
            Vec3f(1, 1, l)
        end
    end
    println(2)

    # rotation for `m`
    rots = lift(bonds) do stub
        normalize(bonds[].direction)
    end
    println(3)

    points = lift(bonds) do stub
        bonds[].origin
    end
    println(4)

    meshscatter!(axis, points; rotation = rots, markersize = scales, color = hex(color), marker = m)
    println(5)
end

function bondscatter2!(axis::LScene, bonds::Observable{BondData}; color::NTuple{3, UInt8}, bond_radius::Number, quality::Integer)

end

function bond_nonpbc(s::AbstractSystem, color_func::Function, colors::Vector{NTuple{3, UInt8}}, bond_radius::Number, quality::Integer)
    if wrapped(s)
        error("bond_nonpbc does not support wrapped system. Use bond_pbc instead.")
    end

    bonds = Dict(cl => BondData() for cl in colors)
    for edge in edges(topology(s))
        atom_id1, atom_id2 = src(edge), dst(edge)
        p1, p2 = position(s, atom_id1), position(s, atom_id2)
        center = p1 + 0.5*(p2 - p1)
        add_bdata!(
            bonds[color_func(s, atom_id1)],
            p1, center-p1, bond_radius
        )
        add_bdata!(
            bonds[color_func(s, atom_id2)],
            center, p2-center, bond_radius
        )
    end

    return bonds
end

function bond_pbc(s::AbstractSystem, color_func::Function, colors::Vector{NTuple{3, UInt8}}, bond_radius::Number, quality::Integer)
    if !wrapped(s)
        error("bond_pbc does not support non-wrapped system. Use bond_nonpbc instead.")
    end

    bonds = Dict(cl => BondData() for cl in colors)
    origin = box(s).origin # box origin
    a, b, c = box(s).axis[:, 1], box(s).axis[:, 2], box(s).axis[:, 3] # box vector
    for edge in edges(topology(s))
        atom_id1, atom_id2 = src(edge), dst(edge)
        p1, p2 = position(s, atom_id1), position(s, atom_id2)
        diff = travel(s, atom_id2) - travel(s, atom_id1)
        bond_vector = (diff[1]*a + diff[2]*b + diff[3]*c + p2) - p1
        @assert all(x -> x ∈ (-1, 0, 1), diff)

        # solve p1 + t(p2 - p1) == αa + βb + γc + origin
        t = Float64[]
        if diff[1] != 0
            α = diff[1]
            (_t, β, γ) = (-p1 + α*a + origin) \ [bond_vector -b -c]
            push!(t, _t)
        elseif diff[2] != 0
            β = diff[2]
            (_t, α, γ) = (-p1 + β*b + origin) \ [bond_vector -a -c]
            push!(t, _t)
        elseif diff[3] != 0
            γ = diff[3]
            (_t, α, β) = (-p1 + γ*c + origin) \ [bond_vector -a -b]
            push!(t, _t)
        else
            t = [0.0, 1.0]
        end
        @assert all(x -> -1 <= x <= 1, t)

        # add center point of bond
        sort!(t)
        ic = searchsortedfirst(t, 0.5)
        insert!(t, ic, 0.5)

        # add bond data
        color1, color2 = color_func(s, atom_id1), color_func(s, atom_id2)
        for i in 1:length(t)-1
            p  = p1 + t[i]*bond_vector
            _p = p1 + t[i+1]*bond_vector
            cl = t[i] < 0.5 ? color1 : color2
            add_bdata!(
                bonds[cl],
                p, _p-p, bond_radius
            )
        end
    end

    return bonds
end
