module HMDVisualize

using Colors
# Choose colorschemes with care! Refer to Peter Kovesi's PerceptualColourMaps package, or to Fabio Crameri's Scientific Colour Maps for more information.
using ColorSchemes: colorschemes, get
using GeometryBasics
using Graphs
using HMD
using HMDPolymer
using LinearAlgebra
#using MakieCore
using GLMakie #cite
using PeriodicTable

import GeometryBasics: Cylinder
import Colors: hex

export visualize, color_scheme

###
###### preset colors
###

const atom_color = Dict(
    elements[:H ].number => (UInt8(170), UInt8(  0), UInt8(240)),
    elements[:C ].number => (UInt8(170), UInt8(  0), UInt8(100)),
    elements[:N ].number => (UInt8(170), UInt8( 75), UInt8(150)),
    elements[:O ].number => (UInt8(  0), UInt8( 75), UInt8(150)),
    elements[:F ].number => (UInt8( 80), UInt8( 75), UInt8(150)),
    elements[:Si].number => (UInt8( 10), UInt8(100), UInt8(190)),
    elements[:S ].number => (UInt8( 38), UInt8(208), UInt8(145)),
    elements[:Cl].number => (UInt8(113), UInt8(124), UInt8(145))
)

function default_color(s::AbstractSystem, atom_id::Integer)
    elem = element(s, atom_id)
    if elem ∉ keys(atom_color)
        error("element $(elements[elem]) is not supported yet. ")
    end
    return if elem >= 1
        atom_color[elem]
    else
        (UInt8(170), UInt8(0), UInt8(100))
    end
end

function hex(num::NTuple{3, UInt8})
    return "#" * hex(HSL(num[1]/255, num[2]/255, num[3]/255))
end

include("bond.jl")

###
###### system visualization functions
###

function visualize(traj::AbstractTrajectory{D, F}; color_func::Function=default_color, atom_radius::Number=0.3, bond_radius::Number=0.275, quality::Integer=8) where {D, F<:AbstractFloat}
    if dimension(traj[1]) != 3
        error("expected dimension 3, found $D")
    end
    wrap_coord = wrapped(traj)

    # enumerate all possible colors
    colors = NTuple{3, UInt8}[]
    for reader in traj
        snapshot = reader.reader
        append!(colors, [color_func(snapshot, atom_id) for atom_id in 1:natom(snapshot) if color_func(snapshot, atom_id) ∉ colors])
    end

    fig = Figure()
    sl_x = Slider(fig[2, 1], range = 1:length(traj), startvalue = 1)
    axis = LScene(fig[1,1]; show_axis = false)
    cam3d!(axis; projectiontype = :orthographic, mouse_translationspeed=0.001f0)

    reader = similar_system(traj)

    # box and atoms
    box_mesh = lift(sl_x.value) do index
        update_reader!(reader, traj, index)
        get_boxmesh(reader)
    end
    atoms = lift(box_mesh) do stub
        Point3f.(all_positions(reader))
    end
    atom_colors = lift(box_mesh) do stub
        [hex(color_func(reader, atom_id)) for atom_id in 1:natom(reader)]
    end

    # bonds
    bonds = lift(box_mesh) do stub
        if wrapped(traj)
            bond_pbc(reader, color_func, colors, bond_radius, quality)
        else
            bond_nonpbc(reader, color_func, colors, bond_radius, quality)
        end
    end
    
    meshscatter!(axis, atoms;
        color = atom_colors,
        markersize = atom_radius
    )
    for cl in colors
        c_bonds = @lift $(bonds)[cl]
        bondscatter!(axis, c_bonds; color=cl, bond_radius=bond_radius, quality=quality)
    end
    mesh!(axis, box_mesh)

    return fig
end

function split_bmesh(bmesh::Cylinder3, order::Integer)
    @assert order ∈ [2,3]
    o = origin(bmesh)
    e = extremity(bmesh)
    radius = radius(bmesh)

    if order == 2
        offset = (radius / 3) .* [1.0, 0.0, 0.0]
        return [Cylinder(o .+ offset, e .+ offset, radius / 6),
                Cylinder(o .- offset, e .- offset, radius / 6)
        ]
    elseif order == 3
        offset = (0.4 * radius) .* [1.0, 0.0, 0.0]
        return [Cylinder(o .+ offset, e .+ offset, radius / 5),
                Cylinder(o          , e          , radius / 5),
                Cylinder(o .+ offset, e .+ offset, radius / 5)
        ]
    end
end

function get_boxmesh(s)
    a, b, c = box(s).axis[:,1], box(s).axis[:,2], box(s).axis[:,3]
    s_origin = box(s).origin
    p(i, j, k) = Point{3, Float32}(s_origin .+ (i .* a) .+ (j .* b) .+ (k .* c))

    #line_bewteen!(axis, p(0,0,0), p(0,1,0))
    #line_bewteen!(axis, p(0,0,0), p(1,0,0))
    #line_bewteen!(axis, p(0,1,0), p(1,1,0))
    #line_bewteen!(axis, p(1,0,0), p(1,1,0))
    #line_bewteen!(axis, p(0,0,1), p(0,1,1))
    #line_bewteen!(axis, p(0,0,1), p(1,0,1))
    #line_bewteen!(axis, p(0,1,1), p(1,1,1))
    #line_bewteen!(axis, p(1,0,1), p(1,1,1))
    #line_bewteen!(axis, p(0,0,0), p(0,0,1))
    #line_bewteen!(axis, p(1,0,0), p(1,0,1))
    #line_bewteen!(axis, p(0,1,0), p(0,1,1))
    #line_bewteen!(axis, p(1,1,0), p(1,1,1))

    return merge(normal_mesh.([
        Cylinder(p(0,0,0), p(0,1,0), Float32(0.06)),
        Cylinder(p(0,0,0), p(1,0,0), Float32(0.06)),
        Cylinder(p(0,1,0), p(1,1,0), Float32(0.06)),
        Cylinder(p(1,0,0), p(1,1,0), Float32(0.06)),
        Cylinder(p(0,0,1), p(0,1,1), Float32(0.06)),
        Cylinder(p(0,0,1), p(1,0,1), Float32(0.06)),
        Cylinder(p(0,1,1), p(1,1,1), Float32(0.06)),
        Cylinder(p(1,0,1), p(1,1,1), Float32(0.06)),
        Cylinder(p(0,0,0), p(0,0,1), Float32(0.06)),
        Cylinder(p(1,0,0), p(1,0,1), Float32(0.06)),
        Cylinder(p(0,1,0), p(0,1,1), Float32(0.06)),
        Cylinder(p(1,1,0), p(1,1,1), Float32(0.06))
    ]))
end

function line_bewteen!(axis, p1, p2)
    lines!(axis, [[p1[i], p2[i]] for i in 1:3]...)
    return nothing
end

function color_scheme(value::AbstractFloat; scheme=:viridis)
    return "#" * hex(get(colorschemes[scheme], value))
end

function update_reader!(reader, traj, index)
    import_dynamic!(reader, traj, index)
    if is_reaction(traj, index)
        import_static!(reader, traj, index)
    end

    return nothing
end


end # module
