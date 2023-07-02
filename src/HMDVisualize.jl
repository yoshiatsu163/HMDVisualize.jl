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

export visualize, color_scheme, default_color

###
###### preset colors
###

const atom_color = Dict(
    elements[:H ].number => colorant"hsla(170,   0%,  66%, 1.0)",
    elements[:C ].number => colorant"hsla(240,  39%,  30%, 1.0)",
    elements[:N ].number => colorant"hsla(240, 100%,  78%, 1.0)",
    elements[:O ].number => colorant"hsla( 80,  29%,  58%, 1.0)",
    elements[:F ].number => colorant"hsla(  0,  29%,  58%, 1.0)",
    elements[:Si].number => colorant"hsla( 10,  39%,  74%, 1.0)",
    elements[:S ].number => colorant"hsla(113,  48%,  56%, 1.0)",
    elements[:Cl].number => colorant"hsla(170,   0%,  94%, 1.0)"
)

function default_color(s::AbstractSystem, atom_id::Integer)
    elem = element(s, atom_id)
    return if elem >= 1
        atom_color[elem]
    else elem < 0
        colorant"hsla(170,  0%, 37%, 1.0)"
    end
end

include("bond.jl")

###
###### system visualization functions
###

function visualize(s::AbstractSystem{D, F, SysType}; color_func::Function=default_color, atom_radius::Number=0.3, bond_radius::Number=0.275, quality::Integer=8) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    return visualize(
        Trajectory(s);
        color_func = color_func,
        atom_radius = atom_radius,
        bond_radius = bond_radius,
        quality = quality
    )
end

# 二重結合，共鳴用のassetsを準備
# overlayで追加の色指定
# datainspector
function visualize(traj::AbstractTrajectory{D, F, SysType}; color_func::Function=default_color, overlay=nothing, atom_radius::Number=0.3, bond_radius::Number=0.275, quality::Integer=8) where {D, F<:AbstractFloat, SysType<:AbstractSystemType}
    if dimension(traj[1]) != 3
        error("expected dimension 3, found $D")
    end
    wrap_coord = wrapped(traj)

    # enumerate all possible colors
    colors = Vector{HSLA{Float32}}(undef, 0)
    for reader in traj
        snapshot = reader.reader
        for atom_id in 1:natom(snapshot)
            color = color_func(snapshot, atom_id)
            if color ∉ colors
                push!(colors, color)
            end
        end
    end

    fig = Figure(; backgroundcolor = :black)
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
        [color_func(reader, atom_id) for atom_id in 1:natom(reader)]
    end

    # bonds
    bonds = lift(box_mesh) do stub
        if wrapped(traj)
            bond_pbc(reader, color_func, colors, bond_radius, quality)
        else
            bond_nonpbc(reader, color_func, colors, bond_radius, quality)
        end
    end

    inspect_labels = map(1:natom(reader)) do i
        elem = element(reader, i)
        p = Float32.(position(reader, i))
        "id: $i\n\
        element: $(elements[elem].symbol)\n\
        x: $(p[1])\n\
        y: $(p[2])\n\
        z: $(p[3])"
    end
    meshscatter!(axis, atoms;
        color = atom_colors,
        markersize = atom_radius*2,
        inspector_label = (self, i, p) -> inspect_labels[i]
    )
    for cl in colors
        c_bonds = @lift $(bonds)[cl]
        bondscatter!(axis, c_bonds; color=cl, bond_radius=bond_radius, quality=quality)
    end
    mesh!(axis, box_mesh; color = :white)

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
