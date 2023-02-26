module HMDVisualize

using Colors
# Choose colorschemes with care! Refer to Peter Kovesi's PerceptualColourMaps package, or to Fabio Crameri's Scientific Colour Maps for more information.
using ColorSchemes: colorschemes, get
using GeometryBasics
using Graphs
using HMD
using LinearAlgebra
#using MakieCore
using GLMakie #cite
using PeriodicTable

import GeometryBasics: Cylinder

export visualize, color_scheme

const atom_color = Dict(
    elements[:H ].symbol => "#" * hex(HSL(170/255, 0  /255, 240/255)),
    elements[:C ].symbol => "#" * hex(HSL(170/255, 0  /255, 100/255)),
    elements[:N ].symbol => "#" * hex(HSL(170/255, 75 /255, 150/255)),
    elements[:O ].symbol => "#" * hex(HSL(0  /255, 75 /255, 150/255)),
    elements[:F ].symbol => "#" * hex(HSL(80 /255, 75 /255, 150/255)),
    elements[:Si].symbol => "#" * hex(HSL(10 /255, 100/255, 190/255)),
    elements[:S ].symbol => "#" * hex(HSL(38 /255, 208/255, 145/255)),
    elements[:Cl].symbol => "#" * hex(HSL(113/255, 124/255, 145/255)),
    "beads"              => "#" * hex(HSL(170/255, 0  /255, 100/255))
)

# 色情報にアルファを追加 8桁hex
# trajectoryの動画出力
# 投影法
# 背景色にあわせた色調整
# マウスイベントによるinspection (特に原子間距離)
function visualize(s::AbstractSystem; wrap_coord::Bool=false, add_color::Dict{<:Integer, String}=Dict{Int64, String}())
    if dimension(s) != 3
        error("expected dimension 3, found $D")
    end

    changed = false
    if wrap_coord && !wrapped(s)
        wrap!(s)
        changed = true
    elseif !wrap_coord && wrapped(s)
        unwrap!(s)
        changed = true
    end

    f = Figure()
    axis = LScene(f[1,1]; show_axis = false
            #xticksvisible = false,
            #yticksvisible = false,
            #zticksvisible = false,
            #xticklabelsvisible = false,
            #yticklabelsvisible = false,
            #zticklabelsvisible = false,
            #xgridvisible = false,
            #ygridvisible = false,
            #zgridvisible = false,
            #aspect = :data,
            #viewmode = :fit
    )

    render_box!(axis, s)

    #mmeshを1つに結合する方法を調べる geometrybasics
    #for (pos, e) in zip(all_positions(s), all_elements(s))
    #    meshscatter!(axis, pos
    #                color = atom_color[string(e)],
    #                markersize = 0.3
    #    )
    #end
    set_color(atom_id) = begin
        if atom_id ∈ keys(add_color)
            add_color[atom_id]
        else
            atom_color[string(element(s, atom_id))]
        end
    end
    meshscatter!(axis, Point3f.(all_positions(s)), color=set_color.(1:natom(s)), markersize = 0.3)

    # mesh単位でメタデータ新納情報を持たせ、プロット直前にatomともmerge
    bmesh = [normal_mesh(Cylinder(zeros(3), ones(3), 0.15))] |> empty
    colors = String[]
    for (edge, n) in zip(edges(topology(s)), 1:nbond(s))
        atom_id1, atom_id2 = src(edge), dst(edge)
        #render_bond!(axis, s, atom_id1, atom_id2; radius=0.15, wrap_coord=wrap_coord)
        bm, each_color = bond_mesh(s, atom_id1, atom_id2; radius=0.15, wrap_coord=wrap_coord, color_func=set_color)
        append!(bmesh, bm)
        append!(colors, each_color)
    end
    for cl in unique(colors)
        color_mesh = [bmesh[i] for i in 1:length(colors) if colors[i] == cl]
        mesh!(axis, merge(color_mesh); color=cl)
    end
    
    if changed && wrapped(s)
        unwrap!(s)
    elseif changed && !wrapped(s)
        wrap!(s)
    end

    return f
end

function render_bond!(axis, s::AbstractSystem, atom_id1::Integer, atom_id2::Integer; radius::AbstractFloat, wrap_coord::Bool)
    # boxとbondの交点
    # [(α, p, bool)]
    points = if wrap_coord
        separate_points(s, atom_id1, atom_id2)
    else
        [(len=0.0, point=zeros(3), jump=false)] |> empty!
    end

    # 始点と終点を追加
    push!(points, (len=0.0, point=position(s, atom_id1), jump=false))
    push!(points, (len=1.0, point=position(s, atom_id2), jump=false))

    #println("$(norm(points[1].point .- points[end].point))")
    #for p in points
    #    println("   $(p.len) $(p.point)")
    #end
    #println()

    # 色分けのため中点を追加
    push!(points, (len=0.5, point=zeros(3), jump=false))
    sort!(points, by=tuple->tuple[1]) # sort by α
    i = findfirst(tuple -> tuple.len==0.5, points)
    α_center = (0.5 .- points[i-1].len) / (points[i+1].len .- points[i-1].len)
    points[i] = (len = points[i].len,
                point = points[i-1].point .+ α_center .* (points[i+1].point .- points[i-1].point),
                jump = false)

    #
    p_ = points[2]
    bmesh = [Cylinder(points[1].point, p_.point, radius)]
    for p in points[3:end]
        if !p.jump
            push!(bmesh, Cylinder(p_.point, p.point, radius))
        end
        p_ = p
    end
    bo = bond_order(s, atom_id1, atom_id2)

    # 描画
    for i in 1:length(bmesh)
        α, bm = points[i].len, bmesh[i]
        bond_color = if α < 0.5
            atom_color[string(element(s, atom_id1))]
        else
            atom_color[string(element(s, atom_id2))]
        end
        if bo == 3//2
            mesh!(axis, dashed_bmesh(bm); color = bond_color)
        elseif bo == 2//1
            for b in split_bmesh(bmesh, 2)
                mesh!(axis, b; color = bond_color)
            end
        elseif bo == 3//1
            for b in split_bmesh(bmesh, 3)
                mesh!(axis, b; color = bond_color)
            end
        else
            mesh!(axis, bm; color = bond_color)
        end
    end
end

function bond_mesh(s::AbstractSystem, atom_id1::Integer, atom_id2::Integer; radius::AbstractFloat, wrap_coord::Bool, color_func::Function)
    # boxとbondの交点
    # [(α, p, bool)]
    points = if wrap_coord
        separate_points(s, atom_id1, atom_id2)
    else
        [(len=0.0, point=zeros(3), jump=false)] |> empty!
    end

    # 始点と終点を追加
    push!(points, (len=0.0, point=position(s, atom_id1), jump=false))
    push!(points, (len=1.0, point=position(s, atom_id2), jump=false))

    #println("$(norm(points[1].point .- points[end].point))")
    #for p in points
    #    println("   $(p.len) $(p.point)")
    #end
    #println()

    # 色分けのため中点を追加
    push!(points, (len=0.5, point=zeros(3), jump=false))
    sort!(points, by=tuple->tuple[1]) # sort by α
    i = findfirst(tuple -> tuple.len==0.5, points)
    α_center = (0.5 .- points[i-1].len) / (points[i+1].len .- points[i-1].len)
    points[i] = (len = points[i].len,
                point = points[i-1].point .+ α_center .* (points[i+1].point .- points[i-1].point),
                jump = false)

    #
    p_ = points[2]
    bmesh = [Cylinder(points[1].point, p_.point, radius)]
    for p in points[3:end]
        if !p.jump
            push!(bmesh, Cylinder(p_.point, p.point, radius))
        end
        p_ = p
    end
    bo = bond_order(s, atom_id1, atom_id2)

    # mesh vector
    meshes = [bmesh[1]] |> empty
    colors = String[]
    for i in 1:length(bmesh)
        α, bm = points[i].len, bmesh[i]
        bond_color = α < 0.5 ? color_func(atom_id1) : color_func(atom_id2)
        if bo == 3//2
            push!(meshes, dashed_bmesh(bm))
            push!(colors, bond_color)
        elseif bo == 2//1
            append!(meshes, split_bmesh(bm, 2))
            append!(colors, fill(bond_color, 2))
        elseif bo == 3//1
            append!(meshes, split_bmesh(bm, 3))
            append!(colors, fill(bond_color, 3))
        else
            push!(meshes, bm)
            push!(colors, bond_color)
        end
    end

    return normal_mesh.(meshes), colors
end

function separate_points(s::AbstractSystem, atom_id1::Integer, atom_id2::Integer)
    axis = box(s).axis
    a, b, c = box(s).axis[:,1], box(s).axis[:,2], box(s).axis[:,3]
    b_origin = box(s).origin
    n = travel(s, atom_id2) .- travel(s, atom_id1)

    # pos: box vectorを基底とした時の成分表示() 片方はunwrapped
    pos1 = position(s, atom_id1)
    pos2_unwrap = position(s, atom_id2) .+ (n[1] .* a) .+ (n[2] .* b) .+ (n[3] .* c)
    bond_vector = pos2_unwrap .- pos1

    # pos1 + α(pos2_unwrap - pos1)の第dim成分が1 or 0のときboxとbondが交差
    # 交点は交差一回ごとに2つできるのでslideで点を複製する
    params = Vector{Tuple{Int16, Float64, Float64, Float64}}(undef, 0)
    for dim in 1:3
        if n[dim] == 0
            continue # axis[:,dim]方向には境界をまたいでいない
        end
        basis = [axis[:,i] for i in 1:3 if i != dim]
        plane_origin = b_origin .+ (n[dim] == 1 ? 1.0 : 0.0) .* axis[:,dim]
        (s, t, α) = cross_plane_line(pos1, bond_vector; basis=basis, plane_origin=plane_origin)
        push!(params, (dim, s, t, α))
        @assert 0.0 < α < 1.0
    end
    sort!(params, by=tuple->tuple[4]) # sort by alpha

    slide = zeros(3)
    cross_points = [(len=0.0, point=zeros(3), jump=false)] |> empty!
    for (dim, s, t, α) in params
        p = pos1 .+ (α .* bond_vector) .+ slide
        slide .-=  n[dim] .* axis[:,dim]
        p_reflected = pos1 .+ (α .* bond_vector) .+ slide
        push!(cross_points, (len=α, point=p, jump=false))
        push!(cross_points, (len=α, point=p_reflected, jump=true))
    end

    return cross_points
end

function cross_plane_line(origin::AbstractVector{<:AbstractFloat}, direction::AbstractVector{<:AbstractFloat}; basis, plane_origin)
    a, b = basis[1], basis[2]
    # plane_origin + s*a + t*b == origin + α*direction
    # [a b -direction] * [s; t; α] == origin - plane_origin
    (s, t, α) = [a b -direction] \ (origin .- plane_origin)

    return (s, t, α)
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

function box_coords(s::AbstractSystem, x::AbstractVector{<:AbstractFloat})
    axis = box(s).axis
    origin = box(s).origin

    # p = c[1] .* axis[:,1] .+ c[2] .* axis[:,2] .+ ...
    p = x .- origin
    e_i_e_j = [dot(axis[:,i], axis[:,j]) for i in 1:3, j in 1:3] |> Symmetric
    c = e_i_e_j \ [dot(p, axis[:,dim]) for dim in 1:3]

    return c
end

function real_coords(s::AbstractSystem, c::AbstractVector{<:AbstractFloat})
    axis = box(s).axis
    return box(s).origin .+ mapreduce(dim -> c[dim] .* axis[:, dim], .+, 1:3)
end

function dashed_bmesh(bmesh::Cylinder3)
    bmesh
end

function render_box!(axis, s)
    a, b, c = box(s).axis[:,1], box(s).axis[:,2], box(s).axis[:,3]
    s_origin = box(s).origin
    p(i, j, k) = s_origin .+ (i .* a) .+ (j .* b) .+ (k .* c)

    line_bewteen!(axis, p(0,0,0), p(0,1,0))
    line_bewteen!(axis, p(0,0,0), p(1,0,0))
    line_bewteen!(axis, p(0,1,0), p(1,1,0))
    line_bewteen!(axis, p(1,0,0), p(1,1,0))
    line_bewteen!(axis, p(0,0,1), p(0,1,1))
    line_bewteen!(axis, p(0,0,1), p(1,0,1))
    line_bewteen!(axis, p(0,1,1), p(1,1,1))
    line_bewteen!(axis, p(1,0,1), p(1,1,1))
    line_bewteen!(axis, p(0,0,0), p(0,0,1))
    line_bewteen!(axis, p(1,0,0), p(1,0,1))
    line_bewteen!(axis, p(0,1,0), p(0,1,1))
    line_bewteen!(axis, p(1,1,0), p(1,1,1))

    return nothing
end

function line_bewteen!(axis, p1, p2)
    lines!(axis, [[p1[i], p2[i]] for i in 1:3]...)
    return nothing
end

function Cylinder(p1::Vector{<:AbstractFloat}, p2::Vector{<:AbstractFloat}, radius::AbstractFloat)
    return Cylinder(Point{3, Float32}(p1), Point{3, Float32}(p2), Float32(radius))
end

function color_scheme(value::AbstractFloat; scheme=:viridis)
    return "#" * hex(get(colorschemes[scheme], value))
end


end # module
