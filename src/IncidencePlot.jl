"""
Module to produce incidence plots

* Developer: Hilding Elmqvist, Mogram AB
* First version: December 2020
* License: MIT (expat)

"""
module IncidencePlot

using GLMakie
using AbstractPlotting
using GeometryBasics
using Colors

# include("...ModiaBase/src/Tearing.jl")
# using .Tearing

export showIncidence

removeBlock(ex) = ex
function removeBlock(ex::Expr)
    if ex.head in [:block]
        ex.args[1]
    else
        Expr(ex.head, [removeBlock(arg) for arg in ex.args]...)
    end
end

AbstractPlotting.inline!(false)

"""
    showIncidence(equations, unknowns, G, vActive, substitutions, assigned, assignedVar, blocks; sortVariables=false)

Shows the incidence matrix using the plotting package Makie. If `assigned` is not [], the assigned variable is marked. `blocks` is a vector of vectors
enabling showing of BLT-blocks. After each incidence matrix is output, the user is promted for a RETURN, not to immediately overwriting the plot.

"""
function showIncidence(equations, unknowns, G, vActive, substitutions, assigned, assignedVar, blocks; sortVariables=false)
    scene = Scene(camera = campixel!, show_axis = false, resolution = (1000, 800))
    margin = 20
    size = 600

    nVar = length(unknowns)
    if sortVariables
        sortedVariables = fill(0, nVar)
        varCount = 0
        for b in blocks
            for e in b
                varCount += 1
                sortedVariables[varCount] = assignedVar[e]
            end
        end
        for i in 1:length(vActive)
            if ! vActive[i]
                varCount += 1
                sortedVariables[varCount] = i
            end
        end
    else
        sortedVariables = 1:nVar
    end

    n=length(equations)
    m=length(unknowns)

    width = min(div(size, n), 20)
    height = min(div(size, m), 20)

    function drawText(s, row, column; align=:left, vertical=false)
        rot = if ! vertical; 0.0 else pi/2 end
        text!(scene, s, position = Vec(margin + column*width + width/2, margin + (n+1-row)*height + height/2),
            color = :green, textsize = 15, align = (align, :center), rotation = rot, show_axis = false)
    end

    # -----------------------------------------------------------------------------

    # Draw matrix grid
    for i in 1:n+1
        poly!([Point2f0(margin + width, margin + i*height), Point2f0(margin + (m+1)*width, margin + i*height)], strokecolor = :black)
    end

    for j in 1:m+1
        poly!([Point2f0(margin + j*width, margin + height), Point2f0(margin + j*width, margin + (n+1)*height)], strokecolor = :black)
    end

    for j in 1:length(unknowns)
        # Draw variable number under matrix
        drawText(string(sortedVariables[j]), n+1, j, align=:right, vertical=true)

        # Draw variable names above matrix
        v = unknowns[sortedVariables[j]]
        var = string(v)
        if v in keys(substitutions)
            subs = replace(replace(replace(string(substitutions[v]), "(" => ""), ")" => ""), "--" => "")
            var *= " = " * subs
        end
        drawText(var, 0, j, vertical=true)
    end

    row = 0
    for b in blocks
        for e in b
            row += 1
            # Draw equation number to the left
            drawText(string(e), row, 0, align=:right)

            # Draw equation to the right
            equ = equations[e]
            Base.remove_linenums!(equ)
            equ = removeBlock(equ)
            equ = string(equ)
            equ = replace(equ, "\n" => ";")
            if equ != :(0 = 0)
                drawText(equ, row, m+1)
            end
        end
    end

    row = 0
    for b in blocks
        if length(b) > 1
#=
            println("\nTearing of BLT block")
            @show b
            es::Array{Int64,1} = b
            vs = assignedVar[b]
            td = EquationTearing(G, length(unknowns))
            (eSolved, vSolved, eResidue, vTear) = tearEquations!(td, (e,v) -> v in G[e], es, vs)
            @show eSolved vSolved eResidue vTear
            @show equations[eSolved] unknowns[vSolved] equations[eResidue] unknowns[vTear]
            b = vcat(eSolved, eResidue)
            @show b
            blocksVars = vcat(vSolved, vTear)
=#
            color = :lightgrey
            poly!(Rect(Vec(margin + (row+1)*width, margin + (n+1-row)*height-length(b)*height),
                Vec(length(b)*width, length(b)*height)), color = color, strokecolor = :black, transperancy=true)
        end
        for e in b
            row += 1
            incidence = fill(false, nVar)
            for v in G[e]
                if sortedVariables[v] > 0
                    incidence[sortedVariables[v]] = true
                end
            end

            col = 0
            for v in sortedVariables
                col += 1
                ins = if v == 0 || ! vActive[v]; 0 elseif sortedVariables[v] > 0 && ! incidence[sortedVariables[v]]; 0 elseif v > 0 && v<= length(assigned) && assigned[v] == e; 1 else 2 end
                if ins > 0
                    color = if ins == 1; :red else :blue end
                    poly!(Rect(Vec(margin + col*width, margin + (n+1-row)*height), Vec(width, height)), color = color, strokecolor = :black)
                end
            end
        end
    end

    display(scene)

    println("Continue?")
    readline()
end

function testShowIncidence()
    # RCRL circuit
    equations = Expr[:(ground.p.v = V.n.v), :(V.p.v = Ri.n.v), :(Ri.p.v = R.p.v), :(R.n.v = C.p.v), :(C.n.v = ground.p.v), :(ground.p.i + V.n.i + C.n.i = 0), :(V.p.i + Ri.n.i = 0), :(Ri.p.i + R.p.i = 0), :(R.n.i + C.p.i = 0), :(Ri.v = Ri.p.v - Ri.n.v), :(0 = Ri.p.i + Ri.n.i), :(Ri.i = Ri.p.i), :(Ri.R * Ri.i = Ri.v), :(R.v = R.p.v - R.n.v), :(0 = R.p.i + R.n.i), :(R.i = R.p.i), :(R.R * R.i = R.v), :(C.v = C.p.v - C.n.v), :(0 = C.p.i + C.n.i), :(C.i = C.p.i), :(C.C * der(C.v) = C.i), :(V.v = V.p.v - V.n.v), :(0 = V.p.i + V.n.i), :(V.i = V.p.i), :(V.v = 10), :(ground.p.v = 0)]
    unknowns = Any[:(Ri.v), :(Ri.i), :(Ri.p.i), :(Ri.p.v), :(Ri.n.i), :(Ri.n.v), :(R.v), :(R.i), :(R.p.i), :(R.p.v), :(R.n.i), :(R.n.v), :(C.v), :(der(C.v)), :(C.i), :(C.p.i), :(C.p.v), :(C.n.i), :(C.n.v), :(V.v), :(V.i), :(V.p.i), :(V.p.v), :(V.n.i), :(V.n.v), :(ground.p.i), :(ground.p.v)]
    G = [[27, 25], [23, 6], [4, 10], [12, 17], [19, 27], [26, 24, 18], [22, 5], [3, 9], [11, 16], [1, 4, 6], [3, 5], [2, 3], [2, 1], [7, 10, 12], [9, 11], [8, 9], [8, 7], [13, 17, 19], [16, 18], [15, 16], [14, 15], [20, 23, 25], [22, 24], [21, 22], [20], [27]]
    vActive = Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    substitutions = Dict{Any,Any}()
    assigned = Any[]
    assignedVar = Any[]
    blocks = [[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12], [13], [14], [15], [16], [17], [18], [19], [20], [21], [22], [23], [24], [25], [26]]
    sortVariables = false

    showIncidence(equations, unknowns, G, vActive, substitutions, assigned, assignedVar, blocks, sortVariables=sortVariables)

    # RCRL
    equations = Expr[:(ground.p.v = V.n.v), :(V.p.v = Ri.n.v), :(Ri.p.v = R.p.v), :(R.n.v = C.p.v), :(C.n.v = ground.p.v), :(ground.p.i + V.n.i + C.n.i = 0), :(V.p.i + Ri.n.i = 0), :(Ri.p.i + R.p.i = 0), :(R.n.i + C.p.i = 0), :(Ri.v = Ri.p.v - Ri.n.v), :(0 = Ri.p.i + Ri.n.i), :(Ri.i = Ri.p.i), :(Ri.R * Ri.i = Ri.v), :(R.v = R.p.v - R.n.v), :(0 = R.p.i + R.n.i), :(R.i = R.p.i), :(R.R * R.i = R.v), :(C.v = C.p.v - C.n.v), :(0 = C.p.i + C.n.i), :(C.i = C.p.i), :(C.C * der(C.v) = C.i), :(V.v = V.p.v - V.n.v), :(0 = V.p.i + V.n.i), :(V.i = V.p.i), :(V.v = 10), :(ground.p.v = 0)]
    unknowns = Any[:(Ri.v), :(Ri.i), :(Ri.p.i), :(Ri.p.v), :(Ri.n.i), :(Ri.n.v), :(R.v), :(R.i), :(R.p.i), :(R.p.v), :(R.n.i), :(R.n.v), :(C.v), :(der(C.v)), :(C.i), :(C.p.i), :(C.p.v), :(C.n.i), :(C.n.v), :(V.v), :(V.i), :(V.p.i), :(V.p.v), :(V.n.i), :(V.n.v), :(ground.p.i), :(ground.p.v)]
    G = [[27, 25], [23, 6], [4, 10], [12, 17], [19, 27], [26, 24, 18], [22, 5], [3, 9], [11, 16], [1, 4, 6], [3, 5], [2, 3], [2, 1], [7, 10, 12], [9, 11], [8, 9], [8, 7], [13, 17, 19], [16, 18], [15, 16], [14, 15], [20, 23, 25], [22, 24], [21, 22], [20], [27]]
    vActive = Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    substitutions = Dict{Any,Any}()
    assigned = [10, 13, 12, 3, 11, 2, 17, 16, 8, 14, 15, 4, 0, 21, 20, 9, 18, 19, 5, 25, 24, 7, 22, 23, 1, 6, 26]
    assignedVar = [25, 6, 4, 12, 19, 26, 22, 9, 16, 1, 5, 3, 2, 10, 11, 8, 7, 17, 18, 15, 14, 23, 24, 21, 20, 27, 0]
    blocks = [[26], [1], [25], [22], [2], [5], [18], [4], [10, 13, 12, 8, 16, 17, 14, 3], [11], [7], [23], [15], [9], [19], [6], [20], [21], [24]]
    sortVariables = true
    showIncidence(equations, unknowns, G, vActive, substitutions, assigned, assignedVar, blocks, sortVariables=sortVariables)

end

# testShowIncidence()

end