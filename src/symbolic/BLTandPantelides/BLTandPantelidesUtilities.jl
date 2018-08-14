"""
Module with utility functions for BLTandPantelides.

* Author: Hilding Elmqvist, Mogram AB  
* Date: July-August 2016
* License: MIT

"""
module BLTandPantelidesUtilities

#using ..BLTandPantelides
using ..ModiaLogging 

export buildExtendedSystem, addDependencies, buildFullIncidence
export invertDer, invertAssign
export createNames, printList, printAssignedEquations, printSortedEquations, printUnassigned, makeList

"""
  function buildExtendedSystem(A)

Extend a system according to Pantelides equation (15), i.e. return the incidence for function h(x, der(x)).

* `A`: A[j] = if V[k] = der(V[j]) then k else 0
* `return G`: bipartite graph

Example:   
    julia> BLTandPantelidesUtilities.buildExtendedSystem([5,6,7,8,0,0,0,0,0])
    4-element Array{Any,1}:
     [1,5]
     [2,6]
     [3,7]
     [4,8]
"""
function buildExtendedSystem(A)
    G = []
    for i in 1:length(A)
        a = A[i]
        if a > 0
            push!(G, [i, a])    # h(x, der(x))
        end
    end
    return G
end

function addDependencies(G, Vindices)
    newG = []
    for g in G
        push!(newG, [g; Vindices])
    end  
    return newG
end

"""
    buildFullIncidence(n,m)

Build a bipartite graph with full incidence, i.e. all of the n E-nodes refer to all of the m V-nodes.

* `n`: number of E-nodes
* `m`: number of V-nodes
* `return G`: bipartite graph

Example: 
    julia> BLTandPantelidesUtilities.buildFullIncidence(2,3) 
    2-element Array{Any,1}: 
     [1,2,3] 
     [1,2,3] 
"""
function buildFullIncidence(n, m)
    G = []
    for i in 1:n
        push!(G, [j for j in 1:m])    
    end
    return G
end


"""
    function invertDer(A)

Invert derivative relationships for variables and equations

* `A`: A[j] = if V[k] = der(V[j]) then k else 0 (or correspondingly for E-nodes)
* `return orgIndex`: index of original variable or equation
* `return derOrder`: derivative order

Note that invertDer can be used to invert from list of E-nodes to list of V-nodes as well.

Example:
 julia> BLTandPantelidesUtilities.invertDer([5,6,7,8,10,11,0,0,0,0,0])
([1,2,3,4,1,2,3,4,9,1,2],[0,0,0,0,1,1,1,1,0,2,2])
"""
function invertDer(A)
    orgIndex = [i for i in 1:length(A)]  # Index of original variable or equation
    derOrder = fill(0, length(A)) # Derivative order
    for i in 1:length(A)
        a = A[i]
        if a > 0
            derOrder[a] = derOrder[i] + 1
            orgIndex[a] = orgIndex[i]
        end
    end
    return orgIndex, derOrder
end


"""
    invertAssign(assign, n=length(assign))
    
Invert assignment relationships for variables and equations.

* `assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `n`: number of E-nodes
* `return invAssign`: invAssign[i] contains the V-node to which E-node i is assigned or 0 if E-node i not assigned
* `return unAssigned`: unassigned V-nodes

Note that invertAssign can be used to invert from list of E-nodes to list of V-nodes as well.

Example:
julia> inv=BLTandPantelidesUtilities.invertAssign([0,0,0,0,1,2,7,4,3,9,8])
([5,6,9,8,0,0,7,11,10,0,0],[1,2,3,4])

julia> BLTandPantelides.invertAssign(inv[1])
([0,0,0,0,1,2,7,4,3,9,8],[5,6,10,11])
"""
function invertAssign(assign, n=length(assign))
    invAssign = fill(0, n)
    unAssigned::Vector{Int} = []
    for j in 1:length(assign)
        i = assign[j]
        if i > 0
            invAssign[i] = j
        else
            push!(unAssigned, j)
        end
    end
    return invAssign, unAssigned
end

"""
    function createNames(infixes, A)

Creates names.

* `infixes`: infix strings for original variable 
* `A`: A[j] = if V[k] = der(V[j]) then k else 0 

Example:
julia> BLTandPantelidesUtilities.createNames(["x", "y", "w", "z", "", "", "", "", "T"], [5,6,7,8,10,11,0,0,0,0,0])
x, y, w, z, der(x), der(y), der(w), der(z), T, der2(x), der2(y)

"""
function createNames(infixes, A)
    names = []
    (orgIndex, derOrder) = invertDer(A)
    for j in 1:length(A)
        if derOrder[j] > 0
            d = "der"
            if derOrder[j] > 1
                d = d * string(derOrder[j])
            end
            d = "$d($(infixes[orgIndex[j]]))"
        else
            d = infixes[orgIndex[j]]
        end
        push!(names, d)
    end
    names
end


"""
    function printList(infixes, indices, A, vertical=false)

Print list of variables or equations.

* `infixes`: infix strings for original variable or equation
* `indices`: indices for the variables or equations to be printed
* `A`: A[j] = if V[k] = der(V[j]) then k else 0 (or correspondingly for E-nodes)
* `vertical`: if vertical then new line separation else comma separation

Example:
julia> BLTandPantelidesUtilities.printList(["x", "y", "w", "z", "", "", "", "", "T"], 1:11, [5,6,7,8,10,11,0,0,0,0,0])
x, y, w, z, der(x), der(y), der(w), der(z), T, der2(x), der2(y)

"""
function printList(infixes, indices, A, vertical=false)
    (orgIndex, derOrder) = invertDer(A)
    for ind in 1:length(indices)
        j = indices[ind]
        if j > 0    
            if ind > 1
                if vertical
                    loglnModia()
                else
                    logModia(", ")
                end
            end
          
            if derOrder[j] > 0
                if vertical
                    logModia("DER")
                else
                    logModia("der")
                end
          
                if derOrder[j] > 1
                    logModia(derOrder[j])
                end
                logModia("(")
            end
          
            logModia(infixes[orgIndex[j]])
          
            if derOrder[j] > 0
                logModia(")")
            end
        end
    end
    loglnModia()
end

function makeList(infixes, indices, A, vertical=false)
    l = String[]
    (orgIndex, derOrder) = invertDer(A)
    for ind in 1:length(indices)
        s = ""
        j = indices[ind]
        if j > 0    
            if derOrder[j] > 0
                if vertical
                    s = "DER"
                else
                    s = "der"
                end
          
                if derOrder[j] > 1
                    s *= string(derOrder[j])
                end
                s *= "_" # "("
                # s *= "this."
            end

            s *= string(infixes[orgIndex[j]])
            if derOrder[j] > 0
                # s *= ")"
            end
            
            push!(l, s)
        end
    end
    return l
end

"""
    printAssignedEquations(equations, variables, indices, assign, A, B)

Print assigned equations.

* `equations`: infix string for original equations
* `variables`: infix string for original variables
* `indices`: indices for the equations to be printed
* `assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `A`: A[j] = if V[k] = der(V[j]) then k else 0
* `B`: B[i] = if E[l] = der(E[l]) then l else 0

Example:
See testBLTandPantelides.testPendulum

"""
function printAssignedEquations(equations, variables, indices, assign, A, B)
    (orgIndexVar, derOrderVar) = invertDer(A)
    (orgIndexEqu, derOrderEqu) = invertDer(B)
    (assignedVar, unAssigned) = invertAssign(assign)
    if unAssigned != []
        # @show assignedVar unAssigned
    end
    
    for i in indices
        logModia(lpad(string(i) * ":", 5, " "))
        j = assignedVar[i]
        if j > 0    
            if derOrderVar[j] == 1
                prefix = "der("
                suffix = ")"
            elseif derOrderVar[j] > 1
                prefix = "der" * string(derOrderVar[j]) * "("
                suffix = ")"
            else
                prefix = ""
                suffix = ""
            end
            logModia(lpad(prefix * string(variables[orgIndexVar[j]]) * suffix, 25, " "))
        else
            logModia(" "^25)
        end
        logModia(":  ")
      
        if derOrderEqu[i] == 1
            prefix = "DER( "
            suffix = " )"
        elseif derOrderEqu[i] > 1
            prefix = "DER" * string(derOrderEqu[i]) * "( "
            suffix = " )"
        else  
            prefix = ""
            suffix = ""
        end
        loglnModia(prefix * string(equations[orgIndexEqu[i]]) * suffix)
    end
end


"""
    printSortedEquations(equations, variables, components, assign, A, B)

Print sorted equations.

* `equations`: infix string for original equations
* `variables`: infix string for original variables
* `components`: cell array of components. Each component is a list of indices to E-nodes
* `assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `A`: A[j] = if V[k] = der(V[j]) then k else 0
* `B`: B[i] = if E[l] = der(E[l]) then l else 0

Example:
See testBLTandPantelides.testPendulum

"""
function printSortedEquations(equations, variables, components, assign, A, B)
    loglnModia("[assigned variable]: [differentiation] equation")
    loglnModia("Strongly connected components are enclosed in []")
    for c in components
        if length(c) > 1
            loglnModia("[")
        end
    
        printAssignedEquations(equations, variables, c, assign, A, B)
    
        if length(c) > 1
            loglnModia("]")
        end
    end
end


"""
    printUnassigned(equations, variables, components, assign, A, B)

Print unassigned variables and equations.

* `equations`: infix string for original equations
* `variables`: infix string for original variables
* `assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `A`: A[j] = if V[k] = der(V[j]) then k else 0
* `B`: B[i] = if E[l] = der(E[l]) then l else 0

Example:
See testBLTandPantelides.testPendulum

"""
function printUnassigned(equations, variables, assign, A, B, vActive=[])
    (invAssign, unAssignedVariables) = invertAssign(assign, length(B))
    (ass, unAssignedEquations) = invertAssign(invAssign, length(assign))
    if false # BLTandPantelides.log
        @show vActive
        @show invAssign unAssignedVariables
        @show ass unAssignedEquations
    end
    if vActive != []
    # Don't print not active variables
        unass = []
        for v in unAssignedVariables
            if vActive[v]
                push!(unass, v)
            end
        end
        unAssignedVariables = unass
    end   
    loglnModia("\nUnassigned variables:")
    printList(variables, unAssignedVariables, A)

    loglnModia("\nUnassigned equations:")
    printList(equations, unAssignedEquations, B, true)
end

end