"""
Module with graph theoretical methods for maximum matching, Pantelides index reduction and Tarjan's strongly connected components.

* Author: Hilding Elmqvist, Mogram AB  
* Date: July-August 2016 (rewritten). Array version: Spring 2017
* License: MIT

A bipartite graph (E, V) is defined by an array of arrays. Each E-entry has an integer array of indices to V-entries. 

Example bipartite graph:
 
    G = [
      [3, 5],
      [4, 6],  
      [1, 7, 9],
      [2, 8, 9],
      [1, 2]  
    ]
"""
module BLTandPantelides

export matching, pantelides!, BLT, checkAssign
export matchingArray, pantelidesArray!, BLTArray, compactify
export expandArray, expandArrayOfArray, expandLengths!, scalarToArray, collectDerivatives, invertInvertDer
export expandStringArray
using ..BLTandPantelidesUtilities

"Controls logging"
const log = false

"""
    function augmentPath!(G, i, assign, vColour, eColour, vActive)
Construction of augmenting path

Reference:
Pantelides, C.: The consistent initialization of differential-algebraic systems. SIAM Journal
of Scientific and Statistical Computing, 9(2), pp. 213–231 (1988). 
"""
function augmentPath!(G, i, assign, vColour, eColour, vActive)
    # returns pathFound
    # assign: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
    # i: E-node
    # vActive: set to false has the same effect as deleting V-node and corresponding edges
    # j: V-node

    if log
        println("augmentPath: equation $i")
    end

    pathFound = false
    eColour[i] = true
    
    # If a V-node j exists such that edge (i-j) exists and assign[j] == 0
    for j in G[i]
        if vActive[j] && assign[j] == 0
            pathFound = true
            assign[j] = i
            return pathFound
        end
    end
  
    # For every j such that edge (i-j) exists and j is uncoloured
    for j in G[i]
        if vActive[j] && !vColour[j]
            vColour[j] = true
            k = assign[j]
            pathFound = augmentPath!(G, k, assign, vColour, eColour, vActive)
    
            if pathFound 
                assign[j] = i
                return pathFound
            end
        end
    end
    return pathFound
end


# function augmentPathArray!(G, i, assign, vColour, eColour, vActive, vLengths, eLengths)
# returns pathFound
# assign: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
# i: E-node
# vActive: set to false has the same effect as deleting V-node and corresponding edges
# j: V-node
function augmentPathArray!(G::Array{Array{Int,1},1}, i::Tuple{Int,Int}, assign::Array{Array{Tuple{Int,Int},1},1}, vColour::Array{Array{Bool,1},1}, eColour::Array{Array{Bool,1},1}, vActive::Array{Bool,1}, vLengths::Array{Int,1}, eLengths::Array{Int,1})::Bool
    if log
        println("\naugmentPath: equation $i")
    end

    pathFound = false
    eColour[i[1]][i[2]] = true
    
    # If a V-node j exists such that edge (i-j) exists and assign[j] == 0
    for j in G[i[1]]
        # for j in reverse(G[i[1]])
        if vActive[j]
            # Check for all elements of variable j
            for jj in 1:vLengths[j]
                if assign[j][jj] == (0, 0)
                    pathFound = true
                    assign[j][jj] = i
                    if log
                        println("Simply assigned: assign[$j][$jj] = $i")
                    end
                    return pathFound
                end
            end
        end
    end
  
    # For every j such that edge (i-j) exists and j is uncoloured
    for j in G[i[1]]
        # for j in reverse(G[i[1]])
        if vActive[j] 
            # Check for all elements of variable j
            for jj in 1:vLengths[j]
                if !vColour[j][jj]
                    vColour[j][jj] = true
                    k::Tuple{Int,Int} = assign[j][jj]
                    pathFound = augmentPathArray!(G, k, assign, vColour, eColour, vActive, vLengths, eLengths)
                    if pathFound 
                        assign[j][jj] = i
                        if log
                            println("Assigned: assign[$j][$jj] = $i")
                        end
                        return pathFound
                    end
                end
            end
        end
    end
    return pathFound
end


function checkAssign(assign, VSizes, VTypes, ESizes, ETypes, equationsInfix, variableNames, A, vActive=[a == 0 for a in A])
    println("Checking assignment")
    assignmentOK = true
    for j in 1:length(assign)
        if vActive[j]
            i = assign[j]
            if i > 0 && VSizes[j] != ESizes[i]
                assignmentOK = false
                print("Error: Variable ") 
                printList(variableNames, [j], A, newLine=false)
                println(" (($j)) with size=$(VSizes[j]) is assigned in equation (($i)) with size $(ESizes[i])")
            end
        end
    end

    if assignmentOK
        println("Assignment is OK")
    else
        # error("Assignment not OK")
  end
end


"""
    function matching(G, M)
Find maximum matching in bipartite graph

* `G`: bipartite graph
* `M`: number of V-nodes
* `vActive`: set to false has the same effect as deleting V-node and corresponding edges
* `return assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned

Reference:
Pantelides, C.: The consistent initialization of differential-algebraic systems. SIAM Journal
of Scientific and Statistical Computing, 9(2), pp. 213–231 (1988). 
"""
function matching(G, M, vActive=fill(true, M))
    assign::Array{Int,1} = fill(0, M)
    for i in 1:length(G)
        eColour::Array{Bool,1} = fill(false, length(G))
        vColour::Array{Bool,1} = fill(false, M)
        pathFound = augmentPath!(G, i, assign, vColour, eColour, vActive)
    end
    return assign
end

function matchingArray(G::Array{Array{Int,1},1}, vActive::Array{Bool,1}, vLengths::Array{Int,1}, eLengths::Array{Int,1})
    local assign::Array{Array{Tuple{Int,Int},1},1} 
 
    assign = fillArrayOfArray((0, 0), vLengths)
    for i in 1:length(G)
        eColour::Array{Array{Bool,1},1} = fillArrayOfArray(false, eLengths)
        vColour::Array{Array{Bool,1},1} = fillArrayOfArray(false, vLengths)
        for ii in 1:eLengths[i]
            pathFound = augmentPathArray!(G, (i, ii), assign, vColour, eColour, vActive, vLengths, eLengths)
        end
    end
    return assign
end


# -------------------------------------------------------

"""
    function pantelides!(G, M, A)
Perform index reduction with Pantelides algorithm.

* `G`: bipartite graph (updated)
* `M`: number of V-nodes
* `A`: A[j] = if V[k] = der(V[j]) then k else 0
* `return assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `return A`: A[j] = if V[k] = der(V[j]) then k else 0
* `return B`: B[i] = if E[l] = der(E[i]) then l else 0

Reference:
Pantelides, C.: The consistent initialization of differential-algebraic systems. SIAM Journal
of Scientific and Statistical Computing, 9(2), pp. 213–231 (1988). 
"""
function pantelides!(G, M, A)
    assign::Array{Int,1} = fill(0, M)
    B::Array{Int,1} = fill(0, length(G))
    N = length(G)
    N2 = N
    for k in 1:N2
        pathFound = false
        i = k
        while !pathFound
            # Delete all V-nodes with A[.] != 0 and all their incidence edges from the graph
            vActive::Array{Bool,1} = [a == 0 for a in A]
            # Designate all nodes as "uncoloured"
            eColour::Array{Bool,1} = fill(false, length(G))
            vColour::Array{Bool,1} = fill(false, M)
            pathFound = augmentPath!(G, i, assign, vColour, eColour, vActive)
            if !pathFound
                if log
                    println("\nDifferentiate:")
                end
                
                # For every coloured V-node j do
                for j in 1:length(vColour)
                    if vColour[j]
                        M += 1
                        if log
                            println("New variable derivative: var($M) = der($j)")
                        end
                        push!(A, 0)
                        A[j] = M
                        push!(assign, 0)
                    end
                end
        
                # For every coloured E-node l do
                for l in 1:N
                    if eColour[l]
                        N += 1
                        if log
                            println("New equation derivative:  equ($N) = DER($l)")
                        end
                        
                        # Create new E-node N
                        push!(G, copy(G[l]))
            
                        # Create edges from E-node N to all V-nodes j and A[j] such that edge (l-j) exists
                        for m in 1:length(G[l])
                            j = G[l][m]
                            if !(A[j] in G[N])
                                push!(G[N], A[j])
                            end
                        end
                        push!(B, 0)
                        
                        # Set B[l] = N
                        B[l] = N
                    end
                end

                # For every coloured V-node j
                for j in 1:length(vColour)
                    if vColour[j]
                        if log
                            println("Assigning derivative of variable $(A[j]) to derivative of equation: $(B[assign[j]])")
                        end
                        assign[A[j]] = B[assign[j]]
                    end
                end
                i = B[i]
            end
        end
    end

    return assign, A, B
end

function fillArrayOfArray(val, Sizes)
    a = Array{Any,1}[]
    for s in Sizes
        push!(a, fill(val, s))
    end
    return a
end
  
function expandArray(a, lengths)
    offsets = cumsum(lengths) # Offsets in expanded array for each unexpanded variable
    aa = []
    scalarToArray = [] #fill((), offsets[end])  # Not needed?
    for v in 1:length(a)
        d = a[v]
        for i in 1:lengths[v]
            if d != 0
                push!(aa, offsets[d] - lengths[d] + i)
            else
                push!(aa, 0)
            end
            push!(scalarToArray, v)
        end
    end
    
    # @show aa
    # @show scalarToArray
    return aa #, scalarToArray
end

function expandStringArray(a, lengths, par=false)
    aa = []
    for v in 1:length(a)
        for i in 1:lengths[v]
            if par
                push!(aa, "( " * a[v] * " ) [" * string(i) * "]")
            else
                push!(aa, a[v] * "[" * string(i) * "]")    
            end
        end
    end
    return aa
end

function expandArrayOfArray(G, vLengths, eLengths)
    vOffsets = cumsum(vLengths)
    GG = []
    scalarToArrayEquation = []
    for e in 1:length(G)
        g = G[e]
        gg = []
    
        for v in g
            for i in 1:vLengths[v]
                push!(gg, vOffsets[v] - vLengths[v] + i)
            end
        end
    
        for j in 1:eLengths[e]
            push!(GG, gg)
            push!(scalarToArrayEquation, (e, j))
        end
    end
    
    return GG # , scalarToArrayEquation
end

function expandLengths!(lengths, AB)
    (orgIndex, derOrder) = invertDer(AB)
    # @show lengths
    # @show AB
    # @show orgIndex derOrder
    
    for i in length(lengths) + 1:length(AB)
        # @show i
        org = orgIndex[i]
        # @show org
        push!(lengths, lengths[scalarToArray(lengths)[org]])
    end
end

function scalarToArray(lengths)
    s2A = []
    for i in 1:length(lengths)
        for j in 1:lengths[i]
            push!(s2A, i)
        end
    end
    return s2A
end

function collectDerivatives(AABB, lengths)
    (orgIndex, derOrder) = invertDer(AABB)
    @show AABB lengths
    @show orgIndex derOrder
    n = count(derOrder .== 0)
    orgInd = orgIndex[n + 1:end]
    derOrd = derOrder[n + 1:end]
    p = sortperm(orgInd)
    orgInd2 = orgInd[p]
    derOrd2 = derOrd[p]
    p2 = sortperm(derOrd2)
    derOrd3 = derOrd2[p2]
    orgInd3 = orgInd2[p2]
    orgIndex2 = [1:n; orgInd3]
    derOrder2 = [fill(0, n); derOrd3]

    @show orgIndex2 derOrder2
    return orgIndex2, derOrder2
end



function collectDerivativesOld(AABB, lengths)
    # Collect derivative variables or equations together
    # AABB updated
    AB = [] # New A or B array for array variables or equations
    scalarToArray = fill(0, length(AABB))
    k = 0 # Position in AB array
    
    # Sample last array element from the AABB list
    # since if one element is differentiated, all other elements are also
    for i in cumsum(lengths)
        index = i
        aabb = AABB[i]
        if aabb == 0
            k += 1
            push!(AB, 0)
            scalarToArray[index] = k
        elseif aabb > 0
            # Place variable and all its derivatives after each other
            while aabb > 0
                k += 1
                push!(AB, k + 1)
                scalarToArray[index] = k
                # (AABB[aabb], aabb) = (-1, AABB[aabb])  # Mark variable as already collected
                index = aabb
                aabb = AABB[aabb]
                AABB[index] = -1
            end

            k += 1
            push!(AB, 0)
            scalarToArray[index] = k
        end
    end
    return AB, scalarToArray
end


"""
    function pantelidesArray!(G, M, A, vLengths, eLengths)
Perform index reduction for array equations and variables with Pantelides algorithm.

* `G`: bipartite graph (updated)
* `M`: number of V-nodes
* `A`: A[j] = if V[k] = der(V[j]) then k else 0 (updated)
* `vLengths`: number of elements of the (array) variables, i.e. product of the dimensions (updated)
* `eLengths`: number of elements of the (array) equations, i.e. product of the dimensions (updated)
* `return assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `return A`: A[j] = if V[k] = der(V[j]) then k else 0
* `return B`: B[i] = if E[l] = der(E[i]) then l else 0
"""
function pantelidesArray!(G::Array{Array{Int,1},1}, M::Int, A::Array{Int,1}, vLengths::Array{Int,1}, eLengths::Array{Int,1})
    local assign::Array{Array{Tuple{Int,Int},1},1} 
    local B::Array{Int,1} 
  
    assign = fillArrayOfArray((0, 0), vLengths)
    B = fill(0, length(G))
    N = length(G)
    N2 = N
    for k in 1:N2  ### Gives two split variables
        # for k in N2:-1:1
        for ii in 1:eLengths[k]
            # for ii in eLengths[k]:-1:1
            pathFound = false
            i = k
            while !pathFound
                # Delete all V-nodes with A[.] != 0 and all their incidence edges from the graph
                vActive::Array{Bool,1} = [a == 0 for a in A]
                
                # Designate all nodes as "uncoloured"
                eColour::Array{Array{Bool,1},1} = fillArrayOfArray(false, eLengths)
                vColour::Array{Array{Bool,1},1} = fillArrayOfArray(false, vLengths)
                pathFound = augmentPathArray!(G, (i, ii), assign, vColour, eColour, vActive, vLengths, eLengths)
                if !pathFound
                    if log
                        println("\nDifferentiate:")
                    end
                    # For every coloured V-node j do
                    for j in 1:length(vColour)
                        for jj in 1:vLengths[j]
                            if A[j] == 0 && vColour[j][jj]
                                # if A[j] == 0 && any(vColour[j])
                                M += 1
                                if log
                                    println("New variable derivative: var($M) = der($j)")
                                end
                                push!(A, 0)
                                A[j] = M
                                push!(vLengths, vLengths[j]) # The new differentiated variable has the same length as the original.
                                push!(assign, fill((0, 0), vLengths[j]))
                            end
                        end
                    end

                    # For every coloured E-node l do
                    for l in 1:N
                        for ll in 1:eLengths[l]
                            if B[l] == 0 && eColour[l][ll]
                                # if B[l] == 0 && any(eColour[l])
                                N += 1
                                if log
                                    println("New equation derivative:  equ($N) = DER($l)")
                                end
                                
                                # Create new E-node N
                                push!(G, copy(G[l]))
                                # push!(G, [])  # Problem with wrong sorting was due to this temporary test :-(
                                
                                # Create edges from E-node N to all V-nodes j and A[j] such that edge (l-j) exists
                                for m in 1:length(G[l])
                                    j = G[l][m]
                                    
                                    # if ! (A[j] in G[N])
                                    if !(A[j] in G[l])  ### Strange difference to function pantelides!
                                        push!(G[N], A[j])
                                    end
                                end
                                push!(B, 0)
                
                                # Set B[l] = N
                                B[l] = N
                                push!(eLengths, eLengths[l]) # The new differentiated equation has the same length as the original.
                            end
                        end
                    end
          
                    # For every coloured V-node j
                    # Assign each variable derivative to the differential of the equation to which the variable was originally assigned
                    for j in 1:length(vColour)
                        for jj in 1:vLengths[j]
                            if vColour[j][jj]
                                # if any(vColour[j])
                                # assign[A[j]] = B[assign[j]]   # Scalar version
                                assign[A[j]][jj] = (B[assign[j][jj][1]], assign[j][jj][2]) # Same element of the derivative array as in the original variable is assigned.
                                if log
                                    println("Assigning derivative of variable ($(A[j]), $jj) to derivative of equation: $((B[assign[j][jj][1]], assign[j][jj][2]))")
                                end
                            end
                        end
                    end
                    i = B[i]
                end
            end
        end
    end
    return assign, A, B
end

const notOnStack = 1000000000

"""
Find minimal systems of equations that have to be solved simultaneously.
Reference: 
Tarjan, R. E. (1972), "Depth-first search and linear graph algorithms", SIAM Journal on Computing 1 (2): 146–160, doi:10.1137/0201010
"""
function strongConnect(G, assign, v, nextnode, stack, components, lowlink, number)
    # println("strongConnect: ", v)
  
    if v == 0 
        return nextnode
    end
    
    nextnode += 1
    lowlink[v] = number[v] = nextnode
    push!(stack, v)

    for w in [assign[j] for j in G[v]] # for w in the adjacency list of v
        if w > 0   # Is assigned
            if number[w] == 0 # if not yet numbered
                nextnode = strongConnect(G, assign, w, nextnode, stack, components, lowlink, number)
                lowlink[v] = min(lowlink[v], lowlink[w])
            else
                if number[w] < number[v]
                    # (v, w) is a frond or cross-link
                    # if w is on the stack of points. Always valid since otherwise number[w]=notOnStack (a big number)
                    lowlink[v] = min(lowlink[v], number[w])
                end
            end
        end
    end
  
    if lowlink[v] == number[v]
        # v is the root of a component
        # start a new strongly connected component
        # println("start a new strongly connected component")
        # @show v lowlink[v] number[v]
        # @show stack 
        # @show number[v]
        comp = []
        repeat = true
        while repeat
            # delete w from point stack and put w in the current component
            # println("delete w from point stack and put w in the current component")
            w = pop!(stack)     
            # @show w number[w]
            number[w] = notOnStack
            push!(comp, w)
            # @show comp
            repeat = w != v
        end 
        push!(components, comp)
    end
    return nextnode
end

const notOnStack = 1000000000

"""
    function strongConnectArray(G, assign, v, nextnode, stack, components, lowlink, number, vLengths, eLengths)

Find minimal systems of equations that have to be solved simultaneously.
Reference: 
Tarjan, R. E. (1972), "Depth-first search and linear graph algorithms", SIAM Journal on Computing 1 (2): 146–160, doi:10.1137/0201010
"""
function strongConnectArray(G::Array{Array{Int,1},1}, assign::Array{Array{Tuple{Int,Int},1},1}, v::Tuple{Int,Int}, nextnode::Int, stack::Array{Tuple{Int,Int},1}, components::Array{Array{Tuple{Int,Int},1},1}, lowlink::Array{Array{Int,1},1}, number::Array{Array{Int,1},1}, vLengths::Array{Int,1}, eLengths::Array{Int,1})

    #  println("\nstrongConnect: ", v)
  
    (i, ii) = v
    # @show v, nextnode
    if v == (0, 0) # 0   ### When can this happen?
        return nextnode
    end
    
    nextnode += 1
    lowlink[i][ii] = number[i][ii] = nextnode
    push!(stack, v)

    for w in [assign[j][jj] for j in G[i] for jj in 1:vLengths[j]] # for w in the adjacency list of v
        if w != (0, 0)   # Is assigned
            k = w[1]
            kk = w[2]
            if number[k][kk] == 0 # if not yet numbered
                nextnode = strongConnectArray(G, assign, w, nextnode, stack, components, lowlink, number, vLengths, eLengths)
                lowlink[i][ii] = min(lowlink[i][ii], lowlink[k][kk])
            else
                if number[k][kk] < number[i][ii]
                    # (v, w) is a frond or cross-link
                    # if w is on the stack of points. Always valid since otherwise number[w]=notOnStack (a big number)
                    lowlink[i][ii] = min(lowlink[i][ii], number[k][kk])
                end
            end
        end
    end
  
    if lowlink[i][ii] == number[i][ii]
        # v is the root of a component
        # start a new strongly connected component
        comp = Tuple{Int,Int}[]
        repeat = true
        while repeat
            # delete w from point stack and put w in the current component
            w = pop!(stack) 
            (k, kk) = w
            number[k][kk] = notOnStack
            push!(comp, w)
            if log
                println("Added equation $w to $comp")
            end
            repeat = w != v
        end 
        push!(components, comp)
        if log
            println("Added component $comp to $components")
        end
    end
    return nextnode
end

"""
    function BLT(G, assign)

Find Block Lower Triangular structure for a bipartite graph `G` with assignment `assign`
    
* `G`: bipartite graph
* `assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `return components`: cell array of components. Each component is a list of indices to E-nodes
"""
function BLT(G, assign)
    nextnode::Int = 0
    stack = []
    components = []
    lowlink = fill(0, length(G))
    number = fill(0, length(G))
  
    for v in 1:length(G)
        if number[v] == 0
            nextnode = strongConnect(G, assign, v, nextnode, stack, components, lowlink, number)
        end    
    end
    return components
end

"""
    function BLT(G, assign)

Find Block Lower Triangular structure for a bipartite graph `G` with assignment `assign`
    
* `G`: bipartite graph
* `assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `return components`: cell array of components. Each component is a list of indices to E-nodes
"""
function BLTArray(G::Array{Array{Int,1},1}, assign::Array{Array{Tuple{Int,Int},1},1}, vLengths::Array{Int,1}, eLengths::Array{Int,1})
    nextnode = 0
    stack = Tuple{Int,Int}[]
    components = Array{Tuple{Int,Int},1}[]
    lowlink::Array{Array{Int,1},1} = fillArrayOfArray(0, eLengths)
    number::Array{Array{Int,1},1} = fillArrayOfArray(0, eLengths)
  
    for i in 1:length(G)
        for ii in 1:eLengths[i]
            if number[i][ii] == 0
                nextnode = strongConnectArray(G, assign, (i, ii), nextnode, stack, components, lowlink, number, vLengths, eLengths)
            end
        end    
    end
    return components
end

function compactify(components, assign, vLengths, eLengths, vNames, equationsInfix)
    # @show vLengths eLengths
    componentsArray = [unique([e[1] for e in c]) for c in components] 
    # @show componentsArray
    
    # println("\nVerify that all array equation elements are within the same strongly connected component")
    for c in components
        es = [ee[1] for ee in c] # Collect all array equations involved in c
        for e in c
            if count(e[1] .== es) != eLengths[e[1]] # For each array equation in c, make sure all elements are involved in c.
                println("\nThe equation $e:")
                println(equationsInfix[e[1]])
                println("is split into several strongly connected components:") 
                println(c) 
            end        
        end
    end
    compArr = collect(Iterators.flatten(componentsArray))
    # @show size(compArr)
    # @show size(unique(compArr))
    @assert size(compArr) == size(unique(compArr)) "Error: array split into different strongly connected components." 

    # Collect all variables assigned in each component
    invAssign = zeros(Int, length(eLengths))
    assignArray = zeros(Int, length(vLengths))
    for c in components
        if length(c) > 9 ###
            println("Block")
        end
        
        # @show c
        # @show length(c)
        vars = []
        for e in c
            # @show e
            # println(equationsInfix[e[1]])
            
            # Search for all variables assigned in all equations of the component
            for i in 1:length(assign)
                for j in 1:length(assign[i])
                    # @show e[1] assign[i][j]
                    if e == assign[i][j]   # All equation elements are within the same component.
                        push!(vars, i)
                        # print("push ", vNames[i])
                        if length(assign[i]) > 1
                            # print("[$j]")
                        end
                        # println()
                    end
                end
            end
        end

        if length(c) > 9 ###
            println("end block")
        end
        # @show length(vars)
        vars = unique(vars) # Remove duplicates
        # @show length(vars)
        # @show vars
        cc = unique([e[1] for e in c])
        # @show cc
        # @show length(cc) length(vars)
        # @show sum(eLengths[e] for e in cc)
        # @show sum(vLengths[v] for v in vars)
        # assert(length(cc) == length(vars))
        scalarV = 0
        scalarE = 0
        for i in 1:length(vars)
            # @show cc[i] vars[i]
            invAssign[cc[i]] = vars[i] 
            assignArray[vars[i]] = cc[i] 
            scalarV += vLengths[vars[i]]
            scalarE += eLengths[cc[i]]
        end
        # @show scalarV scalarE 
        if scalarV != scalarE
            for i in 1:length(vars)
                println(vLengths[vars[i]], " ", vNames[vars[i]])
            end
            for i in 1:length(vars)
                println(eLengths[cc[i]], " ", equationsInfix[cc[i]])
            end
        # assert(false)
        end
    end
    # @show invAssign
  
    for e in 1:length(invAssign)
        if invAssign[e] > 0
            # println("equation ", e, "  variable: ", invAssign[e], " ", vNames[invAssign[e]])  
        else
            println("Error: equation ", e, " not assigned")        
            # assert(false)
        end
    end
    return invAssign, assignArray, componentsArray
end


# Could not use using BLTandPantelidesUtilities!

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

function invertInvertDer(orgIndex, derOrder)
    A = fill(0, length(orgIndex))
    for i in 1:length(orgIndex)
        if derOrder[i] > 0
            e = orgIndex[i]
            for d in 1:derOrder[i] - 1
                e = A[e]
            end
            A[e] = i
        end
    end
    return A
end

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

function printAssignedEquations(equations, variables, indices, assign, A, B)
    (orgIndexVar, derOrderVar) = invertDer(A)
    (orgIndexEqu, derOrderEqu) = invertDer(B)
    (assignedVar,) = invertAssign(assign)
    for i in indices
        j = assignedVar[i]
        if j > 0    
            if derOrderVar[j] == 1
                print("der(")
            elseif derOrderVar[j] > 1
                print("der", derOrderVar[j], "(")
            end
            print(variables[orgIndexVar[j]])
            if derOrderVar[j] > 0
                print(")")
            end
            print("   (($j))   :=   ")
        end
      
        if derOrderEqu[i] == 1
            print("DER()   ")
        elseif derOrderEqu[i] > 1
            print("DER", derOrderEqu[i], "()   ")
        end
        println(equations[orgIndexEqu[i]], "   (($i))")
    end
end

function printList(infixes, indices, A, vertical=false; newLine=true)
    (orgIndex, derOrder) = invertDer(A)
    for ind in 1:length(indices)
        j = indices[ind]
        if j > 0    
            if ind > 1
                if vertical
                    println()
                else
                    print(", ")
                end
            end
            
            if derOrder[j] > 0
                if vertical
                    print("DER")
                else
                    print("der")
                end
            
                if derOrder[j] > 1
                    print(derOrder[j])
                end
                print("(")
            end
            
            print(infixes[orgIndex[j]])
            
            if derOrder[j] > 0
                print(")")
            end
        end
    end
    
    if newLine
        println()
    end
end

end