#=
Entry point of Modia

* Developer: Hilding Elmqvist, Mogram AB
* First version: December 2020
* License: MIT (expat)

=#

function assert(condition, message)
    if ! condition
        println("ASSERTION: ", message)
    end
end

function printArray(array, heading; order=1:length(array), log=false, number=true, indent=false, spacing=" ")
    if length(array) > 0 && log
        println(heading)
        n = 0
        for v in array
            n += 1
            if number
                print(lpad(order[n], 4), ": ")
            end
            if indent
                print(("|"*spacing)^(n-1))
            end
            println(replace(replace(string(v), "\n" => " "), "  " => " "))
        end
    end
end


"""
    newName = moveDerToLeaf(name::String)::String

Move `der(..)` to leaf name part.

# Examples
```
moveDerToLeaf("a.b.c")                     # = "a.b.c"
moveDerToLeaf("der(x)")                    # = "der(x)"
moveDerToLeaf("der(der(x)")                # = "der(der(x))"
moveDerToLeaf("der(a.x)")                  # = "a.der(x)"
moveDerToLeaf("der(a.b.x)")                # = "a.b.der(x)"
moveDerToLeaf("der(der(der(a.b.c.x)))")    # = "a.b.c.der(der(der(x)))"
```
"""
function moveDerToLeaf(name::String)::String
    if length(name) >= 4 && name[1:4] == "der(" && !isnothing(findfirst(".", name))
        # name starts with der(..) and name has a dot (".")
        i = 5
        indexRange_old = 1:4
        while true
            indexRange = findnext("der(", name, i)
            if isnothing(indexRange)
                j1 = indexRange_old.stop
                j2 = findlast(".", name).stop
                newName = name[j1+1:j2] * name[1:j1] * name[j2+1:end]
                return newName
            end
            i = indexRange[2]+1
            indexRange_old = indexRange
        end
        @error("moveDerToLeaf($name): Error should not occur")
    end
    return name
end


function performConsistencyCheck(G, Avar, vActive, parameters, unknowns, states, equations, log=false)
    nUnknowns = length(unknowns) - length(states)
    nEquations = length(equations)

    if nEquations > 0 && nUnknowns != nEquations
        printstyled("The number of unknowns ($nUnknowns) and equations ($nEquations) is not equal!\n", bold=true, color=:red)
        printArray(parameters, "Parameters:", log=true)
        printArray(states, "Potential states:", log=true)
        printArray(setdiff(unknowns,states), "Unknowns:", log=true)
        printArray(equations, "Equations:", log=true)
    end

    ok = nUnknowns == nEquations

    # Check consistency
    EG::Array{Array{Int,1},1} = [G; buildExtendedSystem(Avar)]
    assignEG = matching(EG, length(Avar))
#        @assert all(assignEG .> 0) "  The DAE is structurally singular\n  $assignEG."
    assigned = [a > 0 for a in assignEG]
    if ! (all(assigned)) || ! ok
        if nUnknowns > nEquations
            missingEquations = nUnknowns - nEquations
            printstyled("\nThe DAE has no unique solution, because $missingEquations equation(s) missing!\n", bold=true, color=:red)
        elseif nEquations > nUnknowns
            tooManyEquations = nEquations - nUnknowns
            printstyled("\nThe DAE has no unique solution, because $tooManyEquations equation(s) too many\n", bold=true, color=:red)
        else
            printstyled("\nThe DAE has no unique solution, because it is structurally inconsistent\n" *
                        "  (the number of equations is identical to the number of unknowns)\n", bold=true, color=:red)
        end

        singular = true
        hequ = []
        # Create dummy equations for printing: integrate(x, der(x)) = 0
        for i in 1:length(Avar)
            a = Avar[i]
            if a > 0
                n = unknowns[i]
                der = unknowns[a]
                push!(hequ, "integrate($n, $der) = 0")
            end
        end
        equationsEG = vcat(equations, hequ)
        println("\nAssigned equations after consistency check:")
        printAssignedEquations(equationsEG, unknowns, 1:length(EG), assignEG, Avar, fill(0, length(EG)))
        unassigned = []
        for i in eachindex(assignEG)
            if assignEG[i] == 0
                push!(unassigned, unknowns[i])
            end
        end
        printArray(unassigned, "Not assigned variables after consistency check:", log=true)
        if nUnknowns > nEquations
            missingEquations = nUnknowns - nEquations
            error("Aborting, because $missingEquations equation(s) missing!")
        elseif nEquations > nUnknowns
            tooManyEquations = nEquations - nUnknowns
            error("Aborting, because $tooManyEquations equation(s) too many!")
        else
            error("Aborting, because DAE is structurally inconsistent")
        end
    elseif log
        println("  The DAE is structurally nonsingular.")
    end
end



function assignAndBLT(model, modelName, modelModule, equations, unknowns, modelStructure, Avar, G, states, logDetails, log=false, logTiming=false)
    parameters = modelStructure.parameters
    vActive = [a == 0 for a in Avar]

    if drawIncidence
        showIncidence(equations, unknowns, G, vActive, [], [], [], [[e] for e in 1:length(G)])
    end

    if false # logAssignmentBeforePantelides
        @show G length(Avar) vActive
        assign = matching(G, length(Avar), vActive)
        @show assign

        println("\nAssigned equations:")
        Bequ = fill(0, length(G))
        printAssignedEquations(equations, unknowns, 1:length(G), assign, Avar, Bequ)

        printUnassigned(equations, unknowns, assign, Avar, Bequ, vActive)
    end

    if logTiming
        println("\nCheck consistency of equations by matching extended equation set")
        @timeit to "performConsistencyCheck" performConsistencyCheck(G, Avar, vActive, parameters, unknowns, states, equations, log)
    else
        performConsistencyCheck(G, Avar, vActive, parameters, unknowns, states, equations, log)
    end

    if logTiming
        println("\nIndex reduction (Pantelides)")
        @timeit to "pantelides!" assign, Avar, Bequ = pantelides!(G, length(Avar), Avar)
    else
        assign, Avar, Bequ = pantelides!(G, length(Avar), Avar)
    end
    if logDetails
        @show assign Avar Bequ
    end

    moreVariables = length(Avar) > length(unknowns)
    unknowns = vcat(unknowns, fill(:x, length(Avar) - length(unknowns)))
    for i in 1:length(Avar)
        if Avar[i] > 0
            unknowns[Avar[i]] = derivative(unknowns[i], keys(parameters))
        end
    end
    vActive = [a == 0 for a in Avar]

    moreEquations = length(Bequ) > length(equations)
    equations = vcat(equations, fill(:(a=b), length(Bequ) - length(equations)))
    for i in 1:length(Bequ)
        if Bequ[i] > 0
            if log
                println("Differentiating: ", equations[i])
            end
            equations[Bequ[i]] = simplify(derivative(equations[i], keys(parameters)))
            if log
                println("  Differentiated equation: ", equations[Bequ[i]])
            end
        end
    end

    if moreVariables && log
        printArray(unknowns, "Unknowns after index reduction:", log=logDetails)
    end
    if moreEquations && log
        printArray(equations, "Equations after index reduction:", log=logDetails)
    end

    HG = Array{Array{Int64,1},1}()
    originalEquationIndex = []
    highestDerivativeEquations = []
    for i in 1:length(G)
        if Bequ[i] == 0
            push!(HG, G[i])
            push!(highestDerivativeEquations, equations[i])
            push!(originalEquationIndex, i)
        end
    end

    if logTiming
        println("Assign")
        @timeit to "matching" assign = matching(HG, length(Avar), vActive)
    else
        assign = matching(HG, length(Avar), vActive)
    end
    if logTiming
        println("BLT")
        @timeit to "BLT" bltComponents = BLT(HG, assign)
    else
       bltComponents = BLT(HG, assign)
    end
    if logDetails
        @show HG
        @show bltComponents
    end

    if log
        println("\nSorted highest derivative equations:")
        printSortedEquations(highestDerivativeEquations, unknowns, bltComponents, assign, Avar, Bequ)
    end

    if experimentalTranslation
        return specialTranslation(model, modelName, modelModule, highestDerivativeEquations, unknowns, modelStructure, states, bltComponents, assign, Avar, Bequ)
    end


    if drawIncidence
        assignedVar, unAssigned = invertAssign(assign)
        showIncidence(equations, unknowns, HG, vActive, [], assign, assignedVar, bltComponents, sortVariables=true)
    end

    blt = Vector{Int}[]  # TODO: Change type
    for c in bltComponents
        push!(blt, originalEquationIndex[c])
    end
    for i in eachindex(assign)
        if assign[i] != 0
            assign[i] = originalEquationIndex[assign[i]]
        end
    end
    return (unknowns, equations, G, Avar, Bequ, assign, blt, parameters)
end

function setAvar(unknowns)
    variablesIndices = OrderedDict(key => k for (k, key) in enumerate(unknowns))
    Avar = fill(0, length(unknowns))
    states = []
    derivatives = []
    for i in 1:length(unknowns)
        v = unknowns[i]
        if isexpr(v, :call)
            push!(derivatives, v)
            x = v.args[2]
            ix = variablesIndices[x]
            Avar[ix] = i
            push!(states, unknowns[ix])
        end
    end
    return Avar, states, derivatives
end

mutable struct ModelStructure
    parameters::OrderedDict{Any,Any}
    mappedParameters::OrderedDict{Symbol,Any}
    init::OrderedDict{Any,Any}
    start::OrderedDict{Any,Any}
    variables::OrderedDict{Any,Any}
    flows::OrderedDict{Any,Any}
    inputs::OrderedDict{Any,Any}
    outputs::OrderedDict{Any,Any}
#    components
#    extends
    equations::Array{Expr,1}
    hideResults::OrderedSet{Any}   # Do not store these variables in the result data structure
end

ModelStructure() = ModelStructure(OrderedDict(), OrderedDict{Symbol,Any}(), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), Expr[], OrderedSet())

function printDict(label, d)
    if length(d) > 0
        print(label)
        for (k,v) in d
            print(k)
            if v != nothing
                print(" => ", v)
            end
            print(", ")
        end
        println()
    end
end

function printModelStructure(m, label=""; log=false)
    if log
        println("ModelStructure ", label)
        printDict("  parameters:  ", m.parameters)
        printDict("  mappedParameters:  ", m.mappedParameters)
        printDict("  init:        ", m.init)
        printDict("  start:       ", m.start)
        printDict("  variables:   ", m.variables)
        printDict("  flows:       ", m.flows)
        printDict("  inputs:      ", m.inputs)
        printDict("  outputs:     ", m.outputs)
        println("  equations: ")
        for e in m.equations
            println("    ", e)
        end
    end
end

prependDict(dict, prefix) = OrderedDict([prepend(k, prefix) => prepend(v, prefix) for (k,v) in dict])


function mergeModelStructures(parent::ModelStructure, child::ModelStructure, prefix)
    if length(child.mappedParameters) > 0
        merge!(parent.parameters, child.parameters)
        parent.mappedParameters[prefix] = child.mappedParameters

        merge!(parent.init, child.init)
        parent.mappedParameters[prefix] = child.mappedParameters

        merge!(parent.start, child.start)
        parent.mappedParameters[prefix] = child.mappedParameters
    end

    merge!(parent.variables, child.variables)
    merge!(parent.flows, child.flows)
    union!(parent.hideResults, child.hideResults)

    push!(parent.equations, prepend(child.equations, prefix)...)
end

function replaceLinearIntegerEquations(Gint, eInt, GcInt, unknowns, equations)
    for i = 1:length(Gint)
        e = Gint[i]
        if length(e) > 0
            # Construct equation
            rhs = 0
            for j = 1:length(e)
                v  = e[j]
                vc = GcInt[i][j]
                rhs = add(rhs, mult(vc, unknowns[v]))
            end
            equ = :(0 = $rhs)
            equations[eInt[i]] = equ
        else
            equations[eInt[i]] = :(0 = 0)
        end
    end
end

function performAliasReduction(unknowns, equations, Avar, logDetails, log)

    @timeit to "enumerate(unknowns)" variablesIndices = OrderedDict(key => k for (k, key) in enumerate(unknowns))

    linearVariables = []
    nonlinearVariables = []
    linearEquations = Vector{Int}()
    nonlinearEquations = []
    G = Vector{Int}[] # Array{Array{Int64,1},1}()
    Gint = Array{Array{Int64,1},1}()
    GcInt = Array{Array{Int64,1},1}()

    @timeit to "build graph information" for i in 1:length(equations)
        e = equations[i]
        nameIncidence, coefficients, rest, linear = getCoefficients(e)
        incidence = [variablesIndices[n] for n in nameIncidence if n in keys(variablesIndices)]
        unique!(incidence)
        push!(G, incidence)

        if linear && all([c == 1 || c == -1 for c in coefficients]) && rest == 0 && all(n in keys(variablesIndices) for n in nameIncidence)
            push!(linearVariables, nameIncidence...)
            push!(linearEquations, i)
            push!(Gint, incidence)
            push!(GcInt, coefficients)
        else
            push!(nonlinearVariables, nameIncidence...)
            push!(nonlinearEquations, e)
        end
    end

    if logDetails
        @show G
        @show Avar
    end

    @timeit to "unique!(linearVariables)" unique!(linearVariables)
    @timeit to "unique!(nonlinearVariables)" unique!(nonlinearVariables)
    @timeit to "setdiff" vSolveInLinearEquations = setdiff(linearVariables, nonlinearVariables)
    if logDetails
        printArray(linearVariables, "linearVariables:", number=false)
        printArray(equations[linearEquations], "linearEquations:", number=false)
#        printArray(nonlinearVariables, "nonlinearVariables:", number=false)
#        printArray(nonlinearEquations, "nonlinearEquations:", number=false)
        @show linearEquations vSolveInLinearEquations Gint GcInt
    end

    GCopy = deepcopy(G)
    AvarCopy = deepcopy(Avar)
    @timeit to "simplifyLinearIntegerEquations!" (vEliminated, vProperty, nvArbitrary, redundantEquations) = simplifyLinearIntegerEquations!(GCopy, linearEquations, GcInt, AvarCopy)

    if logDetails
        @show vEliminated vProperty nvArbitrary redundantEquations
    end
    if log
        printSimplifiedLinearIntegerEquations(GCopy, linearEquations, GcInt, vEliminated, vProperty, nvArbitrary, redundantEquations, v->string(unknowns[v]), printTest=false)
    end

    substitutions = OrderedDict()
    @timeit to "build substitutions" if length(vEliminated) - nvArbitrary > 0
        for v in vEliminated[nvArbitrary+1:end]
            if isZero(vProperty, v)
                substitutions[unknowns[v]] = 0.0
            elseif isAlias(vProperty, v)
                substitutions[unknowns[v]] = unknowns[alias(vProperty,v)]
            else
                substitutions[unknowns[v]] = :(- $(unknowns[negAlias(vProperty,v)]))
            end
        end
    end

    @timeit to "replaceLinearIntegerEquations" replaceLinearIntegerEquations(GCopy[linearEquations], linearEquations, GcInt, unknowns, equations)
#    printArray(equations, "new equations:", log=log)

    reducedEquations = []
    @timeit to "substitute" for ei in 1:length(equations)
        e = equations[ei]
        if e != :(0 = 0)
            push!(reducedEquations, substitute(substitutions, e))
        end
    end

    reducedUnknowns = []
    for vi in 1:length(unknowns)
        if isNotEliminated(vProperty, vi)
            push!(reducedUnknowns, unknowns[vi])
        end
    end

    reducedVariablesIndices = OrderedDict(key => k for (k, key) in enumerate(reducedUnknowns))
    reducedG = Vector{Int}[] # Array{Array{Int64,1},1}()
    @timeit to "build reducedG" for i in 1:length(reducedEquations)
        e = reducedEquations[i]
        nameIncidence = Incidence[]
        findIncidence!(e, nameIncidence)
        incidence = [reducedVariablesIndices[n] for n in nameIncidence if n in keys(reducedVariablesIndices)]
        unique!(incidence)
        push!(reducedG, incidence)
    end

    reducedAvar, reducedStates, reducedDerivatives = setAvar(reducedUnknowns)
    return reducedEquations, reducedUnknowns, reducedAvar, reducedG, reducedStates, vEliminated, vProperty
end


function stateSelectionAndCodeGeneration(modStructure, Gexplicit, name, modelModule, buildDict, FloatType, TimeType, init, start, inputs, outputs, vEliminated, vProperty, unknownsWithEliminated, mappedParameters, hideResults;
    unitless=false, logStateSelection=false, logCode=false, logExecution=false, logCalculations=false, logTiming=false, evaluateParameters=false, saveCodeOnFile="")
    (unknowns, equations, G, Avar, Bequ, assign, blt, parameters) = modStructure
    Goriginal = deepcopy(G)
    function getSolvedEquationAST(e, v)
        (solution, solved) = solveEquation(equations[e], unknowns[v])
        unknown = solution.args[1]
        expr = solution.args[2]
###        solution = :(_variables[$(v+1)] = $unknown = $expr)
        equ = string(equations[e])
        @assert solved "Equation not solved for $(unknowns[v]): $(equations[e])"
#        sol = string(solution)
#        solution = prepend(makeDerivativeVar(solution, components), :m)
        solution = substituteForEvents(solution)
        if false
    #        solution = :(begin println("Executing: ", $sol); $solution end)
            if sol != equ
                solution = :(try $solution; catch e; println("\n\nFailure executing: ", $sol); println("Original equation: ", $equ, "\n\n"); rethrow(e) end)
            else
                solution = :(try $solution; catch e; println("\n\nFailure executing: ", $sol, "\n\n"); rethrow(e) end)
            end
    #        solution = :(try $solution; catch e; println("Failure executing: ", $sol); println(sprint(showerror,e)); rethrow(e) end)
    #        solution = :(try $solution; catch e; println("Failure executing: ", $sol); printstyled(stderr,"ERROR: ", bold=true, color=:red);
    #        printstyled(stderr,sprint(showerror,e), color=:light_red); println(stderr); end)
        end
        solution = makeDerVar(solution, parameters, inputs, FloatType, evaluateParameters)
        if logCalculations
            var = string(unknowns[v])
            solutionString = string(solution)
            #return :(println("Calculating: ", $solutionString); $solution; println("  Result: ", $var, " = ", upreferred.($(solution.args[1]))))
            return :(println("Calculating: ", $solutionString); $solution; println("  Result: ", $var, " = ", $(solution.args[1])))
        else
            return solution
        end
    end

    function getResidualEquationAST(e, residualName)
        eq = equations[e] # prepend(makeDerivativeVar(equations[e], components), :m)
        lhs = eq.args[1]
        rhs = eq.args[2]
        if isexpr(rhs, :call) && rhs.args[1] == :_DUPLICATEEQUATION
            return nothing
        end
        if isexpr(lhs, :tuple) && all(a == 0 for a in lhs.args) || lhs == :(0)
            eq_rhs = makeDerVar(:($rhs), parameters, inputs, FloatType, evaluateParameters)
            if unitless
                eqs = eq_rhs
            else
                eqs = :(Modia.Unitful.ustrip.($eq_rhs))
            end
        #elseif isexpr(lhs, :tuple) && isexpr(rhs, :call) && unitless
        #    eq_rhs = makeDerVar(:($rhs), parameters, inputs, evaluateParameters)
        #    eq_lhs = makeDerVar(:($lhs), parameters, inputs, evaluateParameters)
        #    eqs =  :( ($eq_rhs .-= $eq_lhs) )
        else
            eq_rhs = makeDerVar(:($rhs), parameters, inputs, FloatType, evaluateParameters)
            eq_lhs = makeDerVar(:($lhs), parameters, inputs, FloatType, evaluateParameters)
            if unitless
                eqs = :( $eq_rhs .- $eq_lhs )
            else
                eqs = :( Modia.Unitful.ustrip.($eq_rhs) .- Modia.Unitful.ustrip.($eq_lhs))
            end
        end
        residual = :(Modia.appendVariable!(_leq_mode.residuals, $eqs))
        residString = string(eqs)
        if logCalculations
            return :(println("Calculating residual: ", $residString); $residual; $residualName = $eqs; println("  Residual: ", $residualName) )
#            return makeDerVar(:(dump($(makeDerVar(eq.args[2]))); dump($(makeDerVar(eq.args[1]))); $residual; println($residualName, " = ", upreferred.(($(eq.args[2]) - $(eq.args[1]))))))
        else
            return residual
        end
    end

    function isLinearEquation(e_original, v_original)
        equ = equations[e_original]
        var = unknowns[v_original]
        return isLinear(equ, var)
    end

    function isSolvableEquation(e_original, v_original)
        equ = equations[e_original]
        var = unknowns[v_original]
        #(solution, solved) = solveEquation(equ, var)
        #return solved
        (rest, factor, linear) = ModiaBase.linearFactor(equ, var)
        return linear && (typeof(factor) <: Number && factor != 0)
    end

    hasParticles(value) = typeof(value) <: MonteCarloMeasurements.StaticParticles ||
                          typeof(value) <: MonteCarloMeasurements.Particles

    invAvar = ModiaBase.revertAssociation(Avar)

    function var_unit(v)
        int_v = invAvar[v]
        if int_v > 0
            # v is a derivative variable
            var = unknowns[int_v]
        else
            # v is a non-differentiated variable
            var = unknowns[v]
        end

        if var in keys(init)
            value = eval(init[var])
        elseif var in keys(start)
            value = eval(start[var])
        else
            value = 0.0
        end
        if int_v > 0
            value = value / u"s"
        end

        if ! (typeof(value) <: AbstractArray)
            un = unit(value)
        elseif length(value) > 0
            un = unit.(value)
            @assert all([un[i] == un[1] for i in 2:length(un)]) "The unit of all elements of state vector must be equal: $var::$(value)"
            un = un[1]
        else
            un = ""
        end
        return SignalTables.unitAsParseableString(un)
    end

    function var_startInitFixed(v_original)
        v = unknowns[v_original]
        value::Any  = nothing
        fixed::Bool = false
        if v in keys(init)
            value = eval(init[v])
            fixed = true
        elseif v in keys(start)
            value = eval(start[v])
        end
        return (value, fixed)
    end

    function eq_equation(e)
        return replace(replace(string(equations[e]), "\n" => " "), "  " => " ")
    end

    # ----------------------------------------------------------------------------
    solvedAST = []
    Gsolvable = copy(G)
    juliaVariables  = [Symbol(u) for u in unknowns]
    stringVariables = [String(Symbol(v)) for v in unknowns]
    allEquations    = equations
    #    @show stringVariables equations G blt assign Avar Bequ Gsolvable


    stateSelectionFunctions = StateSelectionFunctions(
        var_name               = v -> moveDerToLeaf(stringVariables[v]),
        var_julia_name         = v -> juliaVariables[v],
        var_unit               = var_unit,
        var_startInitFixed     = var_startInitFixed,
        var_nominal            = v_original -> NaN,
        var_unbounded          = v_original -> false,
        equation               = eq_equation,
        isSolvableEquation     = isSolvableEquation,
        isLinearEquation       = isLinearEquation,
        getSolvedEquationAST   = getSolvedEquationAST,
        getResidualEquationAST = getResidualEquationAST,
        showMessage            = (message;severity=0,from="???",details="",variables=Int[],equations=Int[]) ->
                                    Modia.showMessage(message, severity, from, details,
                                                      stringVariables[variables],
                                                      string.(allEquations[equations]))
        )

    if length(equations) == 0
        return nothing
    end

    if logTiming
        println("Get sorted and solved AST")
        @timeit to "getSortedAndSolvedAST" (success,AST,equationInfo) = getSortedAndSolvedAST(
            Goriginal, Gexplicit, Array{Array{Int64,1},1}(blt), assign, Avar, Bequ, stateSelectionFunctions;
            unitless, log = logStateSelection, modelName = name)
    else
        (success,AST,equationInfo) = getSortedAndSolvedAST(
            Goriginal, Gexplicit, Array{Array{Int64,1},1}(blt), assign, Avar, Bequ, stateSelectionFunctions;
            unitless, log = logStateSelection, modelName = name)
    end

    if ! success
        error("Aborting")
    end

    if useNewCodeGeneration
        AST, newFunctions = prepareGenerateGetDerivatives(AST, modelModule, functionSize, juliaVariables)
    else
        newFunctions = Vector{Expr}()
    end

#    @show AST
#    @show equationInfo

#=
    selectedStates = [v.x_name_julia for v in equationInfo.x_info]

    startValues = []
    for v in selectedStates
        v = Meta.parse(string(v))
        if v in keys(init)
            value = eval(init[v])
#            value = [0.0 0.0]  # TODO: Fix
#            @show value
            if hasParticles(value)
                push!(startValues, value)
            else
                push!(startValues, value...)
            end
        elseif v in keys(start)
            value = eval(start[v])
            if hasParticles(value)
                push!(startValues, value)
            else
                push!(startValues, value...)
            end
        else
            push!(startValues, 0.0)
        end
    end
    if logCode
        println("startValues = ", startValues)
    end
=#

    vSolvedWithInit = equationInfo.vSolvedWithFixedTrue
    vSolvedWithInitValuesAndUnit = OrderedDict{String,Any}()
    for name in vSolvedWithInit
        nameAsExpr = Meta.parse(name)
        if haskey(init, nameAsExpr)
            vSolvedWithInitValuesAndUnit[name] = eval( init[nameAsExpr] )
        else
            @warn "Internal issue of Modia: $name is assumed to have an init-value, but it is not found."
        end
    end


    # Generate code
    nCrossingFunctions, nAfter, nClocks, nSamples, previousVars, preVars, holdVars = getEventCounters()
    previousVars = Symbol.(previousVars)
    preVars = Symbol.(preVars)
    holdVars = Symbol.(holdVars)

    # Variables added to result
    timeName = :time
    w_invariant_names = setdiff([Symbol(u) for u in unknowns],
                                [Symbol(h) for h in hideResults],
                                Symbol[Symbol(xi_info.x_name_julia)     for xi_info in equationInfo.x_info],
                                Symbol[Symbol(xi_info.der_x_name_julia) for xi_info in equationInfo.x_info])

    if true # logTiming
#        println("Generate code")
        if useNewCodeGeneration
            @timeit to "generate_getDerivativesNew!" code = generate_getDerivativesNew!(AST, newFunctions, modelModule, equationInfo, [:(_p)], timeName, w_invariant_names, previousVars, preVars, holdVars, :getDerivatives, hasUnits = !unitless)
        else
            @timeit to "generate_getDerivatives!" code = generate_getDerivatives!(FloatType, TimeType, AST, equationInfo, [:(_p)], timeName, w_invariant_names, previousVars, preVars, holdVars, :getDerivatives, hasUnits = !unitless)
        end
    else
#        code = generate_getDerivatives!(AST, equationInfo, Symbol.(keys(parameters)), timeName, w_invariant_names, :getDerivatives, hasUnits = !unitless)
        code = generate_getDerivativesNew!(AST, newFunctions, modelModule, equationInfo, [:(_p)], timeName, w_invariant_names, previousVars, preVars, holdVars, :getDerivatives, hasUnits = !unitless)
    end
    if logCode
        #@show mappedParameters
        showCodeWithoutComments(code)
    end

    # Compile code

#    generatedFunction = @RuntimeGeneratedFunction(modelModule, code)
    #getDerivatives = Core.eval(modelModule, code)

    if saveCodeOnFile != ""
        println("  Save generated code on file ", joinpath(pwd(), saveCodeOnFile))
        open(saveCodeOnFile, "w") do io
            print(io, replace(sprint(show,code), r"( *#= .*=#\n)|(#= .*=#)" => ""))
        end
    end

    if logTiming
        println("eval code")
        @time @timeit to "eval(code)" Core.eval(modelModule, code)
    else
        Core.eval(modelModule, code)
    end
    getDerivatives = modelModule.getDerivatives

    # If generatedFunction is not packed inside a function, DifferentialEquations.jl crashes
#    getDerivatives(derx,x,m,time) = generatedFunction(derx, x, m, time)

#    convertedStartValues = convert(Vector{FloatType}, [ustrip(v) for v in startValues])  # ustrip.(value) does not work for MonteCarloMeasurements
#    @show mappedParameters

#    println("Build InstantiatedModel")
    hideResult_names = [string(h) for h in hideResults]
    model = @timeit to "build InstantiatedModel" InstantiatedModel{FloatType,TimeType}(modelModule, name, buildDict, getDerivatives, equationInfo, previousVars, preVars, holdVars,
                                         mappedParameters, timeName, w_invariant_names, hideResult_names;
                                         vSolvedWithInitValuesAndUnit, vEliminated, vProperty,
                                         var_name = (v)->string(unknownsWithEliminated[v]),
                                         nz=nCrossingFunctions, nAfter=nAfter,  unitless=unitless)

    # Execute code
    if logExecution
        println("\nExecute getDerivatives")
    #    @show startValues
        println("First executions of getDerivatives")
        @timeit to "execute getDerivatives" try
            #@time Base.invokelatest(getDerivatives, derx, x_startValues, model, convert(TimeType, 0.0))
            #@time Base.invokelatest(getDerivatives, derx, x_startValues, model, convert(TimeType, 0.0))
            @time init!(model)  # getDerivatives is called
            #@time invokelatest_getDerivatives_without_der_x!(x_startValues, model, convert(TimeType, 0.0))
            #@time invokelatest_getDerivatives_without_der_x!(x_startValues, model, convert(TimeType, 0.0))
#           @show derx
        catch e
            error("Failed: ", e)
            return nothing
        end
    end
    return model
end


function duplicateMultiReturningEquations!(equations)
    duplicatedEquations = []
    for e in equations
        if isexpr(e, :(=)) && isexpr(e.args[1], :tuple) || typeof(e.args[1]) <: NTuple
            lhs = e.args[1]
            if typeof(lhs) <: NTuple
                nargs = length(lhs)
            else
                nargs = length(lhs.args)
            end

            if e.args[2].args[1] != :implicitDependency
                func = :_DUPLICATEEQUATION
            else
                func = :_DUPLICATIMPLICITDEPENDENCY
            end

            nameIncidence = Incidence[]
            findIncidence!(e, nameIncidence, false)
            bandWidth = length(nameIncidence)-nargs+2
            #@show length(nameIncidence) bandWidth
            for i in 1:(nargs-1)
                if false
                    rhs = Expr(:call, :_DUPLICATEEQUATION, nameIncidence...)
                else
                    indices =[]
                    for j in 1:bandWidth
                        push!(indices, mod(j+i-1,length(nameIncidence))+1)
                    end
                    rhs = Expr(:call, :_DUPLICATEEQUATION, nameIncidence[indices]...)
                end
                newE = :(0 = $rhs)
                push!(duplicatedEquations, newE)
            end
        end
    end
    append!(equations, duplicatedEquations)
end


appendSymbol(path::Nothing, name::Symbol) = name
appendSymbol(path         , name::Symbol) = :( $path.$name )


"""
    modifiedModel = buildSubmodels!(model, modelModule, FloatType, TimeType,
                                    instantiateModelOptions, buildDict::AbstractDict)

Traverse `model` and for every `<submodel>` that is a `Model(..)` and has a key-value pair
`:_buildFunction = Par(functionName = <buildFunction>)` and optionally `:_buildOption=<buildOption>`, call

```
updatedSubmodel = <buildFunction>(submodel, modelModule, FloatType::Type, TimeType::Type, instantiateModelOptions,
                      ID, pathAST::Union{Expr,Symbol,Nothing}, buildOption = <buildOption>)
```

A new `updatedSubmodel` is generated from `submodel` merged with additional code and then returned.
The arguments of `<buildFunction>`are:

- `updatedSubmodel`: A potentially new reference to the updated `submodel`
- `modelModule`: Module in which the model is defined
- `FloatType`, `TimeType`: Types used when instantiating `InstantiatedModel{FloatType,TimeType}`
- `instantiateModelOptions`: Optional arguments of `@instantiateModel` provided as `OrderedDict{Symbol,Any}`.
- `ID`: Unique ID to identify the generated submodel (to be used in the code merged into the submodel)
- `pathAST`: Path upto `<submodel>` as Abstract Syntax Tree, such as: `:( a.b.c )`
             (this path might be used as part of a variable name in the code merged into the submodel).
- `buildOption`: Option used for the generation of `buildCode`.

Note, keys `_buildFunction` and `_buildOption` have been deleted in the returned `updatedSubmodel`.
"""
function buildSubmodels!(model::AbstractDict, modelModule, FloatType::Type, TimeType::Type, instantiateModelOptions::OrderedDict{Symbol,Any},
                         buildDict::OrderedDict{String,Any}; pathAST::Union{Expr,Symbol,Nothing}=nothing)
    if haskey(model, :_buildFunction)
        _buildFunction = model[:_buildFunction]
        if haskey(_buildFunction, :functionName)
            buildFunction = _buildFunction[:functionName]
        else
            @error "Model $pathAST has key :_buildFunction but its value has no key :functionName"
        end
        delete!(model, :_buildFunction)
        ID = modelPathAsString(pathAST)
        quotedPathAST = Meta.quot(pathAST)
        if haskey(model, :_buildOption)
            buildOption = model[:_buildOption]
            delete!(model, :_buildOption)
            (model, instantiatedSubmodelStruct) = Core.eval(modelModule, :($buildFunction($model, $modelModule, $FloatType, $TimeType, $instantiateModelOptions, $ID, $quotedPathAST, buildOption=$buildOption)))
        else
            (model, instantiatedSubmodelStruct) = Core.eval(modelModule, :($buildFunction($model, $modelModule, $FloatType, $TimeType, $instantiateModelOptions, $ID, $quotedPathAST)))
        end
        buildDict[ID] = instantiatedSubmodelStruct
        return model
    end

    for (key,value) in model
        if typeof(value) <: OrderedDict && haskey(value, :_class) && value[:_class] == :Model
            model[key] = buildSubmodels!(value, modelModule, FloatType, TimeType, instantiateModelOptions, buildDict; pathAST=appendSymbol(pathAST,key))
        end
    end
    return model
end



"""
    modelInstance = @instantiateModel(model; FloatType = Float64, aliasReduction=true, unitless=false,
        evaluateParameters=false, saveCodeOnFile="", log=false, logModel=false, logDetails=false, logStateSelection=false,
        logCode=false,logExecution=logExecution, logCalculations=logCalculations, logTiming=false, logFile=true)

Instantiates a model, i.e. performs structural and symbolic transformations and generates a function for calculation of derivatives suitable for simulation.

* `model`: model (declarations and equations)
* `FloatType`: Variable type for floating point numbers, for example: Float64, Measurement{Float64}, StaticParticles{Float64,100}, Particles{Float64,2000}
* `aliasReduction`: Perform alias elimination and remove singularities
* `unitless`: Remove units (useful while debugging models and needed for MonteCarloMeasurements)
* `evaluateParameters`: Use evaluated parameters in the generated code.
* `saveCodeOnFile`: If non-empty string, save generated code in file with name `saveCodeOnFile`.
* `log`: Log the different phases of translation
* `logModel`: Log the variables and equations of the model
* `logDetails`: Log internal data during the different phases of translation
* `logStateSelection`: Log details during state selection
* `logCode`: Log the generated code
* `logExecution`: Log the execution of the generated code (useful for timing compilation)
* `logCalculations`: Log the calculations of the generated code (useful for finding unit bugs)
* `logTiming`: Log timing of different phases
* `logFile`: Log file and line number where @instantiatedModel is called
* `return modelInstance prepared for simulation`
"""
macro instantiateModel(model, kwargs...)
    modelName = string(model)
    source = string(__source__.file)*":"*string(__source__.line)

    code = :( instantiateModel($model; modelName=$modelName, modelModule=@__MODULE__, source=$source, $(kwargs...) ) )
    return esc(code)
end


"""
See documentation of macro [`@instantiateModel`]
"""
function instantiateModel(model; modelName="", modelModule=nothing, source=nothing, FloatType = Float64, aliasReduction=true, unitless=false,
    log=false, logModel=false, logDetails=false, logStateSelection=false, logCode=false,
    logExecution=logExecution, logCalculations=logCalculations, logTiming=false, logFile=true, evaluateParameters=false, saveCodeOnFile="")
    if isMonteCarloMeasurements(FloatType) && !unitless
        unitless=true
        printstyled("  @instantiateModel(...,unitless=true, ..) set automatically, because\n  FloatType=MonteCarloMeasurements often fails if units are involved.\n", color=:red)
    end

    #try
    #    model = JSONModel.cloneModel(model, expressionsAsStrings=false)
        if logFile
            println("\nInstantiating model $modelName\n  in module: $modelModule\n  in file: $source")
        end
        resetEventCounters()
        global to = TimerOutput()

        modelStructure = ModelStructure()

        if isexpr(model, :quote)
            model = eval(model) # model defined with macro
        end

        if typeof(model) <: NamedTuple || typeof(model) <: Dict || typeof(model) <: OrderedDict
            # Traverse model and execute functions _buildFunction(..), to include code into sub-models
            buildDict = OrderedDict{String,Any}()
            instantiateModelOptions = OrderedDict{Symbol, Any}(
                :modelName          => modelName,
                :source             => source,
                :aliasReduction     => aliasReduction,
                :unitless           => unitless,
                :log                => log,
                :logModel           => logModel,
                :logDetails         => logDetails,
                :logStateSelection  => logStateSelection,
                :logCode            => logCode,
                :logExecution       => logExecution,
                :logCalculations    => logCalculations,
                :logTiming          => logTiming,
                :logFile            => logFile,
                :evaluateParameters => evaluateParameters,
                :saveCodeOnFile     => saveCodeOnFile)

            TimeType = if FloatType <: Measurements.Measurement ||
                          FloatType <: MonteCarloMeasurements.AbstractParticles;
                          baseType(FloatType) else FloatType end  # baseType(..) is defined in CodeGeneration.jl
            model = deepcopy(model)
            model = buildSubmodels!(model, modelModule, FloatType, TimeType, instantiateModelOptions, buildDict)

            if logModel
                @showModel(model)
            end
            if false
                println("drawModel")
                jsonDiagram = drawModel(modelName, model)
                println()
#                println("jsonModel")
#                JSON.print(jsonDiagram, 2)
            end

            if logTiming
                println("Flatten")
                @timeit to "flatten" flattenModelTuple!(model, modelStructure, modelName, to; unitless, log)
            else
                flattenModelTuple!(model, modelStructure, modelName, to; unitless, log)
            end

            if experimentalTranslation
                interface = buildInterface(model, modelStructure)
            end

            pars = OrderedDict{Symbol, Any}([(Symbol(k),v) for (k,v) in modelStructure.parameters])
            if length(modelStructure.equations) > 0
                flatModel = Model() | pars | Model(equations = [removeBlock(e) for e in modelStructure.equations])
            else
                flatModel = Model() | pars
            end
            printModelStructure(modelStructure, log=false)
            name = modelName
        else
            @show model
            error("Invalid model format")
        end

        duplicateMultiReturningEquations!(modelStructure.equations)

        allVariables = Incidence[]
        if logTiming
            println("Find incidence")
            @timeit to "findIncidence!" findIncidence!(modelStructure.equations, allVariables)
        else
            findIncidence!(modelStructure.equations, allVariables)
        end
        unique!(allVariables)

        if ! experimentalTranslation
            unknowns = setdiff(allVariables, keys(modelStructure.parameters), keys(modelStructure.inputs), [:time, :instantiatedModel, :_leq_mode, :_x, :_path])
        else
            unknowns = setdiff(allVariables, keys(modelStructure.parameters), [:time, :instantiatedModel, :_leq_mode, :_x, :_path])
        end
        Avar, states, derivatives = setAvar(unknowns)
        vActive = [a == 0 for a in Avar]

        if logStatistics || log
            println("Number of states:    ", length(states))
            if length(unknowns)-length(states) != length(modelStructure.equations)
                println("Number of unknowns:  ", length(unknowns)-length(states))
            end
            println("Number of equations: ", length(modelStructure.equations))
        end
        printArray(modelStructure.parameters, "Parameters:", log=log)
        printArray(modelStructure.inputs, "Inputs:", log=log)
        printArray(modelStructure.outputs, "Outputs:", log=log)
        printArray(states, "Potential states:", log=log)
        printArray(setdiff(unknowns, states), "Unknowns:", log=log)
        printArray(modelStructure.equations, "Equations:", log=log)
        if logDetails
            @show modelStructure.parameters modelStructure.mappedParameters modelStructure.init modelStructure.start modelStructure.variables modelStructure.flows modelStructure.inputs modelStructure.outputs
        end

        unknownsWithEliminated = unknowns
        if aliasReduction
            if log || logTiming
                println("Perform alias elimination and remove singularities")
            end
#            @show unknowns
#            @show keys(modelStructure.inputs) keys(modelStructure.outputs)
###            nonInputAndOutputs = setdiff(unknowns, keys(modelStructure.inputs), keys(modelStructure.outputs))
            if logTiming
                @timeit to "performAliasReduction" equations, unknowns, Avar, G, states, vEliminated, vProperty = performAliasReduction(unknowns, modelStructure.equations, Avar, logDetails, log)
                println("Number of reduced unknowns:  ", length(unknowns)-length(states))
                println("Number of reduced equations: ", length(equations))
            else
                equations, unknowns, Avar, G, states, vEliminated, vProperty = performAliasReduction(unknowns, modelStructure.equations, Avar, logDetails, log)
            end
            printArray(states, "States:", log=log)
            printArray(setdiff(unknowns, states), "Unknowns after alias reduction:", log=log)
            printArray(equations, "Equations after alias reduction:", log=log)
        else
            equations = modelStructure.equations
            G = Vector{Int}[] # Array{Array{Int64,1},1}()
            variablesIndices = OrderedDict(key => k for (k, key) in enumerate(unknowns))

            println("Build incidence matrix")
            unknownsSet = Set(unknowns)
            @timeit to "build incidence matrix" for i in 1:length(equations)
                e = equations[i]
                nameIncidence, coefficients, rest, linear = getCoefficients(e)

                incidence = [] # [variablesIndices[n] for n in nameIncidence if n in unknownsSet]
                for n in nameIncidence
                    if n in keys(variablesIndices)
                        push!(incidence, variablesIndices[n])
                    end
                end
                unique!(incidence)
                push!(G, incidence)
            end
            vEliminated = Int[]
            vProperty   = Int[]
        end

        modStructure = assignAndBLT(model, modelName, modelModule, equations, unknowns, modelStructure, Avar, G, states, logDetails, log, logTiming)
        (unknowns, equations, G, Avar, Bequ, assign, blt, parameters) = modStructure

        Gexplicit = Vector{Int}[]
        variablesIndices = OrderedDict(key => k for (k, key) in enumerate(unknowns))
        @timeit to "build explicit incidence matrix" for i in 1:length(equations)
            e = equations[i]
            if isexpr(e.args[2], :call) && e.args[2].args[1] in [:implicitDependency, :_DUPLICATIMPLICITDEPENDENCY]
                if e.args[2].args[1] == :_DUPLICATIMPLICITDEPENDENCY
                    e.args[2] = e.args[2].args[2]
                    e.args[2].args[1] = :_DUPLICATEEQUATION
                else
                    e.args[2] = e.args[2].args[2]
                end

                nameIncidence = Incidence[]
                findIncidence!(e, nameIncidence)

                incidence = [] # [variablesIndices[n] for n in nameIncidence if n in unknownsSet]
                for n in nameIncidence
                    if n in keys(variablesIndices)
                        push!(incidence, variablesIndices[n])
                    end
                end
                unique!(incidence)
                push!(Gexplicit, incidence)
                equations[i] = e
            else
                push!(Gexplicit, G[i])
            end
        end

        if false
            println("Explicit equations:")
            for e in equations
                println(e)
            end

            @show G Gexplicit
        end

        if ! experimentalTranslation
            inst = stateSelectionAndCodeGeneration(modStructure, Gexplicit, name, modelModule, buildDict, FloatType, TimeType, modelStructure.init, modelStructure.start, modelStructure.inputs, modelStructure.outputs,
                vEliminated, vProperty, unknownsWithEliminated, modelStructure.mappedParameters, modelStructure.hideResults;
                unitless, logStateSelection, logCode, logExecution, logCalculations, logTiming, evaluateParameters, saveCodeOnFile)
        else
            interface[:equations] = modStructure[:equations]
            return interface
        end

        if logTiming
            show(to, compact=true)
            println()
        end

        inst #, flatModel
#=
    catch e
        if isa(e, ErrorException)
            println()
            printstyled("Model error: ", bold=true, color=:red)
            printstyled(e.msg, "\n", bold=true, color=:red)
            printstyled("Aborting @instantiateModel($modelName,...) in $modelModule.\n", bold=true, color=:red)
            println()
#            Base.rethrow()
        else
            Base.rethrow()
        end
        if logTiming
            show(to, compact=true)
            println()
        end
    end
=#

end
