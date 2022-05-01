module TestLinearSystems

using Modia
using Modia.OrderedCollections
using Modia.LinearAlgebra
@usingModiaPlot

"""
    ls = LinearStateSpace(; A, B, C, x_init=nothing)

Define a linear state space system Model:

```math
\begin{aligned}
\frac{dx}{dt} &= A*x + B*u \\
            y &= C*x
\end{aligned}
```

where

- `x_init` is the optional vector of init-values. If not provided, zeros will be used.
- The number of **inputs** and **outputs** is **fixed**, after @instantiateModel(..) was called.
- The number of **states** can be **changed** before simulation starts,
  by providing appropriate ``A,B,C`` matrices and `x_init` vector as `merge` value in `simulate!(..., merge = ...)`.

# Example
```
using Modia
usingModiaPlot
SSTest = Model(ss = ModelLinearStateSpace(A=[-1/0.1;;], B=[2.0/0.1;;], C=[1.0;;])), # one State
               y = Var(start=0.0)
               equations = :[ss.u[1] = 1.0,
                             y = ss.y[1]]
         )
ssTest = @instantiateModel(SSTest, merge=Map(ss = Map(A=[-1/0.1 0.0;
                                                          0.0  -1/0.2], B=[2.0/0.1;2.0/0.2;;], C=[1.0 2.0] ) # two states
simulate!(ssTest, stopTime=5.0)
plot(ssTest, ("ss.x", "u", "y"))
```
"""
ModelLinearStateSpace(; kwargs...) = Model(; _buildFunction = :(buildLinearStateSpace!),         # Called once in @instantiateModel(..) before getDerivatives!(..) is generated
                                             _stateInfoFunction = Par(functionName = :(stateInfoLinearStateSpace!)),  # Called once after new A,B,C values are merged
                                             kwargs...)

mutable struct LinearStateSpaceStruct{FloatType}
    path::String  # Path name of instance
    ix::Int       # Index with respect to equationInfo.x_info    
    A::Matrix{FloatType}
    B::Matrix{FloatType}
    C::Matrix{FloatType}
    x_init::Vector{FloatType}  # Initial values of states
    y::Vector{FloatType}       # Internal memory for y
    x::Vector{FloatType}       # Internal memory for x
    derx::Vector{FloatType}    # Internal memory for derx

    function LinearStateSpaceStruct{FloatType}(; A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
                                                 x_init::Union{AbstractVector,Nothing}=nothing, 
                                                 nu::Int, ny::Int,
                                                 u::AbstractVector, y::AbstractVector,  # Code generated with buildLinearStateSpace! provides start values of u and y.
                                                 path::String, kwargs...) where {FloatType}
        #println("... 4: LinearStateSpaceStruct called for $path")
        if length(kwargs) > 0
            @warn "LinearStateSpaceStruct with path=$path ignores keyword arguments: $(kwargs...)"
        end
        @assert(size(A,2) == size(A,1))
        @assert(size(B,1) == size(A,1))
        @assert(size(C,2) == size(A,1))
        @assert(size(u,1) == size(B,2))
        @assert(size(y,1) == size(C,1))
        @assert(size(u,1) == nu)
        @assert(size(y,1) == ny)
        copyA = Matrix{FloatType}(deepcopy(A))
        copyB = Matrix{FloatType}(deepcopy(B))
        copyC = Matrix{FloatType}(deepcopy(C))
        copy_x_init = if isnothing(x_init); zeros(FloatType, size(A,1)) else Vector{FloatType}(deepcopy(x_init)) end
        new(path, 0, copyA, copyB, copyC, copy_x_init, zeros(FloatType,size(C,1)), zeros(FloatType, size(A,1)), zeros(FloatType,size(A,1)))
    end
end

mutable struct LinearStateSpaceBuild{FloatType}
    # Cannot be changed after @instantiatedModel(..)
    path::String # Path name of instance, e.g. "a.b.c"
    nu::Int      # Number of inputs
    ny::Int      # Number of outputs

    # Can be changed after @instantiateModel
    ls::Union{LinearStateSpaceStruct{FloatType}, Nothing}

    function LinearStateSpaceBuild{FloatType}(path::String, nu::Int, ny::Int) where {FloatType}
        #println("... 2: LinearStateSpaceBuild called with path = $path")
        new(path,nu,ny,nothing)
    end
end


function buildLinearStateSpace!(model::AbstractDict, FloatType::Type, TimeType::Type,
                                buildDict::OrderedCollections.OrderedDict{String,Any},
                                path::Union{Expr,Symbol,Nothing})
    # Called from @instantiatedModel, during instantiation of the model.
    pathAsString = isnothing(path) ? "" : string(path)
    #println("... 1: buildLinearStateSpace! called for path = ", pathAsString)

    # Determine nu,ny from model (must be Integer literals >= 0)
    nu::Int = model[:nu]
    ny::Int = model[:ny]
    @assert(nu >= 0)
    @assert(ny >= 0)
    u_zeros = zeros(FloatType,nu)
    y_zeros = zeros(FloatType,ny)

    # Define code to be generated
    lsCode = Model(ls      = Var(hideResult=true),
                   success = Var(hideResult=true),
                         u = Var(input  = true, start = u_zeros),
                         y = Var(output = true, start = y_zeros),
                    equations = :[
                        ls = getLinearStateSpace!(instantiatedModel, $pathAsString)
                        y = computeOutputs!(instantiatedModel, ls)
                        success = computeStateDerivatives!(instantiatedModel, ls, u)])

    # Store build info in buildDict
    buildDict[pathAsString] = LinearStateSpaceBuild{FloatType}(pathAsString, nu, ny)
    return lsCode
end


function stateInfoLinearStateSpace!(model::AbstractDict, FloatType::Type, TimeType::Type,
                                    buildDict::OrderedCollections.OrderedDict{String,Any},
                                    eqInfo::Modia.EquationInfo,
                                    path::String)::Nothing
    # Called during evaluation of the parameters (before initialization)                                 
    #println("... 3: stateInfoLinearStateSpace! called for $path with model = $model")
    lsBuild::LinearStateSpaceBuild{FloatType} = buildDict[path]
    ls = LinearStateSpaceStruct{FloatType}(; path, model...)
    A = ls.A
    @assert(size(A,2) == size(A,1))
    @assert(size(ls.B,2) == lsBuild.nu)
    @assert(size(ls.C,1) == lsBuild.ny)
    ls.ix = Modia.addOrUpdateStateInfo(eqInfo, path*".x", path*".der(x)", ls.x_init)
    lsBuild.ls = ls    
    return nothing
end
    

function getLinearStateSpace!(instantiatedModel::SimulationModel{FloatType,TimeType}, path::String)::LinearStateSpaceStruct{FloatType} where {FloatType,TimeType}
    ls = instantiatedModel.buildDict[path].ls
    copyState!(instantiatedModel, ls.ix, ls.x)
    return ls
end

function computeOutputs!(instantiatedModel, ls)
    mul!(ls.y, ls.C, ls.x)
    return ls.y
end

function computeStateDerivatives!(instantiatedModel, ls, u)::Bool
    # ls.derx .= ls.A*ls.x + ls.B*u
    mul!(ls.derx, ls.A, ls.x)
    mul!(ls.derx, ls.B, u, 1.0, 1.0)
    Modia.set_hiddenStateDerivative!(instantiatedModel, ls.ix, ls.derx)
    return true
end
       
# T*der(x) + x = u 
T = 0.2;
SSTest = Model(
            ss = ModelLinearStateSpace(A=[-1.0/T;;], B=[1.0/T;;], C=[0.9;;], x_init=[0.2], nu=1, ny=1),
            equations = :[ss.u = [2.0],
                          y = ss.y]
         )
         
ssTest = @instantiateModel(SSTest, logCode=true)
simulate!(ssTest, stopTime=1.0, log=true, logStates=true)
#Modia.printResultInfo(ssTest)
plot(ssTest, ("ss.x", "ss.u", "ss.y"))

end