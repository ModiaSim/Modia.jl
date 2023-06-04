module TestBouncingBall

using Modia
@usingModiaPlot

# Does not work: edge(..) appears in a linear equation system

BouncingBall = Model(
    e = 0.7,
    g = 9.81,
    h = Var(init = 1.0),
    v = Var(init = 0.0),
    flying = Var(start = true),   # due to pre(flying) -> InstantiatedModel(....; pre_startValues = [true], ....)
    equations = :[
        # desired: 
        #   flying = edge(-h) ? reinit(v, -e*v) : pre(flying)
        flying = edge(instantiatedModel, 1, -h, "-h", _leq_mode) ?
                    !reinit(instantiatedModel, _x, 2, -e*v, _leq_mode) : pre(flying),
        der(h) = v,
        der(v) = flying ? -g : 0.0
    ]
)


model = @instantiateModel(BouncingBall, logCode=true)
@show model.equationInfo.x_info

# Temporary code, until edge, reinit and pre are supported in the Modia language
    eh = model.eventHandler
    eh.nz = 1
    eh.z  = ones(eh.nz)
    eh.zPositive = fill(false, eh.nz)
    model.pre = [true]

simulate!(model, Tsit5(), stopTime = 3.0, log=true, logEvents=true)   # requiredFinalStates = [-0.3617373025974107]

plot(model, ["h", "v", "flying"])

end