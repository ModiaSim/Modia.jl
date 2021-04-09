module BouncingBallSimulation

using TinyModia
using DifferentialEquations
using ModiaPlot

function reInit(v,e)
    v = -e*v
    if v < 0.01
        v = 0.0
    end
    return v
end

BouncingBall = Model(
    e = 0.7, 
    g = 9.81, 
    h = Var(init = 1.0),
    v = Var(init = 0.0),
    flying = Var(start = true),
    equations = :[
        contact = edge(instantiatedModel, 1, -h, "-h", _leq_mode),
        if contact
            v = reInit(v,e)
            flying = v > 0.0      
        end,
        der(h) = v,
        der(v) = flying ? -g : 0.0
    ]
)

model = @instantiateModel(BouncingBall, logCode=true)

simulate!(model, Tsit5(), stopTime = 3.0, nz=1, log=true, logEvents=true)   # requiredFinalStates = [-0.3617373025974107]

plot(model, ["h", "v", "flying"])

end