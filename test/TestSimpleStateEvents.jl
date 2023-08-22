module TestSimpleStateEvents

using Modia
@usingModiaPlot


SimpleStateEvents = Model(
    fmax = 1.5,
    m    = 1.0,
    k    = 1.0,
    d    = 0.1,
    s    = Var(init = 2.0),
    v    = Var(init = 0.0),
    equations = :[
        sPos = positive(s)
        f = if sPos; 0.0 else fmax end   # or: f = sPos ? 0.0 : fmax
        v = der(s)
        m*der(v) + d*v + k*s = f
    ]
)

model = @instantiateModel(SimpleStateEvents)

simulate!(model, Tsit5(), stopTime = 10, log=false, logEvents=false, requiredFinalStates = [1.1756992201144405, -0.5351307571761175])

plot(model, ["s", "v", "sPos", "f"])

# Test IDA
simulate!(model, IDA(), nlinearMinForDAE=1, stopTime = 10, log=false, logEvents=false, requiredFinalStates = [1.1756992201144405, -0.5351307571761175])

end