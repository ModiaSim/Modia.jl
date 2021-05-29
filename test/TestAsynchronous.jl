module TestAsynchronous

#println("\nAsynchronousExamples: Demonstrating the ability to simulate models with asynchronous semantics")

using TinyModia, ModiaPlot
using DifferentialEquations

BooleanPulse1 = Model(
	startTime = parameter,
	width = parameter,
	equations = :[
        y = after(startTime) && ! after(startTime + width)
	]
)

TestBooleanPulse1 = Model(
	pulse = BooleanPulse1 | Map(startTime=2u"s", width=2u"s"),
)

model = @instantiateModel(TestBooleanPulse1, log=true, logCode=true)
simulate!(model, Tsit5(), stopTime=10, log=true, logEvents=true)
plot(model, [("pulse.y"), ("pulse.startTime")], heading="TestBooleanPulse1", figure=1)


BooleanPulse = Model(
	width = parameter | Map(value=50, min=0, max=100) | info"Width of pulse in % of period",
	period = parameter | Map(min=0) | info"Time for one period",
	startTime = parameter | info"Time instant of first pulse",
	pulseStart = Var() | info"Start time of pulse",
	y = output,
	equations = :[
		Twidth = period*width/100
#		clock1 = Clock(startTime, period)  # Clock with 4 arguments not available yet.
		clock1 = Clock(period)  # Gives unit error if unitless=false
		pulseStart = sample(time, clock1)
        y = after(pulseStart) && ! after(pulseStart + Twidth)
	]
)

TestBooleanPulse = Model(
	pulse = BooleanPulse | Map(startTime=2u"s", period=2u"s"),
)

model = @instantiateModel(TestBooleanPulse, log=true, logCode=true)
simulate!(model, Tsit5(), stopTime=10, log=true, logEvents=true)
plot(model, [("pulse.y"), ("pulse.pulseStart")], heading="TestBooleanPulse", figure=1)


# References:
# http://osp.mans.edu.eg/cs212/Seq_circuits_Summary_FF-types.htm
# https://www.kth.se/social/files/56124acff27654028dfc0f80/F12asyFSM1_eng.pdf

SRFlipFlop = Model(
    Q = Var(init=false),
    equations = :[
		Q = S || ! R && previous(Q)  # Currently gives: ERROR: LoadError: MethodError: no method matching pre(::Bool, ::TinyModia.SimulationModel{Float64, Float64}, ::Int64)
    ]
) 

TestSRFlipFlop = Model(
	set = BooleanPulse1 | Map(startTime=2u"s", width=2u"s"),
	reset = BooleanPulse1 | Map(startTime=5u"s", width=2u"s"),
	sr = SRFlipFlop,
	equations = :[
		connect(set.y, sr.S)
		connect(reset.y, sr.R)
	]
)

model = @instantiateModel(TestSRFlipFlop, log=true, logCode=true, unitless=true)
simulate!(model, Tsit5(), stopTime=10, log=true, logEvents=true)
plot(model, [("sr.S", "sr.R"), ("sr.Q")], heading="TestSRFlipFlop", figure=1)

end