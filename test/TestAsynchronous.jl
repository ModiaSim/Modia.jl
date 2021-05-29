module TestAsynchronous

#println("\nAsynchronousExamples: Demonstrating the ability to simulate models with asynchronous semantics")

using TinyModia, ModiaPlot
using DifferentialEquations

# References:
# http://osp.mans.edu.eg/cs212/Seq_circuits_Summary_FF-types.htm
# https://www.kth.se/social/files/56124acff27654028dfc0f80/F12asyFSM1_eng.pdf

BooleanPulse = Model(
	startTime = parameter,
	width = parameter,
	equations = :[
 		y = positive(ustrip(time-startTime)) && ! positive(ustrip(time-(startTime+width)))
		#y = time > startTime && ! (time > startTime+width)
	]
)

TestBooleanPulse = Model(
	set = BooleanPulse | Map(startTime=2u"s", width=2u"s"),
)

model = @instantiateModel(TestBooleanPulse, log=true, logCode=true)
simulate!(model, Tsit5(), stopTime=10, log=true, logEvents=true)
plot(model, ("set.y",), heading="TestBooleanPulse", figure=1)


SRFlipFlop = Model(
    Q = Var(init=false),
    equations = :[
		Q = S || ! R && previous(Q)
    ]
) 

TestSRFlipFlop = Model(
	set = BooleanPulse | Map(startTime=2u"s", width=2u"s"),
	reset = BooleanPulse | Map(startTime=5u"s", width=2u"s"),
	sr = SRFlipFlop,
	equations = :[
		connect(set.y, sr.S)
		connect(reset.y, sr.R)
	]
)

model = @instantiateModel(TestSRFlipFlop, log=true, logCode=true)
simulate!(model, Tsit5(), stopTime=10, log=true, logEvents=true)
plot(model, [("sr.S", "sr.R"), ("sr.Q")], heading="TestSRFlipFlop", figure=1)

end