module TestSynchronous

#println("\nSynchronousExamples: Demonstrating the ability to simulate models with synchronous semantics")

using Modia
@usingModiaPlot


MassWithSpringDamper = Model(
    m = 1.0,
    k = 1.0,           # spring constant
    d = 0.1,           # damping coefficient
    x = Var(init=0.0), # Position
    v = Var(init=5.0), # Velocity
    equations = :[
        der(x) = v
        m*der(v) = f - k*x - d*v
    ]
) 

SpeedControl = MassWithSpringDamper | Model(
    k    = 0.0,
    K    = 5.0,   # Gain of speed P controller
    vref = 100.0, # Speed ref.
    vd = Var(start=0.0),
    u  = Var(start=0.0),  
    equations = :[ 
        vd = sample(v, Clock(0.1))    # speed sensor
        u  = K*(vref-vd)              # P controller for speed
        f  = hold(u, Clock(0.1))      # force actuator
    ]
)

speedControl = @instantiateModel(SpeedControl)
simulate!(speedControl, Tsit5(), stopTime=1.5, logEvents=false, requiredFinalStates = [133.3845202465569, 98.03694943140056] )
plot(speedControl, [("v", "vd"), ("u","f")], heading="SpeedControl", figure=1)


SpeedControlPI = MassWithSpringDamper | Model(
    k    = 0.0,
    K    = 5.0,   # Gain of speed P controller
	Ti = 2, 	  # Integral time 
    vref = 100.0, # Speed ref.
	dt = 0.1,     # sampling interval
    vd = Var(start=0.0),
    u  = Var(start=0.0),  
	e = Var(start=0.0),
	intE = Var(start=10.0),  # To check that previous(intE) is corrected initialized
    equations = :[ 
		vd = sample(v, Clock(0.1))     # speed sensor
		
		# PI controller for speed
		e = vref-vd
        previous_intE = previous(intE, Clock(0.1))
		intE = previous_intE + e
		u = K*(e + intE/Ti)

		# force actuator
		f = hold(u, Clock(0.1))
    ]
)

speedControlPI = @instantiateModel(SpeedControlPI)
simulate!(speedControlPI, Tsit5(), stopTime=1.5, logEvents=false, requiredFinalStates = [145.828878491462, 99.95723678839984])
plot(speedControlPI, [("v", "vd"), ("u","f"), ("previous_intE", "intE")], heading="SpeedControlPI", figure=2)


end