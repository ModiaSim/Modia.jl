"""
Modia module with block component models (inspired from Modelica Standard Library).

* Developer: Hilding Elmqvist, Mogram AB, Martin Otter, DLR  
* Copyright (c) 2016-2021: Hilding Elmqvist, Martin Otter
* License: MIT (expat)

"""
#module Blocks

using Modia

#export Gain, FirstOrder, Feedback, PI, Step, Ramp # Sine, Switch, MIMO


# Single-Input-Single-Output continuous control block
SISO = Model(
    u = input, 
    y = output
)

# Gain
Gain = SISO | Model(
    k = 1, # (info = "Gain")
    equations = :[
        y = k*u ]
)

# First-order transfer function block (= 1 pole)
FirstOrder = SISO | Model(
    k = 1.0,
    T = 1.0*u"s",
    x = Var(init=0.0),
    equations = :[
        der(x) = (k * u - x) / T
        y = x ]
)

# Output difference between commanded and feedback input
Feedback = Model(
    u1 = input | info"Input 1",
    u2 = input | info"Input 2",
    y = output | info"Output signal",
    equations = :[
        y = u1 - u2 ]
) 

# Proportional-Integral controller
PI = SISO | Model(
    k = 1.0, # (info = "Gain")
    T = 1.0u"s", # (min = 1E-10, info = "Time Constant (T>0 required)")
    x = Var(init=0.0),
    equations = :[
        der(x) = u / T
        y = k * (x + u) ]
)

# Single-Output continuous control block
SO = Model(
    y = output 
)

# Base class for a continuous signal source
SignalSource = SO | Model(
    offset = 0.0, # info = "Offset of output signal y"
    startTime = 0.0*u"s" # info = "Output y = offset for time < startTime")
)

# Step signal
Step = SignalSource | Model(
    height = 1.0,
    equations = :[
        y = offset + (time < startTime ? 0*height : height) ]  # 0*height is needed, if height has a unit
)

# Ramp signal
Ramp = SignalSource | Model(
    height   = 1.0,
    duration = 2.0u"s",
    equations = :[
        y = offset + (time < startTime ? 0.0*height :   # 0*height is needed, if height has a unit
                     (time < startTime + duration ? (time - startTime)*height/duration :
                                                    height)) ]
)


# Linear state space system
StateSpace = Model(   
    A = parameter | fill(0.0,0,0),
    B = parameter | fill(0.0,0,0),   
    C = parameter | fill(0.0,0,0), 
    D = parameter | fill(0.0,0,0),
    u = input, 
    y = output,
    x = Var(init = zeros(0)),
    equations = :[
        der(x) = A*x + B*u
             y = C*x + D*u
    ]
)




# -------------------------------------------------------

#=
# Sinusoidal signal
@model Sine begin
    amplitude = Parameter(1.0, info = "Amplitude of sine wave")
    freqHz = Parameter(1.0,info = "Frequency of sine wave")
    phase = Parameter(0.0, info = "Phase of sine wave")
    offset = Parameter(0.0, info = "Offset of output signal y")
    startTime = Parameter(0.0, info = "Output y = offset for time < startTime")
    SO()
equations
    y = offset + if time < startTime;  0 else amplitude*sin(2*pi*freqHz*(time - startTime) + phase) end
end

@model Sine2 begin
    # Generate sine signal
    amplitude = Parameter(1.0, info = "Amplitude of sine wave")
    freqHz = Parameter(1.0,info = "Frequency of sine wave")
    phase = Parameter(0.0, info = "Phase of sine wave")
    SignalSource()
equations
    y = offset + if time < startTime;  0 else amplitude * sin(2 * pi * freqHz * (time - startTime) + phase) end
end

# Switch
@model Switch begin
    sw = Input(Boolean(info = "Switch position (if `true`, use `u1`, else use `u2`)"))
    u1 = Input(info = "Input 1")
    u2 = Input(info = "Input 2")
    y = Output(info = "Output signal")
equations
    y = if sw; u1 else u2 end
end

# ABCD model
@model ABCD(A = -1, B = 1, C = 1, D = 0) begin
    u = Input(info = "Input signal"); y = Output(info = "Output signal")
    x = Float(start = 0)
equations
    der(x) = A * x + B * u
    y = C * x + D * u
end
=#
#end
