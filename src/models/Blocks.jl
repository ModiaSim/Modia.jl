"""
Modia module with block component models (inspired from Modelica Standard Library).

* Developer: Hilding Elmqvist, Mogram AB  
* Copyright (c) 2016-2018: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)

"""
module Blocks

#using ..Instantiation
using Modia
export FirstOrder, Feedback, PI, Step, Sine, Switch, MIMO



@model Block begin
end

"""
Single-Input-Single-Output continuous control block
"""
@model SISO begin
    u = Float(info = "Input signal")
    y = Float(info = "Output signal")
    @extends Block()
end

"""
First-order transfer function block (= 1 pole)
"""
@model FirstOrder begin
    k = Parameter(1.0, info = "Gain")
    T = Parameter(1.0, info = "Time Constant")

    @extends SISO()
    @inherits u, y
    @equations begin 
        der(y) = (k * u - y) / T
    end
end 

"""
Output difference between commanded and feedback input
"""
@model Feedback begin
    u1 = Float(info = "Input 1")
    u2 = Float(info = "Input 2")
    y = Float(info = "Output signal")
    @equations begin 
        y = u1 - u2
    end
end 

"""
Proportional-Integral controller
"""
@model PI begin
    k = Parameter(1.0, info = "Gain")
    T = Parameter(1.0, min = 1E-10, info = "Time Constant (T>0 required)")

    @extends SISO()
    @inherits u, y
    x = Float(start=0)   # State of block 
    @equations begin 
        der(x) = u / T
        y = k * (x + u)
    end
end

"""
Single-Output continuous control block
"""
@model SO begin
    y = Variable(info = "Output signal")
    @extends Block()
end

"""
Base class for a continuous signal source
"""
@model SignalSource begin
    offset = Parameter(0.0, info = "Offset of output signal y")
    startTime = Parameter(0.0, info = "Output y = offset for time < startTime")
    @extends SO()
end

"""
Step signal
"""
@model Step begin
    height = Parameter(1.0, info = "Height of step")
    @extends SignalSource()
    @inherits y, offset, startTime
    t = Float(start=0.0)
    @equations begin 
        y = offset + (t < startTime ? 0 : height)
        der(t) = 1
    end
end

"""
Sinusoidal signal
"""
@model Sine begin
    amplitude = Parameter(1.0, info = "Amplitude of sine wave")
    freqHz = Parameter(info = "Frequency of sine wave")
    phase = Parameter(0.0, info = "Phase of sine wave")
    offset = Parameter(0.0, info = "Offset of output signal y")
    startTime = Parameter(0.0, info = "Output y = offset for time < startTime")
    @extends SO()
    @inherits y
    t = Float(start=0.0)
    @equations begin 
#  y = offset + if time < startTime;  0 else amplitude*sin(2*pi*freqHz*(time - startTime) + phase) end
        y = offset + if t < startTime;  0 else amplitude * sin(2 * pi * freqHz * (t - startTime) + phase) end
        der(t) = 1
    end
end

@model Sine2 begin
    # Generate sine signal
    amplitude = 1 # Amplitude of sine wave
    freqHz = Parameter() # Frequency of sine wave
    phase = 0 # Phase of sine wave
    @extends SignalSource()
    @inherits y, offset, startTime
    @equations begin 
        y = offset + if time < startTime;  0 else amplitude * sin(2 * pi * freqHz * (time - startTime) + phase) end
    end
end

"""
Switch
"""
@model Switch begin
    sw = Boolean(info = "Switch position (if `true`, use `u1`, else use `u2`)")
    u1 = Variable(info = "Input 1")
    u2 = Variable(info = "Input 2")
    y = Variable(info = "Output signal")
    @equations begin
        y = if sw; u1 else u2 end
    end
end

"""
ABCD model
"""
@model ABCD begin
    A = -1; B = 1; C = 1; D = 0
    u = Float(info = "Input signal"); y = Float(info = "Output signal")
    x = Float(start = 0)
    @equations begin
        der(x) = A * x + B * u
        y = C * x + D * u
    end
end


end
