"""
Modia module with block component models (inspired from Modelica Standard Library.

* Developer: Hilding Elmqvist, Mogram AB  
* Copyright (c) 2016-2018: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)

"""
module Blocks

using ..Instantiation
using Modia
export FirstOrder, Feedback, PI, Step, Sine, Switch, MIMO



@model Block begin
end

@model SISO begin
    # Single Input Single Output continuous control block
  @extends Block()
  u=Float()
  y=Float()
end

@model FirstOrder begin
    # First order transfer function block (= 1 pole)
  k=1   # Gain
  T=1   # Time Constant

  @extends SISO()
  @inherits u, y
@equations begin 
  der(y) = (k*u - y)/T
  end
end 

@model Feedback begin
    # Output difference between commanded and feedback input
  u1=Float()
  u2=Float()
  y=Float()
@equations begin 
  y = u1 - u2
  end
end 

@model PI begin
    # Proportional-Integral controller  
  k=1   # Gain
  T=Parameter(1, min=1E-10)   # Time Constant (T>0 required)

  @extends SISO()
  @inherits u, y
  x=Float(start=0)   # State of block 
@equations begin 
  der(x) = u/T
  y = k*(x + u)
  end
end

@model SO begin
    # Single Output continuous control block
  @extends Block()
  y=Variable()
end

@model SignalSource begin
    # Base class for continuous signal source
  @extends SO()
  offset=0      # Offset of output signal y
  startTime=0   # Output y = offset for time < startTime
end

@model Step begin
    # Generate step signal of type Real
  height=1      # Height of step
  @extends SignalSource()
  @inherits y, offset, startTime
  t = Float(start=0.0)
@equations begin 
  y = offset + (t < startTime ? 0 : height)
  der(t) = 1
  end
end

@model Sine begin
    # Generate sine signal
  amplitude=1 # Amplitude of sine wave
  freqHz=Parameter() # Frequency of sine wave
  phase=0 # Phase of sine wave
  offset=0      # Offset of output signal y
  startTime=0   # Output y = offset for time < startTime
  @extends SO()
  @inherits y
  t = Float(start=0.0)
@equations begin 
#  y = offset + if time < startTime;  0 else amplitude*sin(2*pi*freqHz*(time - startTime) + phase) end
  y = offset + if t < startTime;  0 else amplitude*sin(2*pi*freqHz*(t - startTime) + phase) end
  der(t) = 1
  end
end

@model Sine2 begin
    # Generate sine signal
  amplitude=1 # Amplitude of sine wave
  freqHz=Parameter() # Frequency of sine wave
  phase=0 # Phase of sine wave
  @extends SignalSource()
  @inherits y, offset, startTime
@equations begin 
  y = offset + if time < startTime;  0 else amplitude*sin(2*pi*freqHz*(time - startTime) + phase) end
  end
end

@model Switch begin
  sw=Boolean()
  u1=Variable(); u2=Variable()
  y=Variable()
@equations begin
  y = if sw; u1 else u2 end
  end
end

@model ABCD begin
  A=-1; B=1; C=1; D=0
  u=Float(); y=Float()
  x=Float(start=0)
@equations begin
  der(x) = A*x + B*u
  y = C*x + D*u
  end
end


end