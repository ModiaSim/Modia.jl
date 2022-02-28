"""
Modia module with electric component models (inspired from Modelica Standard Library).

* Developer: Hilding Elmqvist, Mogram AB  
* Copyright (c) 2016-2021: Hilding Elmqvist
* License: MIT (expat)

"""
#module Electric

using Modia

Var(args...; kwargs...) = (;args..., kwargs...)
Var(value::Union{Float64, Int64, Bool, String, Expr}, args...; kwargs...) = (;value = value, args..., kwargs...)

parameter = :parameter => true
input = :input => true
output = :output => true
potential = :potential => true
flow = :flow => true

v = Var(potential, nominal=10)
@show v

v = Var(5, parameter, min=0)
@show v

v = Var(potential, min=0, flow, nominal=10)
@show v

Pin = Model( v = Var(potential, nominal=10), i = Var(flow) )
@show Pin


Input(; kwargs...) = (;input=true, kwargs...)
Output(; kwargs...) = (;output=true, kwargs...)
Potential(; kwargs...) = (;potential=true, kwargs...)
Flow(; kwargs...) = (;flow=true, kwargs...)


Pin = Model( v = Potential(nominal=10), i = Flow() )

@show Pin



#Pin = Model( potentials = :[v], flows = :[i] )
Pin = Model( v = Var(;pot), i = Var(;flow) )

OnePort = Model( p = Pin, n = Pin, equations = :[
        v = p.v - n.v
        0 = p.i + n.i
        i = p.i ] 
    )

"""
  Resistor(R=1.0u"Ω")

Electrical resistor

`R` - Resistance Ω
"""  
Resistor = OnePort | Model( R = 1.0u"Ω", equations = :[ R*i = v ] )

Capacitor = OnePort | Model( C = 1.0u"F", v = Map(init=0.0u"V"), equations = :[ C*der(v) = i ] ) 

Inductor = OnePort | Model( L = 1.0u"H", init=Map(i=0.0u"A"), equations = :[ L*der(i) = v ] )

ConstantVoltage = OnePort | Model( V = 1.0u"V", equations = :[ v = V ] )

Ground = Model( p = Pin, equations = :[ p.v = 0.0u"V" ] ) 

# Ideal operational amplifier (norator-nullator pair), but 3 pins
IdealOpAmp3Pin = Model(
    in_p = Pin,
    in_n = Pin,
    out = Pin,
    equations = :[ 
        in_p.v = in_n.v
        in_p.i = 0u"A"
        in_n.i = 0u"A" ]
)

# Partial generic voltage source using the input signal as source voltage 
PartialSignalVoltage = Model(
    inputs = :[v],
    p = Pin,
    n = Pin
) 

# Generic voltage source using the input signal (without and with unit) as source voltage 
SignalVoltage = PartialSignalVoltage | Model(
    equations = :[
        p.v - n.v = v
                0 = p.i + n.i
                i = p.i ]
) 
UnitlessSignalVoltage = PartialSignalVoltage | Model(
    equations = :[
        p.v - n.v = v*u"V"
                0 = p.i + n.i
                i = p.i ]
) 


# Partial sensor to measure the current in a branch
PartialCurrentSensor = Model(
    outputs = :[i],
    p = Pin, # (info = "Positive pin")
    n = Pin, # (info = "Negative pin")
)

# Sensor to measure the current in a branch
CurrentSensor = PartialCurrentSensor | Model(
    equations = :[
        p.v = n.v
          0 = p.i + n.i    
          i = p.i]
)
UnitlessCurrentSensor = PartialCurrentSensor | Model(
    equations = :[
        p.v = n.v
        0 = p.i + n.i    
        i = p.i/u"A"]
)


#=
# Step voltage source
@model StepVoltage begin
    V = Parameter(1.0, start = 1.0, info = "Voltage") #, T = Unitful.V)
    startTime = Parameter(0.0, start = 0.0, info = "Start time") # , T = Unitful.s)
    OnePort()
equations
    v = if time < startTime; 0 else V end
end

@model VoltageSource begin
    OnePort()
    offset = Par(0.0) # Voltage offset
    startTime = Par(0.0) # Time offset 
    signalSource = SignalSource(offset=offset, startTime=startTime)
equations 
    v = signalSource.y
end 

#=
@model SineVoltage1 begin
    # Sine voltage source
    V = Parameter() # Amplitude of sine wave
    phase = Par(0.0) # Phase of sine wave
    freqHz = Parameter() # Frequency of sine wave
    VoltageSource(signalSource=Sine(amplitude=V, freqHz=freqHz, phase=phase))
end 
=#

# Sinusoidal voltage source
@model SineVoltage begin
    V = Parameter() # Amplitude of sine wave
    phase = Par(0.0) # Phase of sine wave
    freqHz = Parameter() # Frequency of sine wave
    VoltageSource()
equations
    v = V*sin(10*time)    
end 


# Ideal diode
@model IdealDiode begin # Ideal diode
    OnePort()
    Ron = Par(1.0E-2) # Forward state-on differential resistance (closed diode resistance)
    Goff = Par(1.0E-2) # Backward state-off conductance (opened diode conductance)
    Vknee = Par(0) # Forward threshold voltage
    # off = Variable(start=true) # Switching state
    s = Float(start=0.0) # Auxiliary variable for actual position on the ideal diode characteristic
    #=  
    s = 0: knee point
    s < 0: below knee point, diode conducting
    s > 0: above knee point, diode locking 
    =#
equations
    # off := s < 0
#        v = s * if !positive(s); 1 else Ron end + Vknee
#        i = s * if !positive(s); Goff else 1 end + Goff * Vknee
    v = s * if !(s>0); 1 else Ron end + Vknee
    i = s * if !(s>0); Goff else 1 end + Goff * Vknee
end 

@model Diode begin 
    OnePort()
    Ids=Par(1.e-6) # "Saturation current";
    Vt=Par(0.04) # "Voltage equivalent of temperature (kT/qn)";
    Maxexp = Par(15) # "Max. exponent for linear continuation";
    R=Par(1.e8) # "Parallel ohmic resistance";
equations
    i = if v/Vt > Maxexp; Ids*(exp(Maxexp)*(1 + v/Vt - Maxexp) - 1) else Ids*(exp(v/Vt) - 1) end + v/R
end 
=#
#end
