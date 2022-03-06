"""
Modia module with rotational component models (inspired from Modelica Standard Library).

* Developer: Hilding Elmqvist, Mogram AB, Martin Otter, DLR
* Copyright (c) 2016-2021: Hilding Elmqvist, Martin Otter
* License: MIT (expat)

"""
#module Rotational

#export Flange, Inertia, Spring, SpringDamper, EMF, IdealGear, Torque, CurrentSensor, Fixed, Damper, IdealGear_withSupport, SpeedSensor

using Modia

# Connector for 1D rotational systems
Flange = Model(
    phi = potential,
    tau = flow
)

# Flange fixed in housing at a given angle
#=
Fixed = Model(
    flange = Flange,
    phi0   = 0.0u"rad",
    equations = :[
        flange.phi = 0.0u"rad"]
)
=#
Fixed = Model(
    flange = Flange,
    equations = :[
        flange.phi = 0.0]
)

# 1D-rotational component with inertia
Inertia = Model(
    flange_a = Flange, # (info = "Left flange of shaft")
    flange_b = Flange, # (info = "Right flange of shaft")
    J = 1.0u"kg*m^2", # (0, min=0, info = "Moment of inertia") #, T = u"kg*m^2")    
    phi = Var(init=0.0u"rad"), 
    w = Var(init=0.0u"rad/s"),
    equations = :[
        phi = flange_a.phi
        phi = flange_b.phi
        w = der(phi)
        a = der(w)
        J * a = flange_a.tau + flange_b.tau ]
)

# Partial model for the compliant connection of two rotational 1-dim. shaft flanges
PartialCompliant = Model(
    flange_a = Flange,
    flange_b = Flange,
    equations = :[
        phi_rel = flange_b.phi - flange_a.phi
        flange_b.tau = tau
        flange_a.tau = -tau ]
)

# Linear 1D rotational spring
Spring = PartialCompliant | Model(
    c = 1.0u"N*m/rad",   # (min = 0, info = "Spring constant")
    phi_rel0 = 0.0u"rad", # (info = "Unstretched spring angle") 
    equations = :[
        tau = c * (phi_rel - phi_rel0) ]
) 

# Linear 1D rotational spring with damper
SpringDamper = PartialCompliant | Model(
    c = 1.0*u"N*m/rad", # (min = 0, info = "Spring constant")
    d = 0.0u"N*m*s/rad", # (info = "Damping constant")
    phi_rel0 = 0.0u"rad", # (info = "Unstretched spring angle")
    equations = :[
        tau = c * (phi_rel - phi_rel0) + d * der(phi_rel) ]
)

# Electromotoric force (electric/mechanic) transformer
EMF = Model(
    k = 1.0u"N*m/A", # (info = "Transformation coefficient")
    p = Pin,
    n = Pin,
    flange = Flange,
  
    equations = :[
        v = p.v - n.v
        0 = p.i + n.i
        i = p.i

        phi = flange.phi
        w = der(phi)
        k * w = v
        flange.tau = -k * i ]
)

# Ideal gear
IdealGear = Model(
    flange_a = Flange, # "Left flange of shaft"
    flange_b = Flange, # "Right flange of shaft"
    ratio    = 1.0, # "Transmission ratio (flange_a.phi/flange_b.phi)"
    equations = :[
        phi_a = flange_a.phi
        phi_b = flange_b.phi
        phi_a = ratio * phi_b
        0 = ratio * flange_a.tau + flange_b.tau ]
) 

# Ideal gear with support
IdealGear_withSupport = Model(
    flange_a = Flange, # "Left flange of shaft"
    flange_b = Flange, # "Right flange of shaft"
    support  = Flange, # "Support flange"
    ratio    = 1.0,  # "Transmission ratio"
    equations = :[
        phi_a = flange_a.phi - support.phi
        phi_b = flange_b.phi - support.phi
        phi_a = ratio * phi_b
        0 = ratio * flange_a.tau + flange_b.tau
        0 = flange_a.tau + flange_b.tau + support.tau ]
) 


# Partial input signal acting as external torque on a flange
PartialTorque = Model(
    tau = input,
    flange = Flange
)

# Input signal acting as external torque on a flange
Torque         = PartialTorque | Model(equations = :[flange.tau = -tau])
UnitlessTorque = PartialTorque | Model(equations = :[flange.tau = -tau*u"N*m"])

#=
# Partial model for the compliant connection of two rotational 1-dim. shaft flanges where the relative angle and speed are used as preferred states
PartialCompliantWithRelativeStates = Model(
    flange_a = Flange,
    flange_b = Flange,
    init     = Map(phi_rel=0.0u"rad", w_rel=0.0u"rad/s"),
    equations = :[
        phi_rel      = flange_b.phi - flange_a.phi
        w_rel        = der(phi_rel)
        a_rel        = der(w_rel)
        flange_b.tau = tau
        flange_a.tau = -tau ]
)

# Linear 1D rotational damper
Damper = PartialCompliantWithRelativeStates | Model(
    d = 1.0u"N*m*s/rad", # (info = "Damping constant"),
    equations = :[
        tau = d * w_rel ]
)
=#
# Linear 1D rotational damper
Damper = Model(
    flange_a = Flange,
    flange_b = Flange,
    phi_rel  = Var(start=0.0u"rad"),
    d        = 1.0u"N*m*s/rad", # (info = "Damping constant"),
    equations = :[
        phi_rel      = flange_b.phi - flange_a.phi
        w_rel        = der(phi_rel)
        flange_b.tau = tau
        flange_a.tau = -tau
        tau          = d * w_rel ]
)



# Partial model to measure a single absolute flange variable
PartialAbsoluteSensor = Model(
    flange   = Flange,
    equations = :[flange.tau = 0] 
)

# Ideal sensor to measure the absolute flange angular velocity
AngleSensor         = PartialAbsoluteSensor | Model(phi = output, equations = :[phi = flange.phi])
UnitlessAngleSensor = PartialAbsoluteSensor | Model(phi = output, equations = :[phi = flange.phi*u"1/rad"])

SpeedSensor         = PartialAbsoluteSensor | Model(w = output, equations = :[w = der(flange.phi)])
UnitlessSpeedSensor = PartialAbsoluteSensor | Model(w = output, equations = :[w = der(flange.phi)*u"s/rad"])


#end
