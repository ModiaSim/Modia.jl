"""
Modia module with translational component models (inspired from Modelica Standard Library).

* Developer: Hilding Elmqvist, Mogram AB, Martin Otter, DLR
* Copyright (c) 2021: Hilding Elmqvist, Martin Otter
* License: MIT (expat)

"""
#module Rotational

#export Flange, Inertia, Spring, SpringDamper, EMF, IdealGear, Torque, CurrentSensor, Fixed, Damper, IdealGear_withSupport, SpeedSensor

using Modia

# Connector for 1D translational systems
TranslationalFlange = Model(
    s = potential,  # Absolute position of flange
    f = flow        # Cut force directed into the flange
)

# Flange fixed in housing at a given position
FixedPosition = Model(
    flange = TranslationalFlange,
    equations = :[
        flange.s = 0.0]
)

# Rigid connection of two translational 1D flanges
PartialRigid = Model(
    flange_a = TranslationalFlange,
    flange_b = TranslationalFlange,
    L = 0.0u"m",                # Length of component, from left flange to right flange
    equations = :[
        flange_a.s = s - L/2
        flange_b.s = s + L/2]
)

# 1D-translational mass
TranslationalMass = PartialRigid | Model(
    flange_a = TranslationalFlange, # (info = "Left flange of mass")
    flange_b = TranslationalFlange, # (info = "Right flange of mass")
    m = 1.0u"kg",                   # Mass of sliding mass
    s = Var(init = 0.0u"m"),        # Absolute position of center of component
    v = Var(init = 0.0u"m/s"),
    equations = :[
        v = der(s)
        m*der(v) = flange_a.f + flange_b.f]
)

# Partial input signal acting as external force on a flange
PartialForce = Model(
    f      = input,
    flange = TranslationalFlange
)

# Input signal acting as external torque on a flange
Force         = PartialForce | Model(equations = :[flange.f = -f])
UnitlessForce = PartialForce | Model(equations = :[flange.f = -f*u"N"])


# Partial model to measure a single absolute flange variable
PartialAbsoluteSensor = Model(
    flange    = TranslationalFlange,
    equations = :[flange.f = 0] 
)

# Ideal sensor to measure the absolute flange position and velocity
VelocitySensor         = PartialAbsoluteSensor | Model(v = output, equations = :[v = der(flange.s)])
UnitlessVelocitySensor = PartialAbsoluteSensor | Model(v = output, equations = :[v = der(flange.s)*u"m/s"])

PositionSensor         = PartialAbsoluteSensor | Model(s = output, equations = :[s = flange.s])
UnitlessPositionSensor = PartialAbsoluteSensor | Model(s = output, equations = :[s = flange.s*u"m"])


#end
