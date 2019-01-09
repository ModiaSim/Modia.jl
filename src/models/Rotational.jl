"""
Modia module with rotational component models (inspired from Modelica Standard Library).

* Developer: Hilding Elmqvist, Mogram AB  
* Copyright (c) 2016-2018: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)

"""
module Rotational

export Flange, Inertia, Spring, SpringDamper, EMF, IdealGear, Torque, CurrentSensor 
using ..Electric
using ..Blocks
using Unitful

using Modia

"Rotational angle variable"
Angle(; args...) = Variable(; start = 0.0, size = (), 
                       T = u"rad", info = "Rotational angle", args...)
TorqueVar(; args...) = Variable(; start = 0.0, size = (), 
                           T = u"N*m", info = "Torque", args...)

"""
Connector for 1D rotational systems
"""
@model Flange begin
    phi = Angle() 
    tau = TorqueVar(flow = true)
end

"""
1D-rotational component with inertia
"""
@model Inertia begin
    J = Parameter(0, min=0, info = "Moment of inertia", T = u"kg*m^2")

    flange_a = Flange(info = "Left flange of shaft")
    flange_b = Flange(info = "Right flange of shaft")

    phi = Angle()
    w = Var(start = 0.0, info = "Angular velocity", T = u"rad/s")
    a = Var(info = "Angular acceleration", T = u"rad/s^2")
    @equations begin 
        phi = flange_a.phi
        phi = flange_b.phi
        w = der(phi)
        a = der(w)
        J * a = flange_a.tau + flange_b.tau
    end
end

"""
Partial model for the compliant connection of two rotational 1-dim. shaft flanges
"""
@model PartialCompliant begin
    phi_rel = Angle()
    tau = TorqueVar()
    flange_a = Flange(info = "Left flange of shaft")
    flange_b = Flange(info = "Right flange of shaft")
    @equations begin 
        phi_rel = flange_b.phi - flange_a.phi
        flange_b.tau = tau
        flange_a.tau = -tau
    end
end 

"""
Linear 1D rotational spring
"""
@model Spring begin
    c = Parameter(min = 0, start = 1.0e5, info = "Spring constant", T = u"N*m/rad")
    phi_rel0 = Parameter(0.0, start = 0.0, info = "Unstretched spring angle", T = u"rad")
    @extends PartialCompliant()
    @inherits tau, phi_rel
    @equations begin 
        tau = c * (phi_rel - phi_rel0)
    end
end 

"""
Linear 1D rotational spring with damper
"""
@model SpringDamper begin
    c = Parameter(min = 0, start = 1.0e5, info = "Spring constant", T = u"N*m/rad")
    d = Parameter(info = "Damping constant", T = u"N*m*s/rad")
    phi_rel0 = Parameter(0.0, start = 0.0, info = "Unstretched spring angle", T = u"rad")
    @extends PartialCompliant(phi_rel=Angle(size=(), state=false))
    @inherits tau, phi_rel
    @equations begin 
        tau = c * (phi_rel - phi_rel0) + d * der(phi_rel)
    end
end 

"""
Electromotoric force (electric/mechanic) transformer
"""
@model EMF begin
    k = Parameter(1.0, info = "Transformation coefficient", T = u"N*m/A")

    p = Pin(info = "Positive pin")
    n = Pin(info = "Negative pin")
    flange = Flange(info = "Support/housing of the EMF shaft")
  
    v = Voltage()
    i = Current()
    phi = Angle(state = false)
    w = Variable(T = u"rad/s", info = "Angular velocity")
    @equations begin 
        v = p.v - n.v
        0 = p.i + n.i
        i = p.i

        phi = flange.phi
        w = der(phi)
        k * w = v
        flange.tau = -k * i
    end
end

"""
Ideal gear without inertia
"""
@model IdealGear begin
    # @extends Gear()

    # @extends PartialElementaryTwoFlangesAndSupport2()
    flange_a = Flange(info = "Left flange of shaft")
    flange_b = Flange(info = "Right flange of shaft")
    ratio = Parameter(1.0, info = "Transmission ratio (flange_a.phi/flange_b.phi)")
    phi_a = Angle(info = "Angle of the left flange")
    phi_b = Angle(info = "Angle of the right flange")
    @equations begin 
        phi_a = flange_a.phi
        phi_b = flange_b.phi
        phi_a = ratio * phi_b
        0 = ratio * flange_a.tau + flange_b.tau
    end
end 

"""
Input signal acting as external torque on a flange
"""
@model Torque begin
    tau = TorqueVar()
    flange = Flange()
    @equations begin 
        flange.tau = -tau
    end
end 

@model RotationalSensor begin
   
end 

"""
Sensor to measure the current in a branch
"""
@model CurrentSensor begin
    p = Pin(info = "Positive pin")
    n = Pin(info = "Negative pin")
    i = Current()
    @extends RotationalSensor()
    @equations begin 
        p.v = n.v
        p.i = i
        n.i = -i
    end
end 

end
