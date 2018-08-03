"""
Modia module with rotational component models (inspired from Modelica Standard Library.

* Developer: Hilding Elmqvist, Mogram AB  
* Copyright (c) 2016-2018: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)

"""
module Rotational

export Flange, Inertia, Spring, SpringDamper, EMF, IdealGear, Torque, CurrentSensor 
using ..Instantiation
using ..Electric
using ..Blocks

using Modia

@model Flange begin
  phi=Var(T=Float64, size=())
  tau=Float(flow=true, size=())
end

@model Inertia begin
    # 1D-rotational component with inertia"
  J=Parameter(0, min=0)   # Moment of inertia 

  flange_a=Flange()  # Left flange of shaft
  flange_b=Flange()  # Right flange of shaft

  phi=Float(start=0.0)
  w=Float(start=0.0)
  a=Float()  
@equations begin 
  phi = flange_a.phi
  phi = flange_b.phi
  w = der(phi)
  a = der(w)
  J*a = flange_a.tau + flange_b.tau
  end
end

@model PartialCompliant begin
    # Partial model for the compliant connection of two rotational 1-dim. shaft flanges
  phi_rel=Float()
  tau=Float(size=())
  flange_a=Flange()
  flange_b=Flange()
@equations begin 
  phi_rel = flange_b.phi - flange_a.phi
  flange_b.tau = tau
  flange_a.tau = -tau
  end
end 

@model Spring begin
    # Linear 1D rotational spring
  @extends PartialCompliant()
  @inherits tau, phi_rel
  c=Parameter(min=0, start=1.0e5)       # Spring constant
  phi_rel0=0   # Unstretched spring angle
@equations begin 
  tau = c*(phi_rel - phi_rel0)
  end
end 

@model SpringDamper begin
    # Linear 1D rotational spring
  @extends PartialCompliant(phi_rel=Float(size=(), state=false))
  @inherits tau, phi_rel
  c=Parameter(min=0, start=1.0e5)       # Spring constant
	d=Parameter()
  phi_rel0=0   # Unstretched spring angle
@equations begin 
  tau = c*(phi_rel - phi_rel0) + d*der(phi_rel)
  end
end 

@model EMF begin
    # Electromotoric force (electric/mechanic transformer)
  k=1

  p=Pin()
  n=Pin()
  flange=Flange()
	
  v=Float(size=())
  i=Float(size=())
  phi=Float(state=false, size=())
  w=Float(size=())
@equations begin 
  v = p.v - n.v
  0 = p.i + n.i
  i = p.i

  phi = flange.phi
  w = der(phi)
  k*w = v
  flange.tau = -k*i
  end
end

@model IdealGear begin
    # Ideal gear without inertia
#  @extends Gear()

#  @extends PartialElementaryTwoFlangesAndSupport2()
  flange_a=Flange()  # Left flange of shaft
  flange_b=Flange()  # Right flange of shaft
  ratio=1     # Transmission ratio (flange_a.phi/flange_b.phi)
  phi_a=Float(size=())
  phi_b=Float(size=())
@equations begin 
  phi_a = flange_a.phi
  phi_b = flange_b.phi
  phi_a = ratio*phi_b
  0 = ratio*flange_a.tau + flange_b.tau
  end
end 

@model Torque begin # Input signal acting as external torque on a flange
  flange=Flange()  # Flange of shaft
  tau=Float()
@equations begin 
  flange.tau = -tau
  end
end 

@model RotationalSensor begin
   
end 

@model CurrentSensor begin
    # Sensor to measure the current in a branch
  @extends RotationalSensor()
  p=Pin()
  n=Pin()
  i=Float(size=())
@equations begin 
  p.v = n.v
  p.i = i
  n.i = -i
  end
end 

end