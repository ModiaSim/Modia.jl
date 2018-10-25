#
# Initial author: Martin Otter, DLR-SR (first version: Jan. 14, 2017)
# License: MIT (expat)



"""
   module Media

Test module to evaluate how complex media could be implemented in (pure) Julia and then used in Modia.
"""
module Media

# Definition of medium states
abstract type MediumState end

mutable struct State_pT <: MediumState
   p::Float64
   T::Float64
end

mutable struct State_ph <: MediumState
   p::Float64
   h::Float64
end


# Definition of generic Medium and its functions
abstract type AbstractMedium end
density(               m::AbstractMedium,state::MediumState) = error("... function not defined for $m")
specificEnthalpy(      m::AbstractMedium,state::MediumState) = error("... function not defined for $m")
specificInternalEnergy(m::AbstractMedium,state::MediumState) = error("... function not defined for $m")

setState_pT(           m::AbstractMedium,p,T)                = error("... function not defined for $m")
setState_ph(           m::AbstractMedium,p,h)                = error("... function not defined for $m")

end



"""
    module SimpleWater

Simple medium model of water
"""
module SimpleWater
   import ..Media

   # Constants of medium that have the same value for every medium instance
   const cp_const = 4184
   const cv_const = 4184
   const d_const  = 995.586
   const T0       = 273.15

   # Constants of medium specific to a medium instance
   struct Medium <: Media.AbstractMedium
      h_offset::Float64
      Medium(;h_offset=0.0) = new(h_offset)
   end

   # Functions of medium
   Media.density(               m::Medium, state::Media.State_pT) = d_const
   Media.specificEnthalpy(      m::Medium, state::Media.State_pT) = cp_const*(state.T - T0) + m.h_offset
   Media.specificInternalEnergy(m::Medium, state::Media.State_pT) = cv_const*(state.T - T0)
  
   Media.setState_pT(m::Medium, p, T)::Media.State_pT = Media.State_pT(p,T)
   Media.setState_ph(m::Medium, p, h)::Media.State_pT = Media.State_pT(p,T0+(h-m.h_offset)/cp_const)
end



"""
    module Test_Media

Test media models
"""
module Test_Media

import ..Media
import ..SimpleWater

#=
@model FluidPort begin
  medium
  p=Float()
  m_flow=Float(flow=true)
end
=#

# Set a medium at one connector
medium = SimpleWater.Medium(h_offset=10.0)

# Propagate medium in connector

# Inside a model
state = Media.setState_ph(medium, 1e5, 300.0)   # p,h from connector
d     = Media.density(medium,state)
h     = Media.specificEnthalpy(medium,state)


println("d = ", d)
println("h = ", h)


end
