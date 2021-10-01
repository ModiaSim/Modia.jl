export Object3D, Scene, Visual, Solid
export Box, Beam, Cylinder, Sphere, Ellipsoid
export Cone, Capsule, GearWheel, Grid, SpringShape
export CoordinateSystem, FileMesh
export Font, TextShape
export VisualMaterial
export MassProperties, MassPropertiesFromShape
export MassPropertiesFromShapeAndMass
export UniformGravityField
export RefPath, ptpJointSpace, scheduleReferenceMotion
export calculateRobotMovement
export getRefPathPosition, getRefPathInitPosition, getVariables
export multibodyResiduals!, setModiaJointVariables!
export Fix
export Revolute, RevoluteWithFlange
export Prismatic, PrismaticWithFlange
export J123, J132, FreeMotion

export buildModia3D

import Modia3D
using ModiaLang

# ModiaLang models
include("$(ModiaLang.path)/models/Blocks.jl")
include("$(ModiaLang.path)/models/Electric.jl")
include("$(ModiaLang.path)/models/Rotational.jl")
include("$(ModiaLang.path)/models/Translational.jl")

include("$(Modia3D.path)/src/ModiaInterface/buildModia3D.jl")

Object3D(        ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Object3D), _path = true, kwargs...)
Scene(           ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Scene)                 , kwargs...)
Visual(          ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Visual)                , kwargs...)
Solid(           ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Solid)                 , kwargs...)
Box(             ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Box)                   , kwargs...)
Beam(            ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Beam)                  , kwargs...)
Cylinder(        ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Cylinder)              , kwargs...)
Sphere(          ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Sphere)                , kwargs...)
Ellipsoid(       ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Ellipsoid)             , kwargs...)
Capsule(         ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Capsule)               , kwargs...)
Cone(            ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Cone)                  , kwargs...)
SpringShape(     ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Spring)                , kwargs...)
GearWheel(       ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.GearWheel)             , kwargs...)
Grid(            ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Grid)                  , kwargs...)
VisualMaterial(  ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.VisualMaterial)        , kwargs...)
MassProperties(  ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.MassProperties)        , kwargs...)
CoordinateSystem(; kwargs...) = Par(; _constructor = :(Modia.Modia3D.CoordinateSystem)      , kwargs...)
FileMesh(        ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.FileMesh)              , kwargs...)
Font(            ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Font)                  , kwargs...)
TextShape(       ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.TextShape)             , kwargs...)
Fix(             ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Fix)                   , kwargs...)

MassPropertiesFromShape()  = Par(; _constructor = :(Modia.Modia3D.MassPropertiesFromShape))
MassPropertiesFromShapeAndMass(;mass) = Par(; _constructor = :(Modia.Modia3D.MassPropertiesFromShapeAndMass), mass = mass)
UniformGravityField(; kwargs...) = Par(; _constructor = :(Modia.Modia3D.UniformGravityField), kwargs...)


RefPath(; kwargs...) = Modia3D.ReferencePath(; kwargs...)

ptpJointSpace(; kwargs...) = Modia3D.ptpJointSpace(; kwargs...)
scheduleReferenceMotion(; kwargs...) = Modia3D.scheduleReferenceMotion(; kwargs...)

calculateRobotMovement(args...) = Modia3D.calculateRobotMovement(args...)
getRefPathPosition(args...) = Modia3D.getRefPathPosition(args...)
getRefPathInitPosition(args...) = Modia3D.getRefPathInitPosition(args...)

getVariables(args...) = (args...,)

multibodyResiduals!(args...)     = Modia3D.multibodyResiduals!(args...)
setModiaJointVariables!(args...) = Modia3D.setModiaJointVariables!(args...)

Revolute(; obj1, obj2, axis=3, phi=Var(init=0.0), w=Var(init=0.0), canCollide=true) = Model(; _constructor = Par(value = :(Modia.Modia3D.Revolute), _path = true, ndof = 1),
    obj1 = Par(value = obj1),
    obj2 = Par(value = obj2),
    axis = Par(value = axis),
    canCollide = Par(value = canCollide),
    phi  = phi,
    w    = w,
    equations = :[
        w   = der(phi)
        qdd = der(w)   # standardized name for the generalized joint accelerations
        variables = getVariables(phi, w, 0.0) # standardized name for the generalized joint position, velocity, force
        ]
)

RevoluteWithFlange(; obj1, obj2, axis=3, phi=Var(init=0.0), w=Var(init=0.0), canCollide=true) = Model(; _constructor = Par(value = :(Modia.Modia3D.Revolute), _path = true, ndof = 1),
    obj1   = Par(value = obj1),
    obj2   = Par(value = obj2),
    axis   = Par(value = axis),
    canCollide = Par(value = canCollide),
    flange = Flange,
    phi    = phi,
    w      = w,
    equations = :[
        phi = flange.phi
        w   = der(phi)
        qdd = der(w)   # standardized name for the generalized joint accelerations
        variables = getVariables(phi, w, flange.tau) # standardized name for the generalized joint position, velocity, force
        ]
)

Prismatic(; obj1, obj2, axis=1, s=Var(init=0.0), v=Var(init=0.0), canCollide=true) = Model(; _constructor = Par(value = :(Modia.Modia3D.Prismatic), _path = true, ndof = 1),
    obj1 = Par(value = obj1),
    obj2 = Par(value = obj2),
    axis = Par(value = axis),
    canCollide = Par(value = canCollide),
    s    = s,
    v    = v,
    equations = :[
        v   = der(s)
        qdd = der(v)   # standardized name for the generalized joint accelerations
        variables = getVariables(s, v, 0.0) # standardized name for the generalized joint position, velocity, force
        ]
)

PrismaticWithFlange(; obj1, obj2, axis=1, s=Var(init=0.0), v=Var(init=0.0), canCollide=true) = Model(; _constructor = Par(value = :(Modia.Modia3D.Prismatic), _path = true, ndof = 1),
    obj1   = Par(value = obj1),
    obj2   = Par(value = obj2),
    axis   = Par(value = axis),
    canCollide = Par(value = canCollide),
    flange = TranslationalFlange,
    s      = s,
    v      = v,
    equations = :[
        s   = flange.s
        v   = der(s)
        qdd = der(v)   # standardized name for the generalized joint accelerations
        variables = getVariables(s, v, flange.f) # standardized name for the generalized joint position, velocity, force
        ]
)

"""
    J = J123(rot123::AbstractVector)

Return joint rot. kinematics matrix `J` for Cardan angles `rot123` (rotation sequence x-y-z).
"""
function J123(rot123::AbstractVector)

    (sbe, cbe) = sincos(rot123[2])
    (sga, cga) = sincos(rot123[3])
    return [     cga     -sga  0.0 ;
             cbe*sga  cbe*cga  0.0 ;
            -sbe*cga  sbe*sga  cbe ] / cbe

end

"""
    J = J132(rot132::AbstractVector)

Return joint rot. kinematics matrix `J` for Cardan angles `rot132` (rotation sequence x-z-y).
"""
function J132(rot132::AbstractVector)

    (sga, cga) = sincos(rot132[2])
    (sbe, cbe) = sincos(rot132[3])
    return [ cbe      0.0  sbe     ;
            -sbe*cga  0.0  cbe*cga ;
             cbe*sga  cga  sbe*sga ] / cga

end

FreeMotion(; obj1, obj2, r=Var(init=zeros(3)), rot=Var(init=zeros(3)), v=Var(init=zeros(3)), w=Var(init=zeros(3))) = Model(; _constructor = Par(value = :(Modia.Modia3D.FreeMotion), _path = true, ndof = 6),
    obj1 = Par(value = obj1),
    obj2 = Par(value = obj2),
    r    = r,
    rot  = rot,
    v    = v,
    w    = w,
    equations = :[
        der(r)   = v
        der(rot) = J123(rot) * w
        der(v) = qdd[1:3]
        der(w) = qdd[4:6]
        variables = getVariables(r, rot, v, w)
        ]
)
