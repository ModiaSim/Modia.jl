export Object3D, Scene, Visual, Solid
export Box, Beam, Cylinder, Sphere, Ellipsoid
export Cone, Capsule, GearWheel, Grid, SpringShape
export CoordinateSystem, FileMesh
export Font, TextShape, ModelicaShape
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
export J123, J132, J123or132, singularRem, FreeMotion, change_rotSequenceInNextIteration!

export buildModia3D

import Modia3D
using ModiaLang

# ModiaLang models
include("$(ModiaLang.path)/models/Blocks.jl")
include("$(ModiaLang.path)/models/Electric.jl")
include("$(ModiaLang.path)/models/Rotational.jl")
include("$(ModiaLang.path)/models/Translational.jl")

include("$(Modia3D.path)/src/ModiaInterface/buildModia3D.jl")

Object3D(        ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Object3D{FloatType}), _path = true, kwargs...)
Scene(           ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Scene{FloatType})                 , kwargs...)
Visual(          ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Visual)                           , kwargs...)
Solid(           ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Solid{FloatType})                 , kwargs...)
Box(             ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Box{FloatType})                   , kwargs...)
Beam(            ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Beam{FloatType})                  , kwargs...)
Cylinder(        ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Cylinder{FloatType})              , kwargs...)
Sphere(          ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Sphere{FloatType})                , kwargs...)
Ellipsoid(       ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Ellipsoid{FloatType})             , kwargs...)
Capsule(         ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Capsule{FloatType})               , kwargs...)
Cone(            ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Cone{FloatType})                  , kwargs...)
SpringShape(     ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Spring)                           , kwargs...)
GearWheel(       ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.GearWheel)                        , kwargs...)
Grid(            ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Grid)                             , kwargs...)
VisualMaterial(  ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.VisualMaterial)                   , kwargs...)
MassProperties(  ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.MassProperties{FloatType})        , kwargs...)
CoordinateSystem(; kwargs...) = Par(; _constructor = :(Modia.Modia3D.CoordinateSystem)                 , kwargs...)
FileMesh(        ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.FileMesh)                         , kwargs...)
Font(            ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Font)                             , kwargs...)
TextShape(       ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.TextShape)                        , kwargs...)
ModelicaShape(   ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.ModelicaShape)                    , kwargs...)
Fix(             ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Fix{FloatType})                   , kwargs...)
Bushing(         ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.Bushing{FloatType})               , kwargs...)
SpringDamperPtP( ; kwargs...) = Par(; _constructor = :(Modia.Modia3D.SpringDamperPtP{FloatType})       , kwargs...)

MassPropertiesFromShape()              = Par(; _constructor = :(Modia.Modia3D.MassPropertiesFromShape{FloatType}))
MassPropertiesFromShapeAndMass(; mass) = Par(; _constructor = :(Modia.Modia3D.MassPropertiesFromShapeAndMass{FloatType}), mass = mass)
UniformGravityField(; kwargs...)       = Par(; _constructor = :(Modia.Modia3D.UniformGravityField), kwargs...)

RefPath(; kwargs...)                 = Modia3D.ReferencePath(; kwargs...)
ptpJointSpace(; kwargs...)           = Modia3D.ptpJointSpace(; kwargs...)
scheduleReferenceMotion(; kwargs...) = Modia3D.scheduleReferenceMotion(; kwargs...)
calculateRobotMovement(args...)      = Modia3D.calculateRobotMovement(args...)
getRefPathPosition(args...)          = Modia3D.getRefPathPosition(args...)
getRefPathInitPosition(args...)      = Modia3D.getRefPathInitPosition(args...)

getVariables(args...)            = (args...,)
multibodyResiduals!(args...)     = Modia3D.multibodyResiduals!(args...)
setModiaJointVariables!(args...) = Modia3D.setModiaJointVariables!(args...)

Revolute(; obj1, obj2, axis=3, phi=Var(init=0.0), w=Var(init=0.0), canCollide=true) = Model(; _constructor = Par(value = :(Modia.Modia3D.Revolute{FloatType}), _path = true, ndof = 1),
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

RevoluteWithFlange(; obj1, obj2, axis=3, phi=Var(init=0.0), w=Var(init=0.0), canCollide=true) = Model(; _constructor = Par(value = :(Modia.Modia3D.Revolute{FloatType}), _path = true, ndof = 1),
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

Prismatic(; obj1, obj2, axis=1, s=Var(init=0.0), v=Var(init=0.0), canCollide=true) = Model(; _constructor = Par(value = :(Modia.Modia3D.Prismatic{FloatType}), _path = true, ndof = 1),
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

PrismaticWithFlange(; obj1, obj2, axis=1, s=Var(init=0.0), v=Var(init=0.0), canCollide=true) = Model(; _constructor = Par(value = :(Modia.Modia3D.Prismatic{FloatType}), _path = true, ndof = 1),
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



"""
    next_isrot123 = change_rotSequenceInNextIteration!(rot::AbstractVector, isrot123::Bool, instantiatedModel::SimulationModel, x, rot_name)

Change rotation sequence of `rot` from `x-axis, y-axis, z-axis` to `x-axis, z-axis, y-axis` or visa versa in the next event iteration:

- If `isrot123 = true`, return `next_isrot123 = false` and `x[..] = rot132fromR(Rfromrot123(rot))`

- If `isrot123 = false`, return `next_isrot123 = true` and `x[..] = rot123fromR(Rfromrot132(rot))`
"""
function change_rotSequenceInNextIteration!(rot::AbstractVector, isrot123::Bool, instantiatedModel::SimulationModel, x, rot_name)::Bool
    if isrot123
        #println("        switch $rot_name 123 -> 132")
        next_rot      = Modia3D.rot132fromR(Modia3D.Rfromrot123(rot))
        next_isrot123 = false
    else
        #println("        switch $rot_name 132 -> 123")
        next_rot      = Modia3D.rot123fromR(Modia3D.Rfromrot132(rot))
        next_isrot123 = true
    end

    # Change x-vector with the next_rot values
    eqInfo = instantiatedModel.equationInfo
    startIndex = eqInfo.x_info[ eqInfo.x_dict[rot_name] ].startIndex
    x[startIndex]   = next_rot[1]
    x[startIndex+1] = next_rot[2]
    x[startIndex+2] = next_rot[3]
    return next_isrot123
end


singularRem(ang) = abs(rem2pi(ang, RoundNearest)) - 1.5  # is negative/positive in valid/singular angle range
J123or132(rot, isrot123) = isrot123 ? J123(rot) : J132(rot)


FreeMotion(; obj1, obj2, r=Var(init=zeros(3)), rot=Var(init=zeros(3)), v=Var(init=zeros(3)), w=Var(init=zeros(3))) = Model(; _constructor = Par(value = :(Modia.Modia3D.FreeMotion{FloatType}), _path = true, ndof = 6),
    obj1 = Par(value = obj1),
    obj2 = Par(value = obj2),
    r    = r,
    rot  = rot,
    v    = v,
    w    = w,

    next_isrot123 = Var(start=true),
    _rotName = "???",  # is changed by buildModia3D to the full path name of "rot"

    equations = :[
        der(r) = v

        isrot123 = pre(next_isrot123)
        rot2_singularity = positive(singularRem(rot[2]))
        next_isrot123 = if rot2_singularity; change_rotSequenceInNextIteration!(rot, isrot123, instantiatedModel, _x, _rotName) else isrot123 end
        der(rot) = J123or132(rot,isrot123) * w

        der(v) = qdd[1:3]
        der(w) = qdd[4:6]
        variables = getVariables(r, rot, v, w, isrot123)
        ]
)
