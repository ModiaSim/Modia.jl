module BouncingSphere3D_1

using Modia
using ModiaLang.StaticArrays

BouncingSphere = Model3D(
    boxHeigth = 0.1,
    groundMaterial = VisualMaterial(color="DarkGreen", transparency=0.5),
    gravField = UniformGravityField(g=9.81, n=[0, -1, 0]),
    world = Object3D(feature=Scene(gravityField=:gravField,
                                          visualizeFrames=false,
                                          defaultFrameLength=0.2,
                                          gap=0.2,
                                          enableContactDetection=true,
                                          visualizeContactPoints=false)),
    worldFrame = Object3D(parent=:world, feature=Visual(shape=CoordinateSystem(length=0.5))),
    ground = Object3D(parent=:world,
                      translation=:[0.0,-boxHeigth/2,0.0],
                      feature=Solid(shape=Box(lengthX=4.0, lengthY=:boxHeigth, lengthZ=0.7),
                                    visualMaterial=:groundMaterial,
                                    solidMaterial="Steel",
                                    collision=true)),
    sphere = Object3D(feature=Solid(shape=Sphere(diameter=0.2),
                                    visualMaterial=VisualMaterial(color="Blue"),
                                    solidMaterial="Steel",
                                    massProperties=MassPropertiesFromShapeAndMass(mass=0.001),
                                    collision=true)),
    free = FreeMotion(obj1=:world, obj2=:sphere, r=Var(init=SVector{3,Float64}(0.0, 1.0, 0.0)))
)

bouncingSphere = @instantiateModel(BouncingSphere, unitless=true, log=false, logStateSelection=false, logCode=false)

#@show bouncingSphere.parameterExpressions
#@show bouncingSphere.parameters

stopTime = 2.2
dtmax = 0.1
tolerance = 1e-8
requiredFinalStates = [0.0, 0.10092547226369847, 0.0, 0.0, 0.01950941258387679, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
simulate!(bouncingSphere, stopTime=stopTime, tolerance=tolerance, dtmax=dtmax, log=true, logStates=true, logEvents=true, 
          requiredFinalStates_atol=0.01, requiredFinalStates=requiredFinalStates)

@usingModiaPlot
plot(bouncingSphere, "free.r", figure=1)

end
