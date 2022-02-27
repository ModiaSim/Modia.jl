module BouncingSphere3D_2

using Modia
using ModiaLang.StaticArrays
@usingModiaPlot

BouncingSphere = Model3D(
    boxHeigth = 0.05,
    groundMaterial = VisualMaterial(color="DarkGreen", transparency=0.5),
    sphereMaterial = VisualMaterial(color="Red"),
    gravField = UniformGravityField(g=9.81, n=[0, -1, 0]),
    world = Object3D(feature=Scene(gravityField=:gravField,
                                   visualizeFrames=false, defaultFrameLength=0.2, gap=0.2,
                                   enableContactDetection=true, visualizeContactPoints=false)),
    worldFrame = Object3D(parent=:world, feature=Visual(shape=CoordinateSystem(length=0.5))),
    ground = Object3D(parent=:world, translation=:[0.0, -boxHeigth/2, 0.0],
                      feature=Solid(shape=Box(lengthX=0.8, lengthY=:boxHeigth, lengthZ=0.5),
                                    visualMaterial=:groundMaterial, solidMaterial="Steel",
                                    collision=true)),
    bottom = Object3D(parent=:world, translation=:[0.0, -0.4-boxHeigth/2, 0.5],
                      feature=Solid(shape=Box(lengthX=0.8, lengthY=:boxHeigth, lengthZ=0.5),
                                    visualMaterial=:groundMaterial, solidMaterial="Steel",
                                    collision=true)),
    wall = Object3D(parent=:world, translation=:[0.0, -0.2, 0.75+boxHeigth/2],
                    feature=Solid(shape=Box(lengthX=0.8, lengthY=0.4, lengthZ=:boxHeigth),
                                  visualMaterial=:groundMaterial, solidMaterial="Steel",
                                  collision=true)),
    sphere = Object3D(feature=Solid(shape=Sphere(diameter=0.2),
                                    visualMaterial=:sphereMaterial, solidMaterial="Steel",
                                    collision=true)),
    free = FreeMotion(obj1=:world, obj2=:sphere, r=Var(init=SVector{3,Float64}(0.0, 1.0, 0.0)), w=Var(init=SVector{3,Float64}(10.0, 0.0, -5.0)))
)

bouncingSphere = @instantiateModel(BouncingSphere, unitless=true, log=false, logStateSelection=false, logCode=false)


stopTime = 2.7
tolerance = 1e-8
requiredFinalStates = [0.29302652789657424, -1.2392524758992525, -0.13171402091799628, 0.11357140027502943, -4.417276610857313, -0.8729857092038694, 2.3810648360145152, -0.426312211044484 , 0.13656197405807574, -7.335181005993924, 0.812877179123953, 4.870682494334355]
simulate!(bouncingSphere, stopTime=stopTime, tolerance=tolerance, log=true, logStates=true, logEvents=true, 
          requiredFinalStates_rtol = 0.2, requiredFinalStates_atol = 0.2, requiredFinalStates=requiredFinalStates)

plot(bouncingSphere, ["free.r" "free.rot"; "free.v" "free.w"], figure=1)

end
