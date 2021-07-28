module BouncingSphere3D_2

using Modia
@usingModiaPlot

BouncingSphere = Model(
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
    free = FreeMotion(obj1=:world, obj2=:sphere, r=Var(init=[0.0, 1.0, 0.0]), w=Var(init=[10.0, 0.0, -5.0]))
)

bouncingSphere = @instantiateModel(buildModia3D(BouncingSphere), unitless=true, log=false, logStateSelection=false, logCode=false)


stopTime = 2.7
tolerance = 1e-8
requiredFinalStates = [0.28711931505126853, -0.9780916472966511, -0.20833195055744314, 0.10925258092251605, -3.8123627743313273, -0.9364166852211897, 1.228619049097046, -0.5724547180688829, 0.22566245156307535, -7.694985753129207, 0.3826986158402523, 5.519005772512668]
simulate!(bouncingSphere, stopTime=stopTime, tolerance=tolerance, log=true, logStates=true, logEvents=true, requiredFinalStates=requiredFinalStates)

plot(bouncingSphere, ["free.r" "free.rot"; "free.v" "free.w"], figure=1)

end
