module Pendulum3D_2

println("... Pendulum3D_2 - simple pendulum modeled with multibody components with 3D animation export")

using  Modia
@usingModiaPlot

Pendulum = Model3D(
    world = Object3D(feature=Scene(animationFile="Pendulum2.json")),
    obj1  = Object3D(feature=Solid(shape=Box(lengthX=1.0, lengthY=0.2, lengthZ=0.2),
                solidMaterial="Steel", visualMaterial=VisualMaterial(color="Blue"))),
    obj2  = Object3D(parent=:obj1, feature=Visual(shape=Cylinder(diameter=0.1, length=0.21),
                visualMaterial=VisualMaterial(color="Red")), translation=[-0.5, 0.0, 0.0]),
    rev   = Revolute(obj1=:world, obj2=:obj2)
)

pendulum = @instantiateModel(Pendulum, unitless=true)
simulate!(pendulum, stopTime=3.0)
plot(pendulum, "rev.phi")
end
