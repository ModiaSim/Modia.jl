module Mobile3D

using Modia

# Modia equation-based models
include("$(Modia.modelsPath)/AllModels.jl")


const depthMax = 5
const enableVisualization = true
const animationFile = string("Mobile", depthMax, ".json")

const rodLength       = 0.1
const rodDiameter     = 0.05*rodLength
const rodMaterial     = "DryWood"
const rodColor        = VisualMaterial(color="Brown", transparency=0.5)
const rodTranslation1 = [0.0,  rodLength/2, 0.0]
const rodTranslation2 = [0.0, -rodLength/2, 0.0]

const sphereDiameter  = 0.3*rodLength
const sphereMaterial  = "DryWood"
const sphereColor     = VisualMaterial(color="Blue")

const barWidth       = 0.3*rodLength
const barThickness   = 0.3*rodLength
const barMaterial    = "DryWood"
const barColor       = VisualMaterial(color="DeepSkyBlue3", transparency=0.5)
const leaveBarLength = 0.8*rodLength

const revDiameter = rodDiameter
const revLength   = 1.2*barThickness
const revColor    = VisualMaterial(color="Red")

const damping = 0.006

Rod = Model(
    frame0 = Object3D(feature = Solid(shape = Cylinder(axis=2, length=rodLength, diameter=rodDiameter),
                                    solidMaterial  = rodMaterial,
                                    visualMaterial = rodColor)),
    frame1 = Object3D(parent = :frame0, translation = rodTranslation1),
    frame2 = Object3D(parent = :frame0, translation = rodTranslation2)
)
Bar = Model(
    L = missing,
#    frame0 = Object3D(feature = Solid(shape = Beam(axis = 1, length = :L, width=barWidth, thickness=barThickness),
    frame0 = Object3D(feature = Solid(shape = Box(lengthX=:(L + barWidth), lengthY=barWidth, lengthZ=barThickness),
                                      solidMaterial  = barMaterial,
                                      visualMaterial = barColor)),
    frame1 = Object3D(parent = :frame0, translation = :[-L/2, 0.0, 0.0]),
    frame2 = Object3D(parent = :frame0, translation = :[ L/2, 0.0, 0.0])
)
RevoluteWithDamping(;obj1, obj2, axis=3, phi_start=0.0) = Model(
    rev     = RevoluteWithFlange(obj1=obj1, obj2=obj2, axis=axis, phi=Var(init=phi_start)),
    cyl     = Object3D(parent  = obj1,
                    feature = Visual(shape = Cylinder(axis=axis, diameter=revDiameter, length=revLength),
                                        visualMaterial = revColor)),
    damper  = Damper | Map(d=damping),
    fixed   = Fixed,
    connect = :[(damper.flange_b, rev.flange),
                (damper.flange_a, fixed.flange)]
)
function barLength(depth)
    if depth<=2
        return leaveBarLength
    end
    return 2*barLength(depth-1)
end
function createMobile(depth)
    if depth == 1
        Model(
            rod    = Rod,
            sphere = Object3D(parent = :(rod.frame0), translation = rodTranslation2,
                              feature = Solid(shape = Sphere(diameter=sphereDiameter),
                                              solidMaterial  = sphereMaterial,
                                              visualMaterial = sphereColor))
        )
    else
        Model(
            rod  = Rod,
            bar  = Bar | Map(L=barLength(depth)),
            sub1 = createMobile(depth-1),
            sub2 = createMobile(depth-1),
            rev0 = RevoluteWithDamping(obj1=:(rod.frame2), obj2=:(bar.frame0)),
            rev1 = RevoluteWithDamping(obj1=:(bar.frame1), obj2=:(sub1.rod.frame1)),
            rev2 = RevoluteWithDamping(obj1=:(bar.frame2), obj2=:(sub2.rod.frame1))
        )
    end
end
Mobile = Model3D(
    world      = Object3D(feature=Scene(gravityField=UniformGravityField(g=9.81, n=[0, -1, 0]), enableContactDetection=false, nominalLength=0.2*barLength(depthMax), enableVisualization=enableVisualization, animationFile=animationFile)),
    worldFrame = Object3D(parent = :world, feature = Visual(shape = CoordinateSystem(length=rodLength))),
    top        = createMobile(depthMax),
    rev0       = RevoluteWithDamping(obj1 = :world, obj2 = :(top.rod.frame1), phi_start=0.2)
)

mobile = @instantiateModel(Mobile, unitless=true, log=false, logDetails=false, logModel=false, logStateSelection=false, logCode=false, logExecution=false, logTiming=false)

stopTime = 5.0
tolerance = 1e-5
requiredFinalStates = [0.011686313703318246, 0.01045668408948015, 0.0019009731577780986, -0.01786066760079631, 0.0038231810505064317, -0.02016203637563375, 0.0126793101749162, 0.01691192217141304, 0.0023271830236842326, -0.017388317974414997, 0.004235202194998209, -0.019458842081769283, 0.02369063170860951, 0.15558539617677464, -0.004928234187014311, -0.24357952222975443, -0.0008598315100012872, -0.25778191880090334, 0.011970016614214521, 0.009770151879848948, 0.0019715145725989034, -0.018428234746089393, 0.003927329453913458, -0.020712258229027464, 0.0130019833269834, 0.01646408641685571, 0.002416123119856374, -0.017925382312598832, 0.004356974507478808, -0.01996955076328113, 0.024253081005718227, 0.15740983132693404, -0.00508053867112354, -0.24798117927445137, -0.0008446712535811513, -0.26259974654104284, -0.04547737419362646, 0.262258261366544, 0.06092051704639535, -0.2992009041652458, 0.062163705985459186, -0.30201328017123014, 0.011774841775337119, 0.008484715801774874, 0.0017813540404277002, -0.018359678091449557, 0.003722833439382662, -0.020559247090548637, 0.012811196185191788, 0.014445725473740444, 0.002212787855986135, -0.018108079678312596, 0.00414068450155657, -0.02006172527918578, 0.024682227696403755, 0.15164227948198014, -0.0064617212778372585, -0.2427675425787477, -0.002449157908043746, -0.2581724674241067, 0.012050676208440915, 0.00785647627890249, 0.0018526658835679355, -0.01890009720306999, 0.003825074776328946, -0.02107627115366446, 0.013126442724564754, 0.014047087346111506, 0.0023025518844593472, -0.018620796297829156, 0.004260695947219754, -0.020542119547756133, 0.025232373917168024, 0.15350921040211213, -0.006601205580033179, -0.24702565773703974, -0.0024230921741892745, -0.2628525384776904, -0.0440032853139086, 0.269398657305908, 0.05921435657256203, -0.3091722470155659, 0.0604662087479324, -0.31187411328295117, 0.07674790919904127, 0.17815099357106337, -0.06512247859139587, -0.1963202500586452, -0.06616717655644311, -0.20502267390766496, 0.04770632153276008, -0.1850427647991364]
simulate!(mobile, stopTime=stopTime, tolerance=tolerance, log=true, requiredFinalStates=requiredFinalStates)

@usingModiaPlot
plot(mobile, "rev0.rev.phi")

end
