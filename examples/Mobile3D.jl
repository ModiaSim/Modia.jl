module Mobile3D

using Modia

# Modia equation-based models
include("$(Modia.modelsPath)/AllModels.jl")


const depthMax = 6
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
Mobile = Model(
    world      = Object3D(feature=Scene(gravityField=UniformGravityField(g=9.81, n=[0, -1, 0]), enableContactDetection=false, nominalLength=0.2*barLength(depthMax), enableVisualization=enableVisualization, animationFile=animationFile)),
    worldFrame = Object3D(parent = :world, feature = Visual(shape = CoordinateSystem(length=rodLength))),
    top        = createMobile(depthMax),
    rev0       = RevoluteWithDamping(obj1 = :world, obj2 = :(top.rod.frame1), phi_start=0.2)
)

mobile = @instantiateModel(buildModia3D(Mobile), unitless=true, log=false, logDetails=false, logModel=false, logStateSelection=false, logCode=false, logExecution=false, logTiming=false)

stopTime = 5.0
tolerance = 1e-5
requiredFinalStates = [-0.009569246139092647, 0.06206606324720822, 0.0010219666026627296, 0.024678616817829752, 0.004218575209138456, 0.028468399940849164, -0.009906941842068384, 0.04296850265297842, 0.0007318009858593783, 0.02017149964438393, 0.003687261686657736, 0.023502975932150394, -0.0448394459612976, 0.026723967721919494, 0.0537435564833799, 0.12331025494664288, 0.05248795026088127, 0.11675479912295675, -0.009264844274058073, 0.06728855030873694, 0.0011992360615969712, 0.025406524349364322, 0.004501037229020791, 0.02917639504873484, -0.00965181608855167, 0.0477780102217653, 0.0008958989126344595, 0.020867265917878713, 0.003946827273570862, 0.024201038339941326, -0.04462421281684384, 0.03746102813853381, 0.05486746017952349, 0.11577048656288175, 0.05363002444778805, 0.11033674247191452, -0.02865814302511715, -0.46249723347164023, 0.02588030023351807, 0.562896850566723, 0.02821667734365895, 0.5623903988440904, -0.009403295780267986, 0.06162802647976917, 0.0010701734591814242, 0.024348108279934208, 0.004222449370467353, 0.028064091478322103, -0.00978017118379629, 0.042937738418303645, 0.0007712669619712274, 0.02000389926557959, 0.0036870480729875752, 0.02329559286742148, -0.04450699270608273, 0.027584416181856048, 0.053512955516200826, 0.1203910664037745, 0.052242780723209266, 0.1148698118645492, -0.00908973848103805, 0.06683238262833087, 0.0012488875262191405, 0.025054556524096706, 0.004506535939791427, 0.02874309504494291, -0.009516262851887764, 0.047741100768557164, 0.0009367831965132963, 0.0206820605278372, 0.003948288304262993, 0.02396953329142854, -0.0442740599719771, 0.03835719168973073, 0.054628209624966825, 0.1126450558374131, 0.053378321393267744, 0.10825519520092068, -0.029248882022186243, -0.45868389528398645, 0.02661598555301995, 0.5585200081021191, 0.02895380855666197, 0.5578643900507765, 0.11057838176009542, -0.4328898557101326, -0.1137607979103064, 0.5023169087131035, -0.11293209024062373, 0.49924136851309536, -0.00949264483432839, 0.06269840391598419, 0.0010594569792564434, 0.024662133036232943, 0.004286395100088604, 0.02844710490443639, -0.009858798739817745, 0.04356190104934767, 0.0007615897059759054, 0.020185339490969664, 0.0037442343747962708, 0.023525248558572666, -0.04490190529834802, 0.02861187090459901, 0.054093472471048135, 0.12121222952180354, 0.05282381817075915, 0.11524012664035016, -0.00918252751200148, 0.06794970842859763, 0.0012381041193024645, 0.02538680700446271, 0.004571342380041327, 0.029148124149157, -0.009598363758120855, 0.048403550630416084, 0.0009271317194699239, 0.020879696023068203, 0.00400620925514733, 0.02421928779134631, -0.044678930789653376, 0.039416625110429974, 0.05521797625599808, 0.11354470175688848, 0.05396904771937769, 0.1087130515149207, -0.029386149298980365, -0.4625396185956002, 0.026768483822850988, 0.5631828688184961, 0.02911172224703885, 0.5626015954247693, -0.009324704028879006, 0.062258112033935976, 0.0011079870888003323, 0.024328333168236824, 0.004290387257149793, 0.028037530012991814, -0.009730064979351052, 0.043535748038149275, 0.0008014306498257068, 0.02001640392846844, 0.0037442737008724775, 0.023314890388120756, -0.04456469081149338, 0.0294911482051167, 0.05385917279095154, 0.11825295160221237, 0.052576102903529355, 0.11332114444600068, -0.00900535540075688, 0.0674898715499478, 0.0012880524022445207, 0.02503113089861741, 0.0045769255979494095, 0.02870901282603822, -0.009460770397544161, 0.04837019570880829, 0.0009683769175080957, 0.020692819349648818, 0.004007904268091745, 0.023984347424616213, -0.04432382727176435, 0.04033007578307317, 0.054974659993453495, 0.11037837031847329, 0.05371447828202559, 0.10659649333068186, -0.029979463641433897, -0.4586694882625391, 0.027508572008452217, 0.5587394516064685, 0.02985302711571344, 0.5580073793505101, 0.10994893558074191, -0.43345695851264765, -0.113030926843198, 0.5031313427136914, -0.11219753434123993, 0.5000032205945514, 0.17883870418449752, -0.3672711553933069, -0.17900124607493148, 0.423658663329383, -0.17837270297492822, 0.424290817579418, 0.011456378191880315, 0.3654620146635912]
simulate!(mobile, stopTime=stopTime, tolerance=tolerance, log=true, requiredFinalStates=requiredFinalStates)

@usingModiaPlot
plot(mobile, "rev0.rev.phi")

end
