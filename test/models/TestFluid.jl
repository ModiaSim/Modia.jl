module TestFluid  #################################################################

module MyModiaMedia #################################################################

export AbstractMedium
export density
export density_der_1, density_der_2

abstract type AbstractMedium end

struct Water <: AbstractMedium
    baseProperties::Symbol
    singleState::Bool
    data::Int64
    Water(data) = new(:BaseProperties_Water,false,data)
end

struct Air <: AbstractMedium
    baseProperties::Symbol
    singleState::Bool
    data::Float64
    extra::Float64
    Air(data,extra) = new(:BaseProperties_Air,false,data,extra)
end

struct SimpleMedium <: AbstractMedium
    baseProperties::Symbol
    singleState::Bool
    d0::Float64
    SimpleMedium(d0) = new(:BaseProperties_SimpleMedium,true,d0)
end

density(Medium::AbstractMedium, p) = error("\nError from density: No Medium defined\n")
density(Medium::Water         , p) = Medium.data*p
density(Medium::Air           , p) = Medium.data*p + Medium.extra
density(Medium::SimpleMedium  , p) = Medium.d0
density(Medium::SimpleMedium)      = Medium.d0

density_der_1(Medium::AbstractMedium, p) = 0
density_der_1(Medium::AbstractMedium)    = 0

density_der_2(Medium::Water, p) = Medium.data
density_der_2(Medium::Air  , p) = Medium.data

end
using  .MyModiaMedia

println("TestFluid: Tests features needed for Fluid models")

using  Modia
import ModiaMath


const MediumVariable() = Var(size=())

module BaseProperties
    using  Modia
    using  ..MyModiaMedia
    using  ...TestFluid: MediumVariable

    export BaseProperties_Air, BaseProperties_SimpleMedium

    @model BaseProperties_SimpleMedium begin
        Medium = MediumVariable()  # input
        p      = Float(start=1.0)  # input
        d      = Float(start=1.0)  # output
        der_d  = Float(start=1.0)  # output
    @equations begin
        d      = density(Medium)
        der_d  = 0.0
        end   
    end


    @model BaseProperties_Air begin
        Medium = MediumVariable()               # input
        p      = Float(start=1.0)               # input
        d      = Float(start=1.0, state=false)  # output
        der_d  = Float(start=1.0)               # output
        der_p  = Float(start=1.0)               # output
    @equations begin
        d     = density(Medium, p)
        der_d = der(d)
        der_p = der(p)
        end   
    end
end

using .BaseProperties

function makeDictFromExportedModiaModels!(modul; heading="Keys of model dictionary:")
    dict = Dict{Symbol,Any}()

    println("\n... ", heading)
    for name in names(modul)
        f = getfield(modul, Symbol(name))
        if typeof(f) == Modia.Instantiation.Model
            dict[name] = f
            println("    ", string(name))
        end
    end
    println("")
    return dict
end

const basePropertiesDict = makeDictFromExportedModiaModels!(BaseProperties, heading="Keys of `const baseProperties`:")
getBaseProperties(Medium::AbstractMedium) = basePropertiesDict[Medium.baseProperties]


const water        = MyModiaMedia.Water(100)
const air          = MyModiaMedia.Air(20, 25)
const simpleMedium = MyModiaMedia.SimpleMedium(15.0)


@model Model1 begin
    Medium1 = Var()
    Medium2 = Var()
    d1 = Float()
    d2 = Float()
    @equations begin
        d1 = density(Medium1, 1)
        d2 = density(Medium2, 1)
    end
end

@model Test1 begin
    m = Model1(Medium1=water, Medium2=air)
end


@model Model2 begin
    Medium1 = Var(size=())
    Medium2 = Var(size=())
    d1 = Float(size=())
    d2 = Float(size=())
    @equations begin
        d1 = density(Medium1, 1)
        d2 = density(Medium2, 1)
        Medium1 = Medium2
    end
end

@model Test2 begin
    m = Model2(Medium1=water)
end


@model FluidPort begin
    p      = Float(start=1.0)
    m_flow = Float(start=1.0, flow=true)
    Medium = MediumVariable()
end


@model Source begin
  Medium = AbstractMedium    # Medium MUST be redefined when instanting Source, otherwise error
  P=1
  a=FluidPort()
  d=Float(start=1.0)
  x=Float(start=999.75)
@equations begin
  a.Medium = Medium
  a.p = P
  d = density(a.Medium, a.p)
  der(x) = 1000*(d-x)
  end
end

@model Sink begin
  P=1
  a=FluidPort()
  d=Float()
  x=Float(start=999.75)
@equations begin
  a.p = P
  d = density(a.Medium, a.p)
  der(x) = 1000*(d-x)
  end
end

@model Pipe begin
    R=0.001
    a=FluidPort()
    b=FluidPort()
@equations begin
    b.Medium = a.Medium 
    R*a.m_flow = a.p - b.p
    0 = a.m_flow + b.m_flow
    end
end 

@model PipeSystem begin
  source=Source(P=100.0, Medium=water)
  sink  =Sink(P=10)
  pipe  =Pipe()
@equations begin
  connect(pipe.a, source.a)
  connect(pipe.b, sink.a)
  end
end 

@model ParallelPipeSystem begin
  source=Source(P=100, Medium=water)
  sink=Sink(P=1)
  pipe1=Pipe()
  pipe2=Pipe(R=0.002)
@equations begin
  connect(pipe1.a, source.a)
  connect(pipe1.b, sink.a)

  connect(pipe2.a, source.a)
  connect(pipe2.b, sink.a)
  end
end 


# Test of BaseProperties
@model FluidPort2 begin
    p      = Float(start=1.0)
    Medium = MediumVariable()
end

@model Volume begin
    p_start=1.0
    a = FluidPort2()
    V = Float(start=1.0)
    medium = getBaseProperties(a.Medium)(p=Float(start=p_start) )
@equations begin
    medium.Medium = a.Medium
    medium.p      = a.p
    der(V)        = 2.0
    medium.der_d*V + medium.d*der(V) = 2.0*a.p  # der(m) = der(d*V) = der(d)*V + d*der(V)
    end
end

@model Test_Volume begin
    v = Volume(a=FluidPort2(Medium=air), p_start=10.0)
    t = Float(start=0.0)
@equations begin
    der(t) = 1.0
    end
end


# ----------------------------------------------------

res = simulate(Test1, 1, logTranslation=true)
ModiaMath.plot(res, ("m.d1", "m.d2"), heading="Test1", figure=1)

res = simulate(Test2, 1, logTranslation=true)
ModiaMath.plot(res, ("m.d1", "m.d2"), heading="Test2", figure=2)


res = simulate(PipeSystem, 1, logTranslation=true)
ModiaMath.plot(res, ("sink.d", "source.d"), heading="PipeSystem", figure=3)

res = simulate(ParallelPipeSystem, 1, logTranslation=true)
ModiaMath.plot(res, ("sink.d", "source.d"), heading="ParallelPipeSystem", figure=4)

res = simulate(Test_Volume, 1, logTranslation=true, removeSingularities=false)
ModiaMath.plot(res, ("v.a.p", "v.medium.d", "v.medium.der_d"), heading="Test_Volume", figure=5)

end