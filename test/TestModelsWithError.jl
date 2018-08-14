module TestModelsWithError

using Modia
using Modia.Electric

println("MODELS WITH ERRORS")
println()

Modia.ModiaLogging.resetTestStatus()

@model WrongName begin
    R = Resist()
end 
checkSimulation(WrongName, 1)


@model WrongParameterName begin
    R = Resistor(res=100.0)
end 
# Pantelides loops!!!
# checkSimulation(WrongParameterName, 1)


@model WrongConnect1 begin
    R = Resistor(R=100.0)
    C = Capacitor()
    @equations begin
        connect(R.n, C.pin)
    end
end
checkSimulation(WrongConnect1, 1)


@model WrongConnect2 begin
    R = Resistor(R=100.0)
    C = Capacitor()
    @equations begin
        connect(R.n)
    end
end
checkSimulation(WrongConnect2, 1)


@model WrongEquation1 begin
    R = Resistor(res=100.0)
    y = Var()
    @equations begin
        R.I = y
    end
end 
println("Should complain about I not declared!")
checkSimulation(WrongEquation1, 1)

@model WrongEquation2 begin
    R = Resistor(res=100.0)
    @equations begin
        R.i = 0
    end
end 
checkSimulation(WrongEquation2, 1)

@model Base begin
    v = Var(start=0.0)
end

@model Extended1 begin
    @extends MyBasic # Forgot()
end
checkSimulation(Extended1, 1)


@model Extended2 begin
    @extends Base # Forgot ()
end
checkSimulation(Extended2, 1)

println("UndefinedVariable test disabled since it stops script. See model code.")
#=
@model UndefinedVariable begin
    @equations begin
        v = 0.0 # ERROR: LoadError: UndefVarError: v not defined 
    end
end
checkSimulation(UndefinedVariable, 1)
=#

println("Extended3 test disabled since it stops script. See model code.")
@model Extended3 begin
    @extends Base()
    @equations begin
        # v = 0.0 # ERROR: LoadError: UndefVarError: v not defined 
        this.v = 0.0
    end
end
checkSimulation(Extended3, 1)

println("Does not complain about undefined dummy!")
@model Extended4 begin
    @extends Base()
    @inherits v, dummy
end
checkSimulation(Extended4, 1)

@model Extended5 begin
    @extends Base()
    @inherits v, dummy
    w = Var()
end
checkSimulation(Extended5, 1)

Modia.ModiaLogging.printTestStatus()

end
