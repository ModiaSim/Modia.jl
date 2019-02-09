module TestTearing

println("TestTearing: Tests tearing algorithm of the symbolic handling.")

using Modia


# Desired:
#   using ModiaMath: plot
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia.ModiaMath: plot


# Tearing1 
# Simulates with tearing=false
# No tearing is performed since coefficients of the block equations are not 1 or -1.
@model Tearing1 begin
    x1 = Float(size=())
    x2 = Float(size=())
    x3 = Float(start=1.0)
@equations begin
     2.0*x1 + 3.0*x2 = x3
    -1.4*x1 + 5.0*x2 = 2.0*x3
    der(x3) = -x3
    end
end 
result = simulate(Tearing1, 1.0; logTranslation=true, logSimulation=true, tearing=true, removeSingularities=true)
plot(result, ("x1", "x2", "x3"))

# Tearing1B
# Simulates with tearing=false
# 
@model Tearing1B begin
    x1 = Float(size=())
    x2 = Float(size=())
    x3 = Float(start=1.0)
@equations begin
     2.0*x1 + x2 = x3
    -1.4*x1 + 5.0*x2 = 2.0*x3
    der(x3) = -x3
    end
end 
result = simulate(Tearing1B, 1.0; logTranslation=true, logSimulation=true, tearing=true, removeSingularities=true)
plot(result, ("x1", "x2", "x3"))


# Tearing2
# A correct solution is:
#    x2 := (x3 - 1.0) - 2.0 * (sin)(x1)
#    -1.4*sin(x1) + 5.0*cos(x2) = 2.0*x3
#     der(x3) = -x3
#
@model Tearing2 begin
    x1 = Float(size=())
    x2 = Float(size=())
    x3 = Float(start=1.0)
@equations begin
     2.0*sin(x1) + x2 = x3-1.0
    -1.4*sin(x1) + 5.0*cos(x2) = 2.0*x3
    der(x3) = -x3
    end
end 
result = simulate(Tearing2, 0.3; logTranslation=true, logSimulation=true, tearing=true, removeSingularities=true)
plot(result, ("x1", "x2", "x3"))

@model TearingCombined begin
    x1 = Float(size=())
    x2 = Float(size=())
    x3 = Float(start=1.0)

    x11 = Float(size=())
    x12 = Float(size=())
    x13 = Float(start=1.0)

    x21 = Float(size=())
    x22 = Float(size=())
    x23 = Float(start=1.0)
    
    x31 = Float(size=())
    x32 = Float(size=())
    x33 = Float(start=1.0)
    x34 = Float(start=1.0)
    x35 = Float(size=())
@equations begin
     2.0*x1 + 3.0*x2 = x3
    -1.4*x1 + 5.0*x2 = 2.0*x3
    der(x3) = -x3
    
     2.0*x11 + x12 = x13
    -1.4*x11 + 5.0*x12 = 2.0*x13
    der(x13) = -x13

     2.0*sin(x21) + x22 = x23-1.0
    -1.4*sin(x21) + 5.0*cos(x22) = 2.0*x23
    der(x23) = -x23    
    
     2.0*sin(x31) + x32 = x33-1.0
    -1.4*sin(x35) + 5.0*cos(x34) = 2.0*x33
    x34 = 1.0*x32
    2*x35 = 2*x31
    der(x33) = -x33    
    end
end 
result = simulate(TearingCombined, 1.0; logTranslation=true, logSimulation=true, tearing=true, removeSingularities=true)
plot(result, ("x1", "x2", "x3"))


end

