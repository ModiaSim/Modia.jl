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
# Translation fails with tearing=true: x2 not defined
# The code in the log is:
#    x1 := -1.4 \ (2.0x3 - 5.0x2)
#    x2 := 3.0 \ (x3 - 2.0x1)
# and this is not good, because tearing introduces a division that might result in a division by zero.
# Tearing should only solve for one of the variables, x1 or x2, and keep one residue equation
#
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
result = simulate(Tearing1, 1.0; logTranslation=true, logSimulation=true, tearing=true, removeSingularities=false)
plot(result, ("x1", "x2", "x3"))


# Tearing2
# The code in the log is:
#   (sin)(x1) = -1.4 \ (2.0x3 - 5.0 * (cos)(x2))
#   x2 := (x3 - 1.0) - 2.0 * (sin)(x1)
#   der(x3) = -x3
# and this is not good, because tearing introduces a potential division by zero.
# The correct solution would be:
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
result = simulate(Tearing2, 0.3; logTranslation=true, logSimulation=true, tearing=true, removeSingularities=false)
plot(result, ("x1", "x2", "x3"))


end
