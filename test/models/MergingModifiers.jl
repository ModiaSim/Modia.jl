module MergingModifiers

println("\nDemonstrating merging modifiers")

using Modia

# Desired:
#   using ModiaMath: plot
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia.ModiaMath: plot


@model M begin
  x = Float(start=1.0, info="state")
@equations begin
    der(x) + 2x = 0  
  end
end 

result = simulate(M, 1, logTranslation=true, removeSingularities=false)
plot(result, ("x"), heading="M", figure=1)

# ------------------------------

@model MInstance begin
  @extends M(x=Var(start=2.0))  # info forgotten
end

result = simulate(MInstance, 1, logTranslation=true, removeSingularities=false)
plot(result, ("x"), heading="MInstance", figure=2)

# ------------------------------


@model MInstance2 begin
  @extends M(y=Var(start=2.0)) # Should give error message for missing y
end

result = simulate(MInstance2, 1, logTranslation=true, removeSingularities=false)
plot(result, ("x"), heading="MInstance2", figure=3)

# ------------------------------

#=
merge(c, e, v) = begin @show c; dump(c); @show e; dump(e.initializers[1].name); dump(e.initializers[1].fdef.args[2].args[2].args[2]); @show v return c end

@model MInstance3 begin
# @extends M(x(start=2.0))  # UndefVarError: x not defined
# @extends M(x=Var(start=2.0))
  @extends M(x=merge(Var(start=2.0), M, :x))
end

result = simulate(MInstance3, 1, logTranslation=true, removeSingularities=false)
plot(result, ("x"), heading="MInstance3", figure=4)
=#

end
