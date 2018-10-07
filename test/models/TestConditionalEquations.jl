module ConditionalEquations

println("\nDemonstrating conditional equations")

using Modia

# Desired:
#   using ModiaMath: plot
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia.ModiaMath: plot

@model Conditional begin
    x = Var(start=1.0)
    u = Var()
    steadyState = Boolean(false, variability=parameter)
    cond = false
    y = Var()
@equations begin
    u = sin(time)
    if ! steadyState
      der(x) + 2x = u  
    else
      0 + 2x = u
    end
    if ! cond
      y = 1
    else
      y = 2
    end
  end
end 

result = simulate(Conditional, 2, logTranslation=true, removeSingularities=false)
plot(result, ("x", "y"), heading="Conditional", figure=1)


@model ConditionalInstance1 begin
    @extends Conditional(steadyState=true, cond=true)
end

result = simulate(ConditionalInstance1, 2, logTranslation=true, removeSingularities=false)
plot(result, ("x", "y"), heading="ConditionalInstance", figure=2)


@model ConditionalInstance2 begin
    @extends Conditional(steadyState=false, cond=false)
end

result = simulate(ConditionalInstance2, 2, logTranslation=true, removeSingularities=false)
plot(result, ("x", "y"), heading="ConditionalInstance", figure=3)


@model Conditional2 begin
    x = Var(start=1.0)
    u = Var()
    steadyState = Boolean(false, variability=parameter)
    cond = false
    y = if cond; Var() end
@equations begin
    u = sin(time)
    if ! steadyState #&& cond
      der(x) + 2x = u  
    else
      0 + 2x = u
    end
    if cond
      y = 1
    end
  end
end 

result = simulate(Conditional2, 2, logTranslation=false, removeSingularities=false)
plot(result, ("x"), heading="Conditional2", figure=1)

end
