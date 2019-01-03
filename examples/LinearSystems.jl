module LinearSystems

println("\nLinearSystems: Demonstrates type and size deduction.")

using Modia
using Modia.Blocks: ABCD, Switch


# Desired:
#   using ModiaMath: plot
#
# In order that ModiaMath need not to be defined in the user environment, it is included via Modia:
using Modia.ModiaMath: plot


@model MySISOABCD begin
  @extends ABCD(x=Float(start=[0.0; 0.0]), A=[-1 -2; 0.5 -2], B=[1; 1], C=[1 1], D=[0], u=1.0)
  @inherits u
@equations begin
  # u = 1.0 # Does not work since size deduction finds the wrong "solution" which is not consistent (inv finds minimum norm solution).
  end
end

result = checkSimulation(MySISOABCD, 10, "x", 0.5000001198147023, logTranslation=true) #logOnFile = true)
# plot(result, "x", heading="MySISOABCD", figure=13)

@model MyMIMOABCD begin
  s=Switch(u1=[1; 1], u2=[2; 3])
  lambda=0.5
  omega=10
  @extends ABCD(x=Float(start=[0.0; 0.0]), A=[-2*lambda omega^2; -1 0], B=[1 0; 0 1], C=[1 0; 0 1], D=zeros(2,2))
  @inherits u
@equations begin
  s.sw = time>5
  u = s.y
  end
end

result = checkSimulation(MyMIMOABCD, 10, "x", 0.004778288940817275, storeEliminated=false) # storeEliminated=false needed. Investigate
plot(result, "x", heading="MyMIMOABCD with generic switch", figure=14)

end