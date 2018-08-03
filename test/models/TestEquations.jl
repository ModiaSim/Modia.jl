module TestEquations

using Modia
using ModiaMath.plot
using Base.Test

# @testset "Equations" begin

@model Test begin
  a = -1; b = 1; c = 1
  n=10; 
  M=Float()
  u=Float(); x=Float(start=0.0)
  on=Variable(); s=Variable()
  QR=Float(); Q=Float(); R=Float()
  dummy=Variable(size=())
@equations begin
  u = if s == "on"; 10 else M end
  der(x) - a*x = b*(s == "on" ? u : 10*u[1,1])
  on = time>1 && time<3
  s := on ? "on" : "off"
  dummy = if time < 0.00001 || abs(time-2.0) < 0.01; println("Size of u: ", size(u), " time=", time) else nothing end
  ((9+1)*diagm(1:n)+ones(n,n))^3*M*diagm(1:n)=eye(n)
  QR = qr(M)
  Q=QR[1]; R=QR[2]
  end
end;

result = simulate(Test, 5, storeEliminated=false, removeSingularities=false)
plot(result, ("x"), figure=1, heading="TestEquations", figure=2)
@test result["x"][end] == 1.1776785083961023

end