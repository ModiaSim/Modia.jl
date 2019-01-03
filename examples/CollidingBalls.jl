module TestBalls

println("\nCollidingBalls: Demonstrating the use of allInstances to set up contact force between any number of balls")

using Modia

# Desired:
#   using ModiaMath: plot
#   using Test
#   using LinearAlgebra
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia.ModiaMath: plot

@static if VERSION < v"0.7.0-DEV.2005"
  using Base.Test
else
  using Modia.Test
  using Modia.LinearAlgebra
end

const k=10000
const d=0
const radius=0.05

springDamper(r, v) = k*r + d*v

function getForce(r, v, positions, 
  velocities, contactLaw)
  # Calculate the forces on this body
  force = zeros(2)
  for i in 1:length(positions)
    pos = positions[i]
    vel = velocities[i]
    if r != pos
      delta = r - pos
      deltaV = v - vel
      f = if norm(delta) < 2*radius; 
        -contactLaw((norm(delta)-
          2*radius)*delta/norm(delta), deltaV) else 
        zeros(2) end
      force += f
    end
  end
  return force
end

@model Ball begin
  r = Var()
  v = Var()
  f = Var()
  m = 1.0
  event = Boolean(size=(), start=false) 
@equations begin
  der(r) = v
  m*der(v) = f
  f = getForce(r, v, allInstances(r),
#    allInstances(v), (r,v) -> (k*r + d*v))
    allInstances(v), springDamper)
    event = Modia.Synchronous.positive(f[1])
  end
end


@model Balls1 begin
  b1 = Ball(r = Var(start=[0.0,2]), v = Var(start=[1,0]))
  b2 = Ball(r = Var(start=[0.5,2]), v = Var(start=[-1,0]))
  b3 = Ball(r = Var(start=[1.0,2]), v = Var(start=[0,0]))
end

#=
result = checkSimulation(Balls1, 1.5, "b3.r", 2.0, expandArrayIncidence=true, storeEliminated=false)
plot(result, ("b1.r", "b2.r", "b3.r"), figure=1)
=#

# -----------------------------

@model Balls2 begin
  b1 = Ball(r = Var(start=[0.0,2]), v = Var(start=[1,0]))
  b2 = Ball(r = Var(start=[0.5,2]), v = Var(start=[0,0]))
  b3 = Ball(r = Var(start=[1.0,2]), v = Var(start=[-1,0]))
end

#=
result = checkSimulation(Balls2, 1.5, "b3.r", 2.0, expandArrayIncidence=true, storeEliminated=false)
plot(result, ("b1.r", "b2.r", "b3.r"), figure=2)
=#

# -----------------------------

@model Balls3 begin
  b1 = Ball(r = Var(start=[1.0,0]), v = Var(start=[1,0]))
  b2 = Ball(r = Var(start=[2.0,0]), v = Var(start=[0,0]))
  b3 = Ball(r = Var(start=[3.0,0]), v = Var(start=[0,0]))
  b4 = Ball(r = Var(start=[4.0,0]), v = Var(start=[0,0]))
  b5 = Ball(r = Var(start=[5.0,0]), v = Var(start=[0,0]))
  b6 = Ball(r = Var(start=[10.0,0]), v = Var(start=[-1,0]))
end

result = simulate(Balls3, 8, expandArrayIncidence=true, storeEliminated=false)
plot(result, Tuple(["b$i.r" for i in 1:5]), heading="Colliding balls", grid=false, figure=19)

# -----------------------------

#=
nBalls=4
@model Balls begin
#  balls = [Ball(r = Var(start=[i,2]), v = Var(start=if i==1; [1,0] else [0,0] end)) for i in 1:nBalls]
  balls = [Ball() for i in 1:nBalls]
end

result = simulate(Balls, 1.5, expandArrayIncidence=true, storeEliminated=false, logTranslation=true, logOnFile=false)
=#

end