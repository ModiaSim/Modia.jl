module HeatTransfer2D

using Modia

# Desired:
#   using ModiaMath: plot
#
# In order that ModiaMath need not to be defined in the user environment, it is included via Modia:
using Modia.ModiaMath: plot


# Definition of model  
const N      = 30 # 100 # Number of nodes in x-direction (same number of nodes in y-direction
const L      = 0.2   # Length in x-, and y-direction
const T0     = 290.0 # Initial temperature
const TNx    = 310.0 # Temperature along the y-line outside x-node N
const TNy    = 330.0 # Temperature along the x-line outside y-node N
const cp     = 910.0 # Material Heat Capacity
const lambda = 237.0 # Material thermal conductivity
const rho    = 271.0 # Material Density
const dL     = L / N # Element length in x- and y-direction
const c      = lambda/(cp*rho*dL*dL)
const N2     = N*N
   
geti(i,j)  = (j-1)*N + i

""" der(T) = heatTransfer2D(T)"""  
function heatTransfer2D(T::Matrix{Float64}) 
   der_T::Matrix{Float64} = zeros(N,N)
   for i in 1:N
      for j in 1:N
         qx1 = if i > 1 ; T[i - 1, j] - T[i,j] else 0.0              end
         qx2 = if i < N ; T[i + 1, j] - T[i,j] else 2*(TNx - T[i,j]) end
         qy1 = if j > 1 ; T[i, j - 1] - T[i,j] else 0.0              end
         qy2 = if j < N ; T[i, j + 1] - T[i,j] else 2*(TNy - T[i,j]) end

         der_T[i,j] = c*(qx1 + qx2 + qy1 + qy2)
      end
   end
   return der_T
end  

# Not used presently:
""" Structure of der( heatTransfer2D, T ) """
function jacobian_incidence(::typeof(heatTransfer2D), args...)
   I::Vector{Int} = fill(1,5*N2)
   J::Vector{Int} = fill(1,5*N2)
   
   r = 0
   for i in 1:N
      for j in 1:N
         i0 = geti(i,j)
         i1 = geti(max(i-1,1),j)
         i2 = geti(min(i+1,N),j)
         i3 = geti(i,max(j-1,1))
         i4 = geti(i,min(j+1,N))

         r += 1
         I[r] = i0
         J[r] = i0
         
         r+= 1
         I[r] = i0
         J[r] = i1

         r+= 1
         I[r] = i0
         J[r] = i2

         r+= 1
         I[r] = i0
         J[r] = i3
         
         r+= 1
         I[r] = i0
         J[r] = i4
      end
   end

   return sparse(I,J,1)
end

  
@model HeatTransfer begin
   T = Float(start=fill(T0,N,N), fixed=true) # Temperature of the nodes
@equations begin
   der(T) = heatTransfer2D(T)
  end
end 

result = simulate(HeatTransfer, 30)
#plot(result, "T", figure=22)

res = Dict{AbstractString,Any}()
res["time"] = result["time"]
T = result["T"]
T = reshape(T, size(T,1), N, N)
v1 = "T[$(div(N,2)),1]"
v2 = "T[1, $(div(N,2))]"
v3 = "T[$N,1]"
v4 = "T[1, $N]"
res[v1] = T[:, div(N,2),1]
res[v2] = T[:, 1, div(N,2)]
res[v3] = T[:, N, 1]
res[v4] = T[:, 1, N]

plot(res, (v1, v2, v3, v4), figure=20, heading="HeatTransfer2D: N=$N, temperature [K]")

end
