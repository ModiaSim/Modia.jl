module TestSpatialDiscretization

using Modia
using ModiaMath:plot

const n = 5
@model SpatialDiscretization begin
    u = Var(start=zeros(n))
    c = 1
    @equations begin
        # der(u[2:n]) = ( u[2:n] - u[1:n-1] ) / c
        der(u) = [1; ( u[2:n] - u[1:n - 1] ) / c]
    end
end

result = simulate(SpatialDiscretization, 1)
plot(result, "u")


splice(uodd, ueven) = [if isodd(i); uodd[div(i,2)+1] else ueven[div(i,2)] end for i in 1:length(ueven)+length(uodd)] 

@model SpatialDiscretization2 begin
    u = Var(size=(5,))
    ueven = Var(start=[1,2])
    uodd = Var(size=(3,))
    @equations begin
      u = splice(uodd, ueven)
      der(ueven) = - u[[2,4]] + u[[1,3]]
      uodd = ones(3)*time
    end
end

result = simulate(SpatialDiscretization2, 1, logTranslation=true, storeEliminated=false)
plot(result, "ueven")

end
