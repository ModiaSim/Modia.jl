module TestSpatialDiscretization

using Modia

# Desired:
#   using ModiaMath: plot
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia.ModiaMath: plot


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
plot(result, "u", figure=1)



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
plot(result, "ueven", figure=2)


n = 100
@model SpatialDiscretization3 begin
    u = Var(size=(2*n+1,))
    ueven = Var(start=[i for i in 1:n])
    uodd = Var(size=(n+1,))
    @equations begin
        u = splice(uodd, ueven)
        der(ueven) = - u[2*(1:n)] + u[2(1:n)-1]
        uodd = ones(n+1)*time
    end
end

@static if VERSION < v"0.7.0-DEV.2005"
    @time result = simulate(SpatialDiscretization3, 1, logTranslation=true, storeEliminated=false)
    plot(result, "ueven", figure=3)
end

@model SpatialDiscretization4 begin
    u = Var(size=(2*n+1,))
    ueven = Var(start=[i for i in 1:n])
    uodd = Var(size=(n+1,))
    @equations begin
        u = splice(uodd, ueven)
        der(ueven) = ueven
        uodd = zeros(n+1)
    end
end

@time result = simulate(SpatialDiscretization4, 1, logTranslation=true, storeEliminated=false)
plot(result, "ueven", figure=4)


maskEven(n) = [if isodd(i); 0 else 1 end for i in 1:n]

@model SpatialDiscretization5 begin
    u = Float(start=[if iseven(i); div(i,2) else 0 end for i in 1:2*n+1], state=false)
    @equations begin
        broadcast(*, maskEven(2*n+1), der(u)) = u 
    end
end

# @time result = simulate(SpatialDiscretization5, 1, logTranslation=true, storeEliminated=false, removeSingularities=false)
# plot(result, "u", figure=5)

end
