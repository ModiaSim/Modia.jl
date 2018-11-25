"""
    module ExactlyRemoveSingularities - exactly remove singularities from A*v = 0 with A an Integer matrix.
    
Use this module in the following way:

   using ExactlyRemoveSingularities
   result = removeSingularities(A,ix,iy)    # result = (iya, eqr, ix1, ix2, eqx, A1, A2)
   printRemoveSingularities(A,result,names)   

where A is a dense or sparse Int matrix A[:,:] describing a linear system

    A*v = 0

with

  - v[ix]: Variables used in the derivative operator, der(v[ix]), so these are the potential states.

  - v[iy]: Variables that must be solved from ```A*v = 0```, because these variables do not appear in the remaining part of the model equations (e.g. absolute potentials). If variables v[iy] cannot be uniquely solved from ```A*v = 0```, then in argument iya the variables are listed that must be set arbitrarily so that v[iy] can be uniquely solved from the equations. For example, if v[iy] contains the absolute potentials and no ground object is present, then v[iya] are the absolute potentials that must be explicitly set, e.g., by introducing an equation v[iya]=0, in order that the equations have a unique solution. Additionally, a warning message should be printed that v[iya] was set to 0.)

The function return arguments are all of type Int and have the following meaning:

  - iya: variables v[iya] can be arbitrarily set (iya is a sub-set of iy)

  - eqr: equations ```A[eqr,:]*v = 0``` are redundant and should be removed (they do not change the result)

  - ix1: x[ix1] are the dependent states (ix1 is a sub-set of ix)

  - ix2: i[ix2] are the independent states (ix2 is a sub-set of ix)

  - eqx, A1, A2: equations ```A[eqx,:]*v = 0``` should be replaced by the equations ```A1*x1 + A2*x2 = 0```, where A1 is square, upper triangular and has full rank with A1[i,i] != 0. In other words, x1 can be computed explicitely from these equations given x2.
    
The calling program should do the following:

   - Remove equations ```A[eqr,:]*v = 0```, because these are redundant equations.
   - Replace equations ```A[eqx,:]*v = 0``` by ```A1*x1 + A2*x2 = 0```, in order that index reduction can be performed by a structural algorithm (e.g. Pantelides algorithm).
   - Introduce equations v[iya] = 0 and print a warning message that these variables have been arbitrarily set to zero. This is needed, in order that the equations have a unique solution.
"""
module ExactlyRemoveSingularities

export removeSingularities, upperTrapezoidal
export printobj, printRemoveSingularities

@static if !(VERSION < v"0.7.0-DEV.2005")
    using SparseArrays
end

"""
    (Afac, rk, p1, p2, Bfac) = upperTrapezoidal(A,B) 
    (Afac, rk, p1, p2)       = upperTrapezoidal(A)
    
A and B are **Int** matrices with dimensions A[n1,n2], B[n1,:] or B[n1], so B can be a matrix or a vector. Furthermore A and B can be dense or sparse matrices (SparseMatrixCSC). If B is not present, it is interpreted
as zero matrix. Typically A has more columns as rows (n2 > n1). The function transforms the equation

    A*X = B

to

    Afac*X[p2,:] = Bfac

with

- Afac an **upper trapezoidal** Int matrix with Afac[i,i] != 0 for i<=rk; Afac[rk+1:end,:] = 0

- Bfac an Int matrix or vector (not present if B is not present)

- rk the rank of A and of Afac

- p1 the row permutation vector of A

- p2 the column permutation vector of A

The transformation is performed by a fraction-free Gauss elimination with
exact (integer) computation. It is then possible to solve for X by splitting
X[p2,:] = [X1;X2] with X1[rk,rk], providing any values for X2 and computing X1 from 

    Afac[1:rk,1:rk]*X1 = Bfac[1:rk,:] - Afac[1:rk,rk+1:end]*X2

Since Afac[1:rk,1:rk] is **upper triangular** with non-zeros on the diagonal,
this is just a forward recursion. Note, the (full) pivoting is performed
in such a way that it is tried to utilize the **upper** part of X in X1.
Therefore, variables that must be solved for should be in the upper part of X.

The permutation vector p1 can be used to determine the equations in A and B that
can be removed without changing the result (provided Bfac[rk+1:end,:] = 0):

    p1red = p1[rg+1:end]
    A[p1red,:]*X = B[p1red,:] are redundant equations that can be removed

This implementation is based on:

  - Bareiss E.H. (1968): Sylvester’s identity and multistep integer-preserving Gaussian elimination. Mathematics of Computation, Vol. 22, No. 103, pp. 565 – 578. Download: http://www.ams.org/journals/mcom/1968-22-103/S0025-5718-1968-0226829-0/

  - Turner P.R. (1995): A simplified fraction-free Integer elimination algorithm. Naval Air Warfare Center. Report NAWCADPAX--96-196-TR. Download: http://www.dtic.mil/cgi-bin/GetTRDoc?AD=ADA313755

  - Dureisseix D. (2012): Generalized fraction-free LU factorization for singular systems with kernel extraction. Linear Algebra and its Applications, Volume 436, Issue 1, pp. 27-40. Download: http://www.sciencedirect.com/science/article/pii/S0024379511004617#
"""
function upperTrapezoidal(A::Matrix{Int}, B::Matrix{Int})
    @assert(size(A, 1) == size(B, 1))
    (n1, n2) = size(A)
    nb   = size(B, 2)
    p1   = collect(1:n1)
    p2   = collect(1:n2)
    Afac = copy(A)
    Bfac = copy(B)
    rk   = n1
    p1k  = 0
    p2k  = 0
    Aik  = 0
    temp = 0
    tempRow    = fill(0, 1, nb)
    pivot      = 1   
    oldPivot   = 1
    pivotFound = false

   
    for k = 1:n1
        # Search for a pivot in Afac[k:n1,min(k,n2):n2]
        if k <= n2
            if Afac[k,k] == 0
                pivotFound = false
                for k2 = k:n2
                    for k1 = k:n1
                        if Afac[k1,k2] != 0 
                            pivotFound = true
                            p1k = k1
                            p2k = k2
                            break
                        end
                    end
                    pivotFound && break
                end     
            
                if pivotFound
                    # Row interchange
                    if p1k != k
                        for i = p2k:n2
                            temp = Afac[k,i]
                            Afac[k  ,i] = Afac[p1k,i]
                            Afac[p1k,i] = temp
                        end
                        tempRow = Bfac[k,:]
                        Bfac[k,:]   = Bfac[p1k,:]
                        Bfac[p1k,:] = tempRow
      
                        temp    = p1[k]
                        p1[k]   = p1[p1k]
                        p1[p1k] = temp
                    end
                  
                    # Column interchange
                    if p2k != k
                        for i = 1:n1
                            temp = Afac[i,k]
                            Afac[i,k  ] = Afac[i,p2k]
                            Afac[i,p2k] = temp
                        end                   
                        temp    = p2[k]
                        p2[k]   = p2[p2k]
                        p2[p2k] = temp            
                    end
                else
                    # no pivot found: submatrix Afac[k:n1,:] has only zeros
                    rk = k - 1           
                    return (Afac, rk, p1, p2, Bfac)
                end
            end
        else
            # submatrix Afac[k:n1,:] has only zeros
            rk = k - 1           
            return (Afac, rk, p1, p2, Bfac)  
        end            

      # Afac[k,k] != 0
        pivot = Afac[k,k]          
        for i = k + 1:n1
            Aik  = Afac[i,k]
            Bfac[i,:] = div.(pivot * Bfac[i,:] - Aik * Bfac[k,:], oldPivot)  # guaranteed no remainder        
            for j = k + 1:n2
                Afac[i,j] = div(pivot * Afac[i,j] - Aik * Afac[k,j], oldPivot)  # guaranteed no remainder
            end
            Afac[i,k] = 0
        end
        oldPivot = pivot
    end
    return (Afac, rk, p1, p2, Bfac)
end

toVector(result) = (result[1], result[2], result[3], result[4], vec(result[5]))
upperTrapezoidal(A::Matrix{Int}, b::Vector{Int}) = toVector(upperTrapezoidal(A, reshape(b, size(b, 1), 1))) 
upperTrapezoidal(A::Matrix{Int})                 =  upperTrapezoidal(A, fill(0, size(A, 1), 0))

function upperTrapezoidal(A::SparseMatrixCSC{Int,Int}, B::SparseMatrixCSC{Int,Int})
    @assert(size(A, 1) == size(B, 1))
    (n1, n2) = size(A)
    nb   = size(B, 2)
    p1   = collect(1:n1)
    p2   = collect(1:n2)
    Afac = copy(A)
    Bfac = copy(B)

    rk   = n1
    p1k  = 0
    p2k  = 0
    Aik  = 0
    temp = 0
    pivot      = 1   
    oldPivot   = 1
    pivotFound = false
   
    for k = 1:n1
        # Search for a pivot in Afac[k:n1,min(k,n2):n2] columnwise
        if k <= n2
            if Afac[k,k] == 0      
                Arows = rowvals(Afac)
                pivotFound = false         
                Avals = nonzeros(Afac)
                for k2 = k:n2  
                    for kk1 in nzrange(Afac, k2)
                        k1  = Arows[kk1]
                        val = Avals[kk1]
                        if k1 >= k && val != 0
                            pivotFound = true
                            p1k = k1
                            p2k = k2
                            break
                        end
                    end
                    pivotFound && break
                end
            
                if pivotFound
                    # Row interchange
                    if p1k != k
                        ii = p2k:n2
                        tempAr       = Afac[k,ii]
                        Afac[k  ,ii] = Afac[p1k,ii]
                        Afac[p1k,ii] = tempAr
                  
                        if nb > 0
                            tempRow = Bfac[k,:]
                            Bfac[k,:]   = Bfac[p1k,:]
                            Bfac[p1k,:] = tempRow
                        end
      
                        temp    = p1[k]
                        p1[k]   = p1[p1k]
                        p1[p1k] = temp
                    end
      
                    # Column interchange
                    if p2k != k
                        tempAc      = Afac[:,k]
                        Afac[:,k  ] = Afac[:,p2k]
                        Afac[:,p2k] = tempAc
      
                        temp    = p2[k]
                        p2[k]   = p2[p2k]
                        p2[p2k] = temp            
                    end
                else
                    # no pivot found: submatrix Afac[k:n1,:] has only zeros
                    rk = k - 1           
                    return (Afac, rk, p1, p2, Bfac)
                end
            end
        else
            # no pivot found: submatrix Afac[k:n1,:] has only zeros
            rk = k - 1           
            return (Afac, rk, p1, p2, Bfac)      
        end
      
        # Subtract row k from row i = k+1:n1       
        j = k + 1:n2           
        for ii in nzrange(Afac, k)    # [i,k], i = k:n1
            Arows = rowvals(Afac)
            Avals = nonzeros(Afac)      
            i = Arows[ii]
            if i == k
            # Afac[k,k] != 0    
                pivot = Avals[ii]     
            elseif i > k 
                Aik = Avals[ii]
                Avals[ii] = 0           
                if nb > 0 
                    Bfac[i,:] = div.(pivot * Bfac[i,:] - Aik * Bfac[k,:], oldPivot)  # guaranteed no remainder        
                end
                Afac[i,j] = div.(pivot * Afac[i,j] - Aik * Afac[k,j], oldPivot)  # guaranteed no remainder       
            end
        end
        oldPivot = pivot
    end
    return (Afac, rk, p1, p2, Bfac)
end

upperTrapezoidal(A::SparseMatrixCSC{Int,Int}) =  upperTrapezoidal(A, sparse(fill(0, size(A, 1), 0)))

"""-----------------------------------------------------------------------
    (iya, eqr, ix1, ix2, eqx, A1, A2) = removeSingularities(A, ix, iy)

A is a dense or sparse Int matrix A[:,:] describing a linear system

    A*v = 0

where

  - v[ix]: Variables used in the derivative operator, der(v[ix]), so these are the potential states.

  - v[iy]: Variables that must be solved from ```A*v = 0```, because these variables do not appear in the remaining part of the model equations (e.g. absolute potentials). If variables v[iy] cannot be uniquely solved from ```A*v = 0```, then in argument iya the variables are listed that must be set arbitrarily so that v[iy] can be uniquely solved from the equations. For example, if v[iy] contains the absolute potentials and no ground object is present, then v[iya] are the absolute potentials that must be explicitly set, e.g., by introducing an equation v[iya]=0, in order that the equations have a unique solution. Additionally, a warning message should be printed that v[iya] was set to 0.)

The function return arguments are all of type Int and have the following meaning:

  - iya: variables v[iya] can be arbitrarily set (iya is a sub-set of iy)

  - eqr: equations ```A[eqr,:]*v = 0``` are redundant and should be removed (they do not change the result)

  - ix1: x[ix1] are the dependent states (ix1 is a sub-set of ix)

  - ix2: i[ix2] are the independent states (ix2 is a sub-set of ix)

  - eqx, A1, A2: equations ```A[eqx,:]*v = 0``` should be replaced by the equations ```A1*x1 + A2*x2 = 0```, where A1 is square, upper triangular and has full rank with A1[i,i] != 0. In other words, x1 can be computed explicitely from these equations given x2.
    
The calling program should do the following:

   - Remove equations ```A[eqr,:]*v = 0```, because these are redundant equations.
   - Replace equations ```A[eqx,:]*v = 0``` by ```A1*x1 + A2*x2 = 0```, in order that index reduction can be performed by a structural algorithm (e.g. Pantelides algorithm).
   - Introduce equations v[iya] = 0 and print a warning message that these variables have been arbitrarily set to zero. This is needed, in order that the equations have a unique solution.
"""
function removeSingularities(A::T, ix::Vector{Int}, iy::Vector{Int}) where T <: AbstractMatrix{Int}
    # Split equation system A*v = 0 into V*[vy,vr] = X*(-x), where
    #   x  = v[ix]
    #   vy = v[iy]
    #   vr = the remaining variables
    (ne, nv) = size(A)   # number of equations (ne) and variables (nv)
    nx = length(ix)     # number of potential states
    set_iy = Set(iy)
    ivr = collect(setdiff(setdiff(Set(1:nv), Set(ix)), set_iy))  # indices of vr
    ivall = [iy;ivr]
    V = [ A[:,iy] A[:,ivr] ]
    X = A[:,ix];

    # Transform V to upper trapezoidal form
    #    Vfac*v[pv2] = Xfac*(-x)
    (Vfac, rk, pv1, pv2, Xfac) = upperTrapezoidal(V, X)

    # iya are the independent variables that are in set_iy
    ivind = ivall[ pv2[rk + 1:end] ]
    iya   = collect(intersect(Set(ivind), set_iy))
  
    if nx == 0
        eqr = pv1[rk + 1:ne]
        eqx = fill(0, 0)
        ix1 = fill(0, 0)
        ix2 = fill(0, 0)
        A1  = fill(0, 0, 0)
        A2  = fill(0, 0, 0)
   
    else
        # Transform state constraints to upper trapezoidal form
        #    Xc*x = 0     // Xc = Xfac[rk+1:end,:]
        #
        #    [Xcfac1 XcFac2;
        #       0       0  ] * [x1,x2] = 0
        Xc = Xfac[rk + 1:end,:]
        (Xcfac, rk2, px1, px2) = upperTrapezoidal(Xc)
        eqx2 = px1[1:rk2]
        eqr2 = px1[rk2 + 1:end]
     
        eqx  = pv1[rk .+ eqx2]
        eqr  = pv1[rk .+ eqr2]
        ix1  = ix[ px2[1:rk2] ]
        ix2  = ix[ px2[rk2 + 1:end] ]
        A1   = Xcfac[1:rk2, 1:rk2]
        A2   = Xcfac[1:rk2, rk2 + 1:end]
    end
  
    return (iya, eqr, ix1, ix2, eqx, A1, A2)
end


"""---------------------------------------------------------------
    str = displayAsString(x)

Return the result of the Julia built-in command display(x) as string 
instead of displaying it on the terminal.
"""
displayAsString(x) = sprint((io, x) -> display(TextDisplay(io), x), x)

@static if VERSION < v"0.7.0-DEV.2005"
else
    searchindex(s,t) = first(something(findfirst(t, s), 0:-1))
    replace(s::String, t::String, u::String) = Base.replace(s, t=>u)
end

"""---------------------------------------------------------------
    str = displayAsString2(x)

Return the result of the Julia built-in command display(x) as string 
instead of displaying it on the terminal, remove the preceding type
information and shift the whole output to the right
"""
function displayAsString2(obj)
    str  = displayAsString(obj)
    str2 = str[searchindex(str, "\n") + 1:end]
    str3 = "   " * replace(str2, "\n", "\n   ")
end


"""---------------------------------------------------------------
    printobj(name,obj)

Pretty print an instance obj as "<name> = <display(obj)>". 
"""
function printobj(name, obj)
    str = displayAsString2(obj)
    println("\n$name =\n$str")
end

function printEquation(coeff, names::Vector{String})
    @assert(length(coeff) == length(names))
    print("   ")
    first = true
    for i = 1:length(names)
        if coeff[i] != 0
            if first
                first = false
                if coeff[i] < 0 
                    print("-")
                end
            else
                if coeff[i] < 0 
                    print(" - ")
                elseif coeff[i] > 0
                    print(" + ")
                end
            end
         
            c = abs(coeff[i])
            if c != 1
                print(c, "*")
            end
            print(names[i])
        end
    end
    println(" = 0")
end

   
"""
    printRemoveSingularities(A, r, names)
    
Prints the result of r=removeSingularities(..), where

  - A: matrix of removeSingularities(..)
  
  - r: the tuple returned by removeSingularities. 
  
  - names: string vector containing the names of the variables v from A*v = 0.
"""
function printRemoveSingularities(A, r, names)
    (iya, eqr, ix1, ix2, eqx, A1, A2) = r
   
    if length(iya) > 0
        iya_names = displayAsString2(names[iya])
        println("\n  Variables that can be arbitrarily set:\n", iya_names)  
    end
   
    if length(eqr) > 0
        println("\n  Equations that can be removed, since redundant: ", eqr)
        for i = 1:length(eqr)
            printEquation(A[eqr[i],:], names)
        end
    end
   
    if length(ix1) > 0
        ix1_names = displayAsString2(names[ix1])
        ix2_names = displayAsString2(names[ix2])
        println("\n  Dependent state variables x1:\n", ix1_names)
        println("\n  Independent state variables x2:\n", ix2_names)
        println("\n  Equations: ", eqx)
        for i = 1:length(eqx)
            printEquation(A[eqx[i],:], names)
        end
        println("\n  that shall be replaced by A1*x1 + A2*x2 = 0: ")
        printobj("   A1", A1)
        printobj("   A2", A2)
        print("\n")
        for i = 1:size(A1, 1)
            printEquation([A1 A2][i,:], [names[ix1];names[ix2]])
        end
    end   
end

end