using OffsetArrays
using Polynomials
using Dates

function compare_spins(sp,pos1,pos2)   # set to 1 when spin on pos1 is equal to pos2, -1 otherwise
    a = (sp >>> pos1) & 1 # true if spin is set on pos1
    b = (sp >>> pos2) & 1
    if a == b
        return 1
    else
        return -1
    end
end
function add_spin(p,nconf,k_max,pos,hx,hy,hz)
    dx    = OffsetArray{Int}(undef, 0:nconf-1)
    dy    = OffsetArray{Int}(undef, 0:nconf-1)
    sflip = OffsetArray{Int}(undef, 0:nconf-1)
    preg  = OffsetArray{Int128}(undef, 0:nconf-1)

    # coordinates of the new spin within the 3d lattice
    x  = pos % hx
    y  = (pos ÷ hx) % hy
    z  = pos ÷ (hx*hy)
    bp = x+hx*y # bitpos? Check real partition function!!!!

    # print("Pos=",pos," X=",x," Y=",y," Z=",z, " bp=",bp,"\n")
    
    if (x>0) # init spin functions as soon as lower x-neighbor is there
        @simd for spin = 0:nconf-1
            dx[spin] = compare_spins(spin,bp,bp-1)
        end
    else # no lower x-neighbor, use open boundary condition
        dx .= 0
    end
    if (y>0) # init spin functions as soon as lower y neighbors is there
        @simd for spin = 0:nconf-1
            dy[spin] = compare_spins(spin,bp,bp-hx)
        end
    else # no lower y-neighbor, use open boundary condition
        dy .= 0
    end
    if (z>0) # init spin functions as soon as lower y neighbors is there
        dz = 1
    else # no lower z-neighbor, use open boundary condition
        dz = 0
    end       

    if (x==hx-1) && (y==hy-1)
        @simd for s = 0:nconf-1
            sflip[s] =  s ⊻ (2^bp-1) # use Z2 symmetry
        end
    else
        @simd for s = 0:nconf-1
            sflip[s] = s ⊻ (2^bp) # flipped config
        end
    end     

    # update series expansion via recursion-step:
    for k in k_max:-1:1 # update all orders k>=1 
        @simd for s = 0:nconf-1
            pnew      = p[s,k  ] + p[s,k-2]*(dx[s]*dy[s] + dx[s]*dz + dy[s]*dz)
            term      = p[s,k-1]*(dx[s]+dy[s]+dz) + p[s,k-3]*(dx[s]*dy[s]*dz)
            preg[s]   = pnew - term
            p[s,k]    = pnew + term
        end
        @simd for s = 0:nconf-1
            p[s,k] = (p[s,k] + preg[sflip[s]]) ÷ 2
        end   
    end
end
function partition_function(k_max,hx,hy,hz)
    # calculate partition function of rectangular lattice (hx,hy,hz)
    nconf = 2^(hx*hy-1) # number of spin configs to keep trace
    p = OffsetArray{Int128}(undef, 0:nconf-1, -2:k_max) # partition function with two additional rows -2,-1 to simplify recursion
    p      .= 0 # reset partition function
    p[:,0] .= 1 # init partition function
    
    zk  = zeros(Int128,k_max+1) # Coefficients of Partition Function

    #  Now build up lattice spin by spin, layer by layer:
    for pos = 0:hx*hy*hz-1
        if hx*hy>=20
            print("Add Spin at pos: ", pos, " from", hx*hy*hz-1,"\n")
        end
        add_spin(p,nconf,k_max,pos,hx,hy,hz)
    end

    # sum into 1st partition function:
    for k = 0:k_max
        zk[k+1] = sum(p[:,k])/nconf
    end
    z = Polynomial(zk,:t)
    return z
end
function partition_simple(k_max,hx,hy,hz)
    # brute force calculate partition function of rectangular lattice (hx,hy,hz) to use as cross Check.
    # This routine will die miserably if asked to compute more than hx*hy*hz around 15
    
    nconf = 2^(hx*hy*hz) # number of spin configs to keep trace
    
    pf = Polynomial{Int128}([0],:t)
    t = Polynomial{Int128}([0,1],:t)

    for spin = 0:nconf-1
        p = Polynomial{Int128}([1],:t)
        for z = 0:hz-1
            for y = 0:hy-1
                for x = 0:hx-1
                    bp = x + hx*y + (hx*hy)*z # bitpos
                    if (x>0) # init spin functions as soon as lower x-neighbor is there
                        dx = compare_spins(spin,bp,bp-1)
                    else # no lower x-neighbor, use open boundary condition
                        dx = 0
                    end
                    if (y>0) # init spin functions as soon as lower y-neighbor is there
                        dy = compare_spins(spin,bp,bp-hx)
                    else # no lower x-neighbor, use open boundary condition
                        dy = 0
                    end
                    if (z>0) # init spin functions as soon as lower y-neighbor is there
                        dz = compare_spins(spin,bp,bp-(hx*hy))
                    else # no lower x-neighbor, use open boundary condition
                        dz = 0
                    end                
                    p = p * (1+dx*t) * (1+dy*t) * (1+dz*t)
                end
            end
        end
        coe = coeffs(p)
        if length(coe) > (k_max+1) # truncate polynomial to order k_max
            p = Polynomial(coe[1:k_max+1],:t)
        end
        # print("p=", p, "\n")
        pf = pf + p
    end
    return pf/nconf

end
function taylor_log(k_max, z)
    logz  = Polynomial{Rational{Int128}}([0],:t)
    for k=1:k_max
        logz = logz + (-1)^(k+1)//k * (z-1)^k # compute series expansion of log
        coe = coeffs(logz)
        if length(coe) > (k_max+1) # truncate polynomial to order k_max
            logz = Polynomial(coe[1:k_max+1],:t) 
        end
    end
    return logz
end

#----------------------------------------------------------------------
# Compute the High Temperature Expansion of the free energy per spin
# for the 3D Ising Model using a recursive Transfer Algorithm that
# generates exact partition functions over finite simple cubic lattices.
# Use convention that lattice (1,1,1) contains 2x2x2 spins.
#----------------------------------------------------------------------
const k_max = 26   # Set maximum order for series expansion
const sum_l = k_max÷2 # Maximum side-len of a single sublattice

z =    OffsetArray{Polynomial{Int128,:t}}(undef,0:sum_l,0:sum_l,0:sum_l)
logz = OffsetArray{Polynomial{Rational{Int128},:t}}(undef,0:sum_l,0:sum_l,0:sum_l)
phi =  OffsetArray{Polynomial{Rational{Int128},:t}}(undef,0:sum_l,0:sum_l,0:sum_l)

# 1. Compute list of required sublattices with (hx+hy+hz) <= k_max/2
lattices = []
for hx=0:sum_l
    for hy=0:sum_l
        for hz=0:sum_l
            if (hx+hy+hz)<=sum_l # sublattice is needed for the required order
                append!(lattices, [[hx,hy,hz]])
            end
        end
    end
end

# 2. Compute Z and Log(Z) for all sublatttices using cubic symmetry where possible
for lat in lattices
    (hx,hy,hz) = lat
    if (hx<=hy) && (hy<=hz) #  compute partition function if hx<=hy<=hz
        print("Compute Lattice: (",hx,",",hy,",",hz,")\n")
        z[hx,hy,hz] = partition_function(k_max,hx+1,hy+1,hz+1)
    else # copy partition function using the cubic symmetry.
        nx,ny,nz = sort([hx,hy,hz])
        #print("Copy Lattice: (",hx,",",hy,",",hz,") from (",nx,",",ny,",",nz,")\n")
        z[hx,hy,hz] = z[nx,ny,nz]
    end
    logz[hx,hy,hz] = taylor_log(k_max, z[hx,hy,hz]) # compute Taylor-Polynomial of logz
end

# 3. Compute phi recursively and combine into free energy density f according to Arisue et.al.
f = Polynomial{Rational{Int128}}([0],:t)
for lat in lattices
    (hx,hy,hz) = lat
    #print("Compute Phi: (",hx,",",hy,",",hz,")\n")
    phi[hx,hy,hz] = logz[hx,hy,hz]
    # Now loop over all smaller sublattices of (hx,hy,hz)
    for nx=0:hx
        for ny=0:hy
            for nz=0:hz
                if ((nx+ny+nz) != (hx+hy+hz)) # (nx,ny,nz) is a smaller sublattice to be embedded in (hx,hy,hz)
                    phi[hx,hy,hz] = phi[hx,hy,hz] - (hx-nx+1)*(hy-ny+1)*(hz-nz+1)*phi[nx,ny,nz]
                end
            end
        end
    end
    global f = f + phi[hx,hy,hz]
end

# Compute free energy density per spin using Guttmann/Enting weights:
g = Polynomial{Rational{Int128}}([0],:t)
for lat in lattices
    (hx,hy,hz) = lat
    sum = hx+hy+hz+3
    s = (k_max+6)/2
    if sum == s
        w = 1
    elseif sum == s-1
        w = -5
    elseif sum == s-2
        w = 10
    elseif sum == s-3
        w = -10
    elseif sum == s-4
        w = 5
    elseif sum == s-5
        w = -1
    else
        w = 0
    end
    #print("lattice (",hx,",",hy,",",hz,") f=f+ ",w,"*",logz[hx,hy,hz],"\n")
    global g = g + w*logz[hx,hy,hz]  # Guttmann / Enting
end