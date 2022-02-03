include("ising_routines.jl")

#----------------------------------------------------------------------
# Compute the High Temperature Expansion of the free energy per spin
# for the 3D Ising Model using a recursive Transfer Algorithm that
# generates exact partition functions over finite simple cubic lattices.
# Here I use the Arisue et.al. notation of a lattice (1,1,1) having 2x2x2 spins.
#----------------------------------------------------------------------
const k_max = 16   # Set maximum order for series expansion
const sum_l = k_max÷2 # Maximum side-len of a single sublattice

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

# 2. Compute Log(Z) for all sublatttices using cubic symmetry where possible
for lat in lattices
    (hx,hy,hz) = lat
    if (hx<=hy) && (hy<=hz) #  compute partition function if hx<=hy<=hz
        print("Compute Lattice: (",hx,",",hy,",",hz,")\n")
        logz[hx,hy,hz] = logz_3d_sc_rtm(k_max,hx+1,hy+1,hz+1) # use +1 to translate arisue notation
    else # copy partition function using the cubic symmetry.
        nx,ny,nz = sort([hx,hy,hz])
        #print("Copy Lattice: (",hx,",",hy,",",hz,") from (",nx,",",ny,",",nz,")\n")
        logz[hx,hy,hz] = logz[nx,ny,nz]
    end
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
    global f += phi[hx,hy,hz]
end
print("Arisue/Fujiwara: f=", f,"\n")

# Compute free energy density per spin using Guttmann/Enting weights:
f = Polynomial{Rational{Int128}}([0],:t)
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
    global f += w*logz[hx,hy,hz]  # Guttmann / Enting
end
print("Enting/Guttmann: f=", f, "\n")