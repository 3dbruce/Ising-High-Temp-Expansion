using OffsetArrays
using Printf
using Dates
using Base.Threads

function compare_spins(sp,pos1,pos2)   # set to 1 when spin on pos1 is equal to pos2, -1 otherwise
    a = (sp >>> pos1) & 1 # true if spin is set on pos1
    b = (sp >>> pos2) & 1
    if a == b
        return 1
    else
        return -1
    end
end

function add_spin(p,nconf,k_max,pos,hx,hy,hz,wgt)      
    dx    = OffsetArray{Int}(undef, 0:nconf-1)
    dy    = OffsetArray{Int}(undef, 0:nconf-1)
    sflip = OffsetArray{Int}(undef, 0:nconf-1)
    preg  = OffsetArray{Int128}(undef, 0:nconf-1)

    j=hz-1-(pos%hz) # Bitpos. at position pos-1

    if (pos>=hx) # calculate spin contribution from lower x-neigbor 
        @simd for spin = 0:nconf-1
            dx[spin] = wgt[1] * compare_spins(spin,j,(j+hx)%hz)
        end
    else # lower x-neighbor still outside of lattice, set to 0 due to open boundary conditions
        dx   .= 0
    end
    if (pos>=hy) # calculate spin contribution from lower y-neigbor 
        @simd for spin = 0:nconf-1
            dy[spin] = wgt[2] * compare_spins(spin,j,(j+hy)%hz)
        end
    else # lower y-neighbor still outside of lattice, set to 0 due to open boundary conditions
        dy   .= 0
    end
    if (pos>=hz) # calculate spin contribution from lower z-neigbor 
        dz = wgt[3] # does not depend on spin-state
    else # lower z-neighbor still outside of lattice, set to 0 due to open boundary conditions
        dz = 0
    end
  
    if j!=(hz-1) # calculate spin-config with flipped lower z-neighbor spin on pos j
        @simd for spin = 0:nconf-1
            sflip[spin] = spin⊻(2^j) 
        end
    else # highest spin, flip all other spins to use Z(2) Symm.
        for spin = 0:nconf-1
            sflip[spin] = (2^hz)-1 - spin⊻(2^j)  
        end
    end

    # update series expansion via recursion-step:
    for k in k_max:-1:1 # update all orders k>=1 
        @simd for s = 0:nconf-1
            pnew      = p[s,k] + p[s,k-2]*(dx[s]*dy[s]+dx[s]*dz+dy[s]*dz)
            term      = p[s,k-1]*(dx[s]+dy[s]+dz) + p[s,k-3]*(dx[s]*dy[s]*dz)
            preg[s]   = pnew - term
            p[s,k]    = pnew + term
        end
        @simd for s = 0:nconf-1
            p[s,k] = (p[s,k] + preg[sflip[s]]) ÷ 2
        end     
    end
end

# Set Parameter for Lattice
const k_max  = 24 # maximum order for series expansion
const layers = 12 # number of layers of hz-spins
const hx,hy,hz = (9,11,13) #
const nconf = 2^(hz-1) # half of required spin-configs due to use of Z(2) Symmetrie

part1  = OffsetArray{Int128}(undef,8,0:k_max)
part2  = OffsetArray{Int128}(undef,8,0:k_max)

p = OffsetArray{Int128}(undef, 0:nconf-1, -2:k_max) # partition function with two additional rows -2,-1 to simplify recursion

# define weight factors to eliminate unphysical loops
weight = [[1,1,1],[-1,1,1],[1,-1,1],[1,1,-1],[-1,-1,1],[1,-1,-1],[-1,1,-1],[-1,-1,-1]]

print("Lattice = (",hx,",",hy,",",hz,")\n")
for w=1:8 # loop over weigths to be able to cancel unphysical FS-loops
 
    tim = Dates.format(now(), "HH:MM:SS")
    print("w=",w," at ",tim,"\n")
    
    # init partition function
    p .= 0
    p[:,0] .= 1   
    
    #  Now build up lattice spin by spin, layer by layer:
    for pos = 0:layers*hz-1
        # print("pos=",pos,"\n")
        add_spin(p,nconf,k_max,pos,hx,hy,hz,weight[w])
    end

    # sum into 1st partition function:
    for k = 0:k_max
        part1[w,k] = sum(p[:,k])
    end 

    # add one additional spin to eliminate surface effects later:
    add_spin(p,nconf,k_max,layers*hz,hx,hy,hz,weight[w])

    # sum into 2nd partition function:
    for k = 0:k_max
        part2[w,k] = sum(p[:,k])
    end    
end
part1