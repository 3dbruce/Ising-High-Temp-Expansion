using OffsetArrays
using Printf
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
    pnew  = OffsetArray{Int128}(undef, 0:nconf-1)

    # coordinates of the new spin within the 3d lattice
    x  = pos % hx
    y  = (pos ÷ hx) % hy
    z  = pos ÷ (hx*hy)
    bp = x+hx*y # bitpos?

    print("Pos=",pos," X=",x," Y=",y," Z=",z, " bp=",bp,"\n")

    
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

    @simd for s = 0:nconf-1
        sflip[s] = s ⊻ (2^bp) # flipped config
        #print("nconf=",nconf," sflip=",sflip[s],"\n")
    end
    if (x==hx-1) && (y==hy-1)
        @simd for s = 0:nconf-1
            sflip[s] =  s ⊻ (2^bp-1) # use Z2 symmetry
            #print("nconf=",nconf," sflip=",sflip[s],"\n")
        end
    end     

    # update series expansion via recursion-step:
    for k in k_max:-1:1 # update all orders k>=1 
        @simd for s=0:nconf-1
            pnew[s] = p[s,k]+p[s,k-1]*(dx[s]+dy[s]+dz)+p[s,k-2]*(dx[s]*dy[s]+dx[s]*dz+dy[s]*dz)+p[s,k-3]*(dx[s]*dy[s]*dz)
            preg[s] = p[s,k]-p[s,k-1]*(dx[s]+dy[s]+dz)+p[s,k-2]*(dx[s]*dy[s]+dx[s]*dz+dy[s]*dz)-p[s,k-3]*(dx[s]*dy[s]*dz)                     
        end
        @simd for s=0:nconf-1
            p[s,k] = ( pnew[s] + preg[sflip[s]] ) ÷ 2
        end
    end
end

# Set Parameter for Simple Cubic Lattice
const k_max    = 24 # maximum order for series expansion
const hx,hy,hz = (4,5,4) # cube size
const nconf    = 2^(hx*hy-1) # number of spin configs to keep track

part1  = OffsetArray{Int128}(undef,0:k_max)
part1 .= 0

p = OffsetArray{Int128}(undef, 0:nconf-1, -2:k_max) # partition function with two additional rows -2,-1 to simplify recursion
    
# init partition function
p      .= 0 # reset partition function
p[:,0] .= 1   
    
#  Now build up lattice spin by spin, layer by layer:
for pos = 0:hx*hy*hz-1
    add_spin(p,nconf,k_max,pos,hx,hy,hz)
end

# sum into 1st partition function:
for k = 0:k_max
    part1[k] = sum(p[:,k])/2
end 

part1