using OffsetArrays
using Polynomials
using Dates

function compare_spins(sp,pos1,pos2)
    # set to 1 when spin on pos1 is equal to pos2, -1 otherwise
    a = (sp >>> pos1) & 1 # true if spin is set on pos1
    b = (sp >>> pos2) & 1
    if a == b
        return 1
    else
        return -1
    end
end
function trunc_poly(k_max, p)
    coe = coeffs(p)
    if length(coe) > (k_max+1) # truncate polynomial to order k_max
        p = Polynomial(coe[1:k_max+1],:t)
    end
    return p
end
function taylor_log(k_max, z)
    # computes the Taylor expansion of log Z up to the order k_max 
    logz = Polynomial{Rational{Int128}}([0],:t)
    zpow = Polynomial{Int128}([1],:t)
    for k=1:k_max
        #logz = logz + (-1)^(k+1)//k * (z-1)^k # compute series expansion of log
        zpow = zpow * (z-1)
        logz += (-1)^(k+1)//k * zpow
        logz = trunc_poly(k_max, logz)
        zpow = trunc_poly(k_max, zpow)
    end
    return logz
end
function add_3d_sc_spin(p,nconf,k_max,pos,hx,hy,hz)
    # adds new spin to a 3d single cubic lattice with size (hx,hy,hz) at position pos and returns
    # the logZ of the partition function Z calculated by an recursive transfer-matrix algorithm.
    # Required memory is proportional to 2^(hx*hy-1) and currently about 8GB for (hx*hy)=(5*5)
    
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
function add_3d_hel_spin(p,nconf,k_max,pos,hx,hy,hz,wgt)
    # adds new spin to a 3d helical lattice with structure (hx,hy,hz) at position pos and returns
    # the logZ of the partition function Z calculated by an recursive transfer-matrix algorithm.
    # Required memory is proportional to 2^(hz-1) and currently about 1GB for hz=22
    # weight structure wgt is used to cancel unphysical FS loops with an odd number of links.

    dx    = OffsetArray{Int}(undef, 0:nconf-1)
    dy    = OffsetArray{Int}(undef, 0:nconf-1)
    sflip = OffsetArray{Int}(undef, 0:nconf-1)
    preg  = OffsetArray{Int128}(undef, 0:nconf-1)

    j=hz-1-(pos%hz) # Bitpos. at position pos-1

    if (pos>=hx) # calculate spin contribution from lower x-neigbor 
        for spin = 0:nconf-1
            dx[spin] = wgt[1] * compare_spins(spin,j,(j+hx)%hz)
        end
    else # lower x-neighbor still outside of lattice, set to 0 due to open boundary conditions
        dx   .= 0
    end
    if (pos>=hy) # calculate spin contribution from lower y-neigbor 
        for spin = 0:nconf-1
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
        for spin = 0:nconf-1
            sflip[spin] = spin ⊻ (2^j) 
        end
    else # highest spin, flip all other spins to use Z(2) Symm.
        for spin = 0:nconf-1
            sflip[spin] = (2^hz)-1 - spin ⊻ (2^j)
        end
    end

    # update series expansion via recursion-step:
    for k in k_max:-1:1 # update all orders k>=1 
        for s = 0:nconf-1
            pnew      = p[s,k] + p[s,k-2]*(dx[s]*dy[s]+dx[s]*dz+dy[s]*dz)
            term      = p[s,k-1]*(dx[s]+dy[s]+dz) + p[s,k-3]*(dx[s]*dy[s]*dz)
            preg[s]   = pnew - term
            p[s,k]    = pnew + term
        end
        for s = 0:nconf-1
            p[s,k] = (p[s,k] + preg[sflip[s]]) ÷ 2
        end     
    end
end
function logz_3d_sc_simple(k_max,hx,hy,hz)
    # brute force calculation of the log of the partition function for rectangular lattices (hx,hy,hz) to
    # use as cross-check. This routine will fail miserably if asked to compute more than around hx*hy*hz=12
    
    nconf = 2^(hx*hy*hz) # number of spin configs to keep track
    
    pf   = Polynomial{Int128}([0],:t)
    t    = Polynomial{Int128}([0,1],:t)
    logz = Polynomial{Rational{Int128}}([0],:t)

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
        pf = pf + trunc_poly(k_max, p)
    end
    coe  = coeffs(pf) .÷ nconf # divide by trivial factor 2^(hx*hy*hz)
    pf   = Polynomial(coe,:t)
    logz = taylor_log(k_max, pf)
    return logz
end
function logz_3d_sc_rtm(k_max,hx,hy,hz)
    # calculate the log of the partition function of a 3d simple cubic rectangular lattice (hx,hy,hz)
    # using a recursive transfer matrix algorithm that keeps track of the last exposed hx*hy spins only.
    # Memory requirement scales with 2^(hx*hy) 
    nconf = 2^(hx*hy-1) # number of spin configs to keep track (a factor of 2 is saved because the Z2 symmetry is used)
    
    p = OffsetArray{Int128}(undef, 0:nconf-1, -2:k_max) # partition function with two additional rows -2,-1 to simplify recursion
    p      .= 0 # reset partition function
    p[:,0] .= 1 # init partition function
    
    zk  = zeros(Int128,k_max+1) # Coefficients of Partition Function

    #  Now build up lattice spin by spin, layer by layer:
    for pos = 0:hx*hy*hz-1
        if hx*hy>=20 # Add progress indicator output
            print("Add Spin at pos: ", pos, " from", hx*hy*hz-1,"\n")
        end
        add_3d_sc_spin(p,nconf,k_max,pos,hx,hy,hz)
    end
    
    # sum into partition function:
    for k = 0:k_max
        zk[k+1] = sum(p[:,k])/nconf
    end
    z = Polynomial(zk,:t)

    # Calculate and return logZ
    logz = taylor_log(k_max,z) 
    return logz
end
function logz_3d_hel_rtm(k_max, layers, hx,hy,hz)
    # compute logZ of the 3d Ising model using helical lattices with structure (hx,hy,hz)
    nconf = 2^(hz-1) # half of required spin-configs due to use of Z(2) Symmetrie

    p = OffsetArray{Int128}(undef, 0:nconf-1, -2:k_max) # partition function with two additional rows -2,-1 to simplify recursion
    zk1  = zeros(Int128,k_max+1) # Coefficients of Partition Function1
    zk2  = zeros(Int128,k_max+1) # Coefficients of Partition Function2
    f = Polynomial{Rational{Int128}}([0],:t)

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
            if hz>= 20
                print("pos=",pos,"\n") # Add progress indicator for large lattices
            end
            add_3d_hel_spin(p,nconf,k_max,pos,hx,hy,hz,weight[w])
        end

        # sum into 1st partition function:
        for k = 0:k_max
            zk1[k+1] = sum(p[:,k]) ÷ nconf
        end
        z1 = Polynomial(zk1,:t)

        # add one additional spin to eliminate surface effects later:
        add_3d_hel_spin(p,nconf,k_max,layers*hz,hx,hy,hz,weight[w])

        # sum into 2nd partition function:
        for k = 0:k_max
            zk2[k+1] = sum(p[:,k]) ÷ nconf
        end
        z2 = Polynomial(zk2,:t)

        # compute free energy as difference of log Z
        f += taylor_log(k_max, z2) - taylor_log(k_max, z1)
    end
    coe  = coeffs(f) // 8 # divide by trivial weight factor
    f   = Polynomial(coe,:t)
    return f
end