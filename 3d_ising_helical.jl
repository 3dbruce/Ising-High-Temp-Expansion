include("ising_routines.jl")

#----------------------------------------------------------------------
# Compute the High Temperature Expansion of the free energy per spin
# for the 3D Ising Model using a recursive Transfer Algorithm that
# generates exact partition functions over finite helical lattices.
#----------------------------------------------------------------------

const k_max = 24
const layers = 12
const lattices = [[-3,9,11,13],[3,1,12,14],[-3,9,14,16],[-3,5,15,16],[3,7,15,16],[-3,10,13,17],[3,5,15,17],[-3,14,15,17],[3,11,16,17],[3,14,16,17],[-1,9,17,19],[-2,9,16,20],[-1,5,17,20],[1,5,19,20],[-2,16,17,21],[5,10,19,21],[2,16,20,21],[-2,1,18,22],[2,17,21,22]]
#const lattices =[[2,9,11,13]]

f = Polynomial{Rational{Int128}}([0],:t)
for lat in lattices
    (wgt,hx,hy,hz) = lat
    logz = logz_3d_hel_rtm(k_max,layers,hx,hy,hz)
    global f += wgt*logz
end
coe  = coeffs(f) // 2 # divide by total weight
f   = Polynomial(coe,:t)
print("f=", f)