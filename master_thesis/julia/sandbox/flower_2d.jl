using Gridap
using GridapDistributed
using GridapEmbedded

# Embedded geometry
L = 1.51
H = 1.51

r0 = 0.5
r1 = 0.15
function ls_flower(x)
    theta = atan(x[1], x[2])
    (x[1]^2 + x[2]^2)^0.5 - r0 - r1*cos(5.0*theta)
end

flower = AnalyticalGeometry(x-> ls_flower(x))

# Total background mesh
n = 50
pmin = 1.1*Point(-L/2.0,-H/2.0)
pmax = pmin+1.1*Point(L, H)
partition = (n, n)
dim = length(partition)
bgmodel = CartesianDiscreteModel(pmin, pmax, partition) |> simplexify

# Physical domain bgmodel - ellipse
cutgeo = cut(bgmodel, !flower)

#  Active/interpolation mesh
Ω_act = Triangulation(cutgeo, ACTIVE)

# Set up integration meshes
Ω = Triangulation(cutgeo, PHYSICAL)
Γ = EmbeddedBoundary(cutgeo)
Fg = GhostSkeleton(cutgeo)
Fint = SkeletonTriangulation(Ω_act)

# Write them out for inspection
writevtk(bgmodel, "bgmodel")
writevtk(Ω_act, "Omega_act")
writevtk(Ω, "Omega")
writevtk(Γ, "Gamma")
writevtk(Fg, "Gamma_g")
writevtk(Fg, "Gamma_g")
writevtk(Fint, "Gamma_int")