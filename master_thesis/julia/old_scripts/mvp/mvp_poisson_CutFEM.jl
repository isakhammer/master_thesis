# %% Import required modules
using Gridap
using GridapEmbedded


# %% Manufactured solution
# Provides a manufactured solution which is 0 on the unit circle
# u(x) = (x[1]^2 + x[2]^2  - 1)*sin(2π*x[1])*cos(2π*x[2])
u_ex(x) = 1 - x[1]^2 - x[2]^2
f(x) = 4
∇u_ex(x) = VectorValue(-2*x[1], -2*x[2])

# Note that we have used the constructor `VectorValue` to build the
# vector that represents the gradient. However, we still need a final
# trick. We need to tell the Gridap library that the gradient of the
# function `u` is available in the function `∇u` (at this moment `u`
# and `∇u` are two standard Julia functions without any connection
# between them). This is done by adding an extra method to the
# function `gradient` (aka `∇`) defined in Gridap:
import Gridap: ∇
∇(::typeof(u_ex)) = ∇u_ex
∇(u_ex) === ∇u_ex

# %% Cartesian background model and embedded boundary

# Background model
domain = (-1.11, 1.11, -1.11, 1.11)
pmin = Point(-1.11, -1.11)
pmax = Point(1.11, 1.11)
partition = (10,10)
bgmodel = CartesianDiscreteModel(pmin, pmax, partition)

# Implicit geometry
# TODO: Define own level set function via AnalyticalGeometry
R  = 1.0
geo = disk(R)

# Cut the background model
cutgeo = cut(bgmodel, geo)

# Set up interpolation mesh and function spaces
Ω_act = Triangulation(cutgeo, ACTIVE)

# Construct function spaces
order = 1
V = TestFESpace(Ω_act, ReferenceFE(lagrangian, Float64, order),conformity=:H1)
U = TrialFESpace(V)

# Set up integration meshes, measures and normals
Ω = Triangulation(cutgeo, PHYSICAL)
Γ = EmbeddedBoundary(cutgeo)
Fg = GhostSkeleton(cutgeo)

# Set up integration measures
degree = 2*order
dΩ   = Measure(Ω, degree)
dΓ   = Measure(Γ, degree)
dFg  = Measure(Fg, degree)

# Set up normal vectors
n_Γ = get_normal_vector(Γ)
n_Fg = get_normal_vector(Fg)

# Write out models and computational domains for inspection
dirname = "figures/poission_cutfem"
mkpath(dirname)
writevtk(bgmodel, dirname*"/bgmodel")
writevtk(Ω, dirname*"/Omega")
writevtk(Ω_act, dirname*"/Omega_act")
writevtk(Γ, dirname*"/Gamma")
writevtk(Fg, dirname*"/Fg")

# Define weak form
# Nitsche parameter
γd = 10.0

# Ghost penalty parameter
γg = 0.1

# Mesh size
h = (pmax - pmin)[1]/partition[1]

# Define bilinear form
a(u,v) =
    ∫( ∇(u)⋅∇(v) ) * dΩ  +
    ∫( (γd/h)*u*v  - u*(n_Γ⋅∇(v)) - (n_Γ⋅∇(u))*v ) * dΓ +
    ∫( (γg*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg

# Define linear form
l(v) =
    ∫( f*v ) * dΩ +
    ∫( u_ex*( (γg/h)*v - (n_Γ⋅∇(v)) )  ) * dΓ

# FE problem
op = AffineFEOperator(a,l,U,V)
uh = solve(op)

# Postprocess
outputfile = dirname*"/PoissonCutFEM"
if outputfile !== nothing
    writevtk(Ω,outputfile,cellfields=["uh"=>uh, "u_ex"=>u_ex])
end
