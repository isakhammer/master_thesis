
using Test
using Gridap

# Analytical manufactured solution
α = 1

# u(x) = x[1]*(x[1]-1)*x[2]*(x[2]-1)
u(x) = cos(x[1])*cos(x[2])

f(x) = Δ(Δ(u))(x)+ α*u(x)
g(x) = Δ(u)(x)


# Domain
L = 2*π

@test f(VectorValue(0.5,0.5)) == ( 4+α )*u(VectorValue(0.5,0.5))
@test g(VectorValue(0.5,0.5)) == -2*u(VectorValue(0.5,0.5))

domain = (0,L,0,L)
partition = (90,90)
model = CartesianDiscreteModel(domain,partition)

# FE space
order = 2
V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order),dirichlet_tags="boundary")
U = TrialFESpace(V,u)

# Triangulation
Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(model)
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
dΛ = Measure(Λ,degree)
nΓ = get_normal_vector(Γ)
nΛ = get_normal_vector(Λ)

# Weak form
const h = (domain[2]-domain[1]) / partition[1]
const γ = 1
a(u,v) = ∫( Δ(u)*Δ(v) + α* u⋅v )dΩ +
         ∫( - mean(Δ(u))*jump(∇(v)⋅nΛ) - jump(∇(u)⋅nΛ)*mean(Δ(v)) + γ/h*jump(∇(u)⋅nΛ)*jump(∇(v)⋅nΛ) )dΛ
l(v) = ∫( v*f )dΩ + ∫( g*(∇(v)⋅nΓ) )dΓ
op = AffineFEOperator(a,l,U,V)

uₕ = solve(op)

# Error
e = u - uₕ
l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))
el2 = l2(e)
eh1 = h1(e)
tol = 1.0e-10
@test el2 < tol
@test eh1 < tol

