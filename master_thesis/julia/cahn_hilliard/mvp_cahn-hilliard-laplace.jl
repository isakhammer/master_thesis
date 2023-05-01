## 
using Gridap
using Gridap.Algebra

## Cahn-hilliard
# ε = 1/100
ε = 1
# Gibb's potential
ψ(u) = u*(1-u^2)

u_ex(x, t::Real) = cos(x[1])*cos(x[2])*exp(-4*t/ε^2)
u_ex(t) = x -> u_ex(x,t)

f_ex(x, t::Real) = 0 


## 
L=2π
n = 10
h = L/n
domain2D = (0, L, 0, L)
partition2D = (n, n)
model = CartesianDiscreteModel(domain2D,partition2D) |> simplexify
resultdir = "cahn-hilliard-results/"
mkpath(resultdir)
writevtk(model, joinpath(resultdir,"model"))

## Function spaces
order = 2
V = TestFESpace(model, ReferenceFE(lagrangian,Float64, order), conformity=:H1)
U = TrialFESpace(V)
Ω = Triangulation(model)
Λ = SkeletonTriangulation(model)
Γ = BoundaryTriangulation(model)

## Forms
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
dΛ = Measure(Λ,degree)

n_Λ = get_normal_vector(Λ)
n_Γ = get_normal_vector(Γ)

# M+ dt A
γ = 1.5*order*( order+1)
τ = ε^2/10
a(u,v) = ∫(u*v)*dΩ + τ/ε^2*( ∫(Δ(v)⊙Δ(u))dΩ
            + ∫(-mean(Δ(v))⊙jump(∇(u)⋅n_Λ) - mean(Δ(u))⊙jump(∇(v)⋅n_Λ) + (γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
            + ∫(-Δ(v)⊙∇(u)⋅n_Γ - Δ(u)⊙∇(v)⋅n_Γ + (γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
        )

# l(u, v) = ∫(τ*f*v + u*v)*dΩ
l(u, v) = ∫(u*v)*dΩ
l(u) = v -> l(u,v)

createpvd(resultdir*"ch-solution") do pvd

    ## Set up linear algebra system
    A = assemble_matrix(a, U, V)
    lu = LUSolver()

    ## time loop
    t0 = 0.0
    T = 100*τ
    t = t0

    println(t)
    # u_dof_vals = rand(Float64, num_free_dofs(U))
    uh = interpolate_everywhere(u_ex(0),U)
    pvd[t] = createvtk(Ω, resultdir*"ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

    # Solve for first time step
    b = assemble_vector(l(uh), V)
    op = AffineOperator(A, b)

    u_dof_vals = get_free_dof_values(uh)
    cache = solve!(u_dof_vals, lu, op)
    uh = FEFunction(U, u_dof_vals)

    t += τ
    pvd[t] = createvtk(Ω, resultdir*"ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

    # Remaining while loop
    while t < T
        b = assemble_vector(l(uh), V)

        u_dof_vals = get_free_dof_values(uh)
        op = AffineOperator(A, b)
        cache = solve!(u_dof_vals, lu, op, cache, false)
        uh = FEFunction(U, u_dof_vals)

        t += τ
        pvd[t] = createvtk(Ω, resultdir*"ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

    end
end 
