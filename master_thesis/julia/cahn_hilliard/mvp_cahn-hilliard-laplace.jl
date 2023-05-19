## 
using Gridap
using Gridap.Algebra

## Cahn-hilliard
ε = 1/30
# ε = 1
# Gibb's potential
ψ(u) = mean(u)*(1-mean(u)*mean(u))
# ψ(u) = u*(1-u^2)
# ψ(u) = u

u_ex(x, t::Real) = cos(x[1])*cos(x[2])*exp(-(4*ε^2 + 2)*t)
u_ex(t) = x -> u_ex(x,t)
f_ex(x, t::Real) = 0 

## 
L=2π
n = 64
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
a(u,v) = ∫(u*v)*dΩ + τ*ε^2*( 
              ∫(Δ(v)⊙Δ(u))dΩ
            + ∫(-mean(Δ(v))⊙jump(∇(u)⋅n_Λ) - mean(Δ(u))⊙jump(∇(v)⋅n_Λ) + (γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
            + ∫(-Δ(v)⊙∇(u)⋅n_Γ - Δ(u)⊙∇(v)⋅n_Γ + (γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
        )

# l(u, v) = ∫(τ*f*v + u*v)*dΩ
l(u, v) =  ∫(u*v)*dΩ + τ * ( 
            ∫(ψ(u)*Δ(v))*dΩ 
          - ∫(ψ(mean(u))*jump(∇(v)⋅n_Λ))*dΛ
          - ∫(ψ(u)*∇(v)⋅n_Γ )*dΓ
          )
l(u) = v -> l(u,v)

createpvd(resultdir*"ch-solution") do pvd

    ## time loop
    t0 = 0.0
    T = 10000*τ
    Nt_max = convert(Int64, ceil((T - t0)/τ))
    Nt = 0
    t = t0

    # Maximal number of Picard iterations
    kmax = 1

    # Initial data
    u_dof_vals = (rand(Float64, num_free_dofs(U)) .-0.5)*2.0
    # uh = interpolate_everywhere(u_ex(0),U)
    uh = FEFunction(U, u_dof_vals)
    pvd[t] = createvtk(Ω, resultdir*"ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

    println("========================================")
    println("Solving Cahn-Hilliard with t0 = $t0, T = $T and time step τ = $τ with Nt_max = $Nt_max timesteps")
    println("========================================")

    ## Set up linear algebra system
    A = assemble_matrix(a, U, V)
    lu = LUSolver()
    cache = nothing

    # Time loop
    while t < T
        Nt += 1
        t += τ
        println("----------------------------------------")
        println("Solving Cahn-Hilliard for t = $t, step $(Nt)/$(Nt_max)")
        k = 0
        while k < kmax
            k += 1
            println("Iteration k = $k")
            b = assemble_vector(l(uh), V)
            op = AffineOperator(A, b)
            cache = solve!(u_dof_vals, lu, op, cache, isnothing(cache))
            uh = FEFunction(U, u_dof_vals)
        end
        println("----------------------------------------")
        pvd[t] = createvtk(Ω, resultdir*"ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])
    end
end 
