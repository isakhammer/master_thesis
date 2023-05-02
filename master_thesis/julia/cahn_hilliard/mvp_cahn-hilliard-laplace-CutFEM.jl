##
using Gridap
using Gridap.Algebra
using GridapEmbedded

## Cahn-hilliard
ε = 1/20
# ε = 1
# Gibb's potential
ψ(u) = mean(u)*(1-mean(u)*mean(u))
# ψ(u) = u*(1-u^2)
# ψ(u) = u

u_ex(x, t::Real) = cos(x[1])*cos(x[2])*exp(-(4*ε^2 + 2)*t)
u_ex(t) = x -> u_ex(x,t)
f_ex(x, t::Real) = 0

##
L=1.11
n = 2^7
h = L/n
# domain2D = (0, L, 0, L)
# partition2D = (n, n)
# model = CartesianDiscreteModel(domain2D,partition2D) |> simplexify

pmin = Point(-L, -L)
pmax = Point(L, L)
partition = (n,n)
bgmodel = CartesianDiscreteModel(pmin, pmax, partition)

resultdir = "cahn-hilliard-results/"
if isdir(resultdir)
    rm(resultdir, force=true, recursive=true)
end
mkpath(resultdir)
writevtk(bgmodel, joinpath(resultdir,"model"))

# Implicit geometry
R  = 1.0
geo = disk(R)

# Cut the background model
cutgeo = cut(bgmodel, geo)
cutgeo_facets = cut_facets(bgmodel,geo)

# Set up interpolation mesh and function spaces
Ω_act = Triangulation(cutgeo, ACTIVE)

## Function spaces
order = 2
# Construct function spaces
V = TestFESpace(Ω_act, ReferenceFE(lagrangian, Float64, order), conformity=:H1)
U = TrialFESpace(V)


# Set up integration meshes, measures and normals
Ω = Triangulation(cutgeo, PHYSICAL)
Γ = EmbeddedBoundary(cutgeo)
Λ = SkeletonTriangulation(cutgeo_facets)
Fg = GhostSkeleton(cutgeo)

# Set up integration measures
degree = 2*order
dΩ   = Measure(Ω, degree)
dΓ   = Measure(Γ, degree)
dΛ   = Measure(Λ, degree) # F_int
dFg  = Measure(Fg, degree)

n_Λ = get_normal_vector(Λ)
n_Γ = get_normal_vector(Γ)

# Set up normal vectors
n_Γ = get_normal_vector(Γ)
n_Λ = get_normal_vector(Λ)
n_Fg = get_normal_vector(Fg)

function jump_nn(u,n)
    return ( n.plus⋅ ∇∇(u).plus⋅ n.plus - n.minus ⋅ ∇∇(u).minus ⋅ n.minus )
end


# M+ dt A
γ = 1.5*order*( order+1)
a_CIP(u,v) = ∫(u*v)*dΩ + τ*ε^2*( ∫(Δ(v)⊙Δ(u))dΩ
            + ∫(-mean(Δ(v))⊙jump(∇(u)⋅n_Λ) - mean(Δ(u))⊙jump(∇(v)⋅n_Λ) + (γ/h)⋅jump(∇(u)⋅n_Λ)⊙jump(∇(v)⋅n_Λ))dΛ
            + ∫(-Δ(v)⊙∇(u)⋅n_Γ - Δ(u)⊙∇(v)⋅n_Γ + (γ/h)⋅ ∇(u)⊙n_Γ⋅∇(v)⊙n_Γ )dΓ
        )

γg1 = 10/2
γg2 = 0.1
g(u,v) = h^(-2)*( ∫( (γg1*h)*jump(n_Fg⋅∇(u))*jump(n_Fg⋅∇(v)) ) * dFg +
                   ∫( (γg2*h^3)*jump_nn(u,n_Fg)*jump_nn(v,n_Fg) ) * dFg)

a(u,v) = a_CIP(u,v) + g(u,v)

# l(u, v) = ∫(τ*f*v + u*v)*dΩ
l(u, v) =  ∫(u*v)*dΩ + τ * (
            ∫(ψ(u)*Δ(v))*dΩ
          - ∫(ψ(mean(u))*jump(∇(v)⋅n_Λ))*dΛ
          - ∫(ψ(u)*∇(v)⋅n_Γ )*dΓ
          )
l(u) = v -> l(u,v)

createpvd(resultdir*"CutFEM-ch-solution") do pvd

    τ = ε^2/10^13

    ## time loop
    t0 = 0.0
    # T = 1000*τ
    # Nt_max = convert(Int64, ceil((T - t0)/τ))
    T = 10^-12
    t = t0

    # Initial data
    u_dof_vals = (rand(Float64, num_free_dofs(U)) .-0.5)*2.0
    # uh = interpolate_everywhere(u_ex(0),U)
    uh = FEFunction(U, u_dof_vals)
    pvd[t] = createvtk(Ω, resultdir*"CutFEM-ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

    println("========================================")
    Nt = 1
    t += τ
    println("Solving CutFEM Cahn-Hilliard with t0 = $t0, T = $T and time step τ = $τ")
    println("----------------------------------------")

    ## Set up linear algebra system
    A = assemble_matrix(a, U, V)
    lu = LUSolver()

    # Solve for first time step
    b = assemble_vector(l(uh), V)
    op = AffineOperator(A, b)
    u_dof_vals = get_free_dof_values(uh)
    cache = solve!(u_dof_vals, lu, op)
    uh = FEFunction(U, u_dof_vals)

    pvd[t] = createvtk(Ω, resultdir*"CutFEM-ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])

    # Remaining while loop
    while t < T
        Nt += 1
        t += τ
        # println("Solving Cahn-Hilliard for t = $t, step $(Nt)/$(Nt_max) and timestep $τ ")
        println("Solving Cahn-Hilliard for t = $t of T=$T, step $(Nt) and timestep $τ ")
        println("----------------------------------------")
        if Nt%50==0
            τ *=2
        end
        b = assemble_vector(l(uh), V)
        u_dof_vals = get_free_dof_values(uh)
        op = AffineOperator(A, b)
        cache = solve!(u_dof_vals, lu, op, cache, false)
        uh = FEFunction(U, u_dof_vals)

        pvd[t] = createvtk(Ω, resultdir*"CutFEM-ch-solution_$t"*".vtu",cellfields=["uh"=>uh, "u_ex"=>u_ex(t)])
    end
end
