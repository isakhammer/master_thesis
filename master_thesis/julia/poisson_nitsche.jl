##
using Gridap

# Define manufactured solution
u(x) = 100*cos(x[1])*cos(x[2])

# -Δu = f
f(x) = 2*u(x)
g(x) = u(x)

# ∂Δu/∂n
# q(x) = 0

function run_exp(n=8, k=2, use_quads=false)
    ## Mesh
    pmin = Point(0.,0.0)
    pmax = Point(2π, 2π)
    partition = (n, n)

    if !use_quads
        model = CartesianDiscreteModel(pmin, pmax, partition) |> simplexify
    else
        model = CartesianDiscreteModel(pmin, pmax, partition)
    end

    # Define triangulation
    Ω = Triangulation(model)
    Γ = BoundaryTriangulation(model)

    # Write out meshes
    folder = "figures/poisson_nitsche/"
    mkpath(folder)
    writevtk(model, folder*"model_n_$(n)")
    writevtk(Ω, folder*"mesh_n_$(n)")

    ## Function spaces
    order = k
    reffe = ReferenceFE(lagrangian, Float64, order)
    V = TestFESpace(Ω, reffe, conformity=:H1)
    U = TrialFESpace(V)

    ## Define the weak form
    degree = 2*order
    dΩ = Measure(Ω, degree)
    dΓ= Measure(Γ, degree)

    # Normal vectors
    n_Γ  = get_normal_vector(Γ)

    # Define mesh and stabilization parameters
    h = norm((pmax-pmin)./VectorValue(partition))

    # γ = 5.0*order*(order+1)  # Penalty parameter
    γ = 5  # Penalty parameter
    μ = γ/h

    a_Ω(u,v) =∫( ∇(v)⋅∇(u) )dΩ
    a_Γ(u,v) =∫( - ( ∇(u)⋅n_Γ )⊙v - u⊙( ∇(v)⋅n_Γ ) + μ*u⊙v )dΓ
    a(u,v) = a_Ω(u,v) + a_Γ(u,v)

    l_Ω(v) = ∫( v⊙f )dΩ
    l_Γ(v) = ∫( -(( ∇(v)⋅n_Γ )⊙g) + μ*(g⊙v) )dΓ
    l(v) = l_Ω(v) + l_Γ(v)

    ## Solve and postprocess
    op = AffineFEOperator(a, l, U, V)
    uh = solve(op)

    # Compute interpolation of exact solution for visual comparison
    u_inter = interpolate_everywhere(u, V)

    # Compute error
    e = u - uh

    outputfile = folder*"possion_nitsche_n_$(n)_k_$(k)_use_quads_$(use_quads)"
    writevtk(Ω, outputfile, cellfields=["uh" => uh, "u_inter" => u_inter, "error" => e])

    l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
    h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))
    el2 = l2(e)
    eh1 = h1(e)

    println("h = $h, L2 error = $el2")
    println("h = $h, H1 error = $eh1")

    (h, el2, eh1)
end

function conv_test(ns, k, use_quads)

    hs = Float64[]
    el2s = Float64[]
    eh1s = Float64[]

    for n in ns
        h, el2, eh1 = run_exp(n, k, use_quads)
        push!(hs, h)
        push!(el2s, el2)
        push!(eh1s, eh1)
    end

    println("Mesh sizes = $hs")
    println("L2 errors  = $el2s")
    println("H1 errors  = $eh1s")

    # Compute eoc
    (hs, el2s, eh1s)
end

function compute_eoc(hs, errs)
    eoc = log.(errs[1:end-1]./errs[2:end])./log.(hs[1:end-1]./hs[2:end])
end

function main()
    use_quads = false
    ks  = [2, 3, 4]
    ns = [4, 8, 16, 32, 64, 128]
    for k in ks
        println("Running EOC test for k = $k")
        println("========================================")
        hs, el2s, eh1s = conv_test(ns, k, use_quads)
        eoc_l2 = compute_eoc(hs, el2s)
        eoc_h1 = compute_eoc(hs, eh1s)
        println("L2 EOC = $eoc_l2")
        println("H1 EOC = $eoc_h1")
        println("========================================")
    end
end

main()
