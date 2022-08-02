using Gridap
import GridapMakie
import GLMakie

function main(;dirname::String ,n=10, L=1, m=1, r=1 )

    prename=dirname*"/L_"*string(round(L, digits=2 ))*"_m_"*string(m)*"_r_"*string(r)*"n_"*string(n)
    u(x) = cos(m*( 2π/L )*x[1])*cos(r*( 2π/L )*x[2])
    order = 2
    domain2D = (0, L, 0, L)
    partition2D = (n,n)
    model = CartesianDiscreteModel(domain2D,partition2D) |> simplexify


    # Spaces
    V = TestFESpace(model, ReferenceFE(lagrangian,Float64,order), conformity=:H1)
    U = TrialFESpace(V)
    Ω = Triangulation(model)
    Λ = SkeletonTriangulation(model)
    Γ = BoundaryTriangulation(model)


    fig = GLMakie.Figure(resolution = (1000, 1000))
    GLMakie.Axis(fig[1,1])
    GLMakie.plot!(Λ)
    GLMakie.wireframe!( Λ, color=:black, linewidth=2)
    GLMakie.wireframe!( Γ, color=:black, linewidth=2)
    GLMakie.save(prename*"_grid.png", fig)

    u_inter = interpolate(u, V)
    fig, _ , plt = GLMakie.plot(Ω, u_inter)
    GLMakie.Colorbar(fig[1,2], plt)
    GLMakie.save(prename*"_sol.png", fig)
end



dirname ="model"
if (isdir(dirname))
    rm(dirname, recursive=true)
end
mkdir(dirname)

main(dirname=dirname, n=2^4, L=1, m=1,r=1)
# main(dirname=dirname, n=2^4, L=2π, m=1,r=1 )
# main(dirname=dirname, n=2^4, L=1, m=7,r=3)

