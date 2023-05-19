##
using Gridap

# 
ε = 1/100


## Define geometry and mesh
N = 50
domain = (0,2π, 0, 2π)
partition = (N, N)
𝒯 = CartesianDiscreteModel(domain, partition) |> simplexify

vtkdirname = "cahn-hilliard"
mkpath(vtkdirname)
modelfname = "model"
writevtk(𝒯, joinpath(vtkdirname, modelfname))



# Spaces

V
