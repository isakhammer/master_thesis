from ngsolve import *
from ngsolve.webgui import Draw
from netgen.geom2d import unit_square
mesh = Mesh (unit_square.GenerateMesh(maxh=0.1))
order = 3

V1 = H1(mesh, order=order, dirichlet="left|bottom|right|top")
V2 = NormalFacetFESpace(mesh, order=order-1, dirichlet="left|bottom|right|top")
V = V1*V2

w,what = V.TrialFunction()
v,vhat = V.TestFunction()

n = specialcf.normal(2)
h = specialcf.mesh_size

def jumpdn(v,vhat):
    return n*(grad(v)-vhat)
def hesse(v):
    return v.Operator("hesse")
def hessenn(v):
    return InnerProduct(n, hesse(v)*n)

dS = dx(element_boundary=True)
a = BilinearForm(V)
a += InnerProduct (hesse(w), hesse(v)) * dx \
     - hessenn(w) * jumpdn(v,vhat) * dS \
     - hessenn(v) * jumpdn(w,what) * dS \
     + 3*order*order/h * jumpdn(w,what) * jumpdn(v,vhat) * dS
a.Assemble()

f = LinearForm(V)
f += 1*v*dx
f.Assemble()


u = GridFunction(V)
u.vec.data = a.mat.Inverse(V.FreeDofs()) * f.vec

Draw (u.components[0], mesh, "disp_DG")
Draw (grad (u.components[0]), mesh, "grad")
Draw (hesse (u.components[0]), mesh, "hesse")


