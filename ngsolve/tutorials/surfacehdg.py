from ngsolve import *
import netgen.gui
from math import pi

from netgen.csg import *
geo = CSGeometry()
geo.Add(Sphere(Pnt(0,0,0),1).bc("sphere"))

from netgen.meshing import MeshingParameters
from netgen.meshing import MeshingStep
mp = MeshingParameters(maxh=0.5, perfstepsend = MeshingStep.MESHSURFACE)

mesh = Mesh(geo.GenerateMesh(mp=mp))
order = 4
if order > 0:
    mesh.Curve(order)
Draw(mesh)

t = Parameter(0.)
b = CoefficientFunction((y,-x,0))
Draw (b, mesh, "b")
eps = 5e-5

n = specialcf.normal(mesh.dim)
h = specialcf.mesh_size
local_tang = specialcf.tangential(mesh.dim)
con = Cross(n,local_tang) #co-normal pointing outside of a surface element
bn = InnerProduct(b,con)


fesl2 = SurfaceL2(mesh, order=order)
fesedge = FacetSurface(mesh, order=order)

fes = FESpace([fesl2,fesedge])
u,uE = fes.TrialFunction()
v,vE = fes.TestFunction()

gradu = u.Trace().Deriv()
gradv = v.Trace().Deriv()
jumpu = u - uE
jumpv = v - vE
dudn = InnerProduct(gradu,con)
dvdn = InnerProduct(gradv,con)


a = BilinearForm(fes)
# diffusion

dr = ds(element_boundary=True)

alpha = 10*(order+1)**2
a += eps*gradu * gradv * ds
a += - eps*dudn * jumpv * dr
a += - eps*dvdn * jumpu * dr
a += eps*alpha/h * jumpu * jumpv * dr


# CONVECTION FORMULATION
# a = BilinearForm(fes)
# # convection (surface version of Egger+Schöberl formulation):
# a += - u * InnerProduct(b,gradv) * ds
# a += IfPos(bn, bn*u, bn*uE) * v * dr + IfPos(bn, bn, 0) * (uE - u) * vE * dr
# a.Assemble()
m = BilinearForm(fes)
m += u*v*ds
m.Assemble()

dt = 0.02
mstar = m.mat.CreateMatrix()
mstar.AsVector().data = m.mat.AsVector() + dt * a.mat.AsVector()
mstarinv = mstar.Inverse()

from math import pi
T = 8*pi


gfu = GridFunction(fes)
gfu.components[0].Set(1.5*exp(-20*(x*x+(y-1)**2+(z)**2)),definedon=mesh.Boundaries("sphere"))
Draw (gfu.components[0], mesh, "u")
Draw (gfu.components[0]*n, mesh, "un", sd=4, autoscale=False)
res = gfu.vec.CreateVector()
SetVisualization(deformation=True)
t.Set(0)

for i in range(int(T/dt)):
    print ("t=", dt*i, end="\r")
    res.data = m.mat * gfu.vec
    gfu.vec.data = mstarinv * res
    Redraw(blocking=True)
    t.Set(t.Get()+dt)

