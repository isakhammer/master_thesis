#!/usr/bin/env python

from ngsolve import *
# from ngsolve.webgui import Draw
from netgen.geom2d import unit_square

def poission_DG():
    # https://docu.ngsolve.org/v6.2.1808/i-tutorials/unit-2.8-DG/DG.html
    mesh = Mesh (unit_square.GenerateMesh(maxh=0.1))
    order=4
    fes = L2(mesh, order=order, dgjumps=True)
    u,v = fes.TnT()
    gfu = GridFunction(fes)  # solution

    jump_u = u-u.Other()
    jump_v = v-v.Other()
    n = specialcf.normal(2)
    mean_dudn = 0.5*n * (grad(u)+grad(u.Other()))
    mean_dvdn = 0.5*n * (grad(v)+grad(v.Other()))

    alpha = 4
    h = specialcf.mesh_size
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u)*grad(v))
    a += SymbolicBFI(alpha*order**2/h*u*v, BND, skeleton=True)
    a += SymbolicBFI(-mean_dudn*jump_v -mean_dvdn*jump_u, skeleton=True)
    a += SymbolicBFI(-n*grad(u)*v-n*grad(v)*u, BND, skeleton=True)
    a.Assemble()

    f = LinearForm(fes)
    f += SymbolicLFI(1*v)
    f.Assemble()

    print(f.vec)
    gfu = GridFunction(fes, name="uDG")
    gfu.vec.data = a.mat.Inverse() * f.vec
    Draw (gfu)

def biharmonic_DG():

    mesh = Mesh (unit_square.GenerateMesh(maxh=0.1))
    order = 4
    fes = L2(mesh, order=order, dgjumps=True)
    wh,vh = fes.TnT()
    gfu = GridFunction(fes)  # solution

    alpha = 4
    gamma = 3
    n = specialcf.normal(2)
    h = specialcf.mesh_size

    def hesse(v):
        return v.Operator("hesse")
    def hesse_nn(v):
        return InnerProduct(n, hesse(v)*n)
    def mean_nn(v):
        return 0.5*n * (hessenn(v)+hessenn(v.Other()))
    def jump_n(v):
        return n*(grad(v)-vhat)

    dS = dx(element_boundary=True)
    A = BilinearForm(fes)
    A += ( alpha*u*v )*dx \
    A += InnerProduct(hesse(u),hesse(v)), VOL
    A += SymbolicBFI( mean_nn(wh)*jump_n(vh), VOL )
    A += SymbolicBFI(( mean_nn(vh)*jump_n(wh) ), VOL)
    A += SymbolicBFI( ( gamma/h )*jump_n(wh)*jump_n(vh), VOL)
    A.Assemble()

    F = LinearForm(fes)
    F += SymbolicBFI(f*v)
    F += SymbolicBFI(g_2*vh, BND) + SymbolicBFI(n*g_2*vh, BND) +
    F.Assemble()

    gfu = GridFunction(fes, name="uDG")
    gfu.vec.data = a.mat.Inverse() * f.vec
    Draw (gfu)


if __name__ == "__main__":
    biharmonic_DG()
    # poission_DG()


