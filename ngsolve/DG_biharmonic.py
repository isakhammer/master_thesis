#!/usr/bin/env python

from ngsolve import *
# from ngsolve.webgui import Draw
from netgen.geom2d import unit_square

def poission_DG():
    # https://docu.ngsolve.org/v6.2.1802/i-tutorials/unit-2.8-DG/DG.html

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
    order=4
    fes = L2(mesh, order=order, dgjumps=True)
    u,v = fes.TnT()
    gfu = GridFunction(fes)  # solution

    jump_u = u-u.Other()
    jump_v = v-v.Other()
    n = specialcf.normal(2)

    def jumpdn(v,vhat):
        return n*(grad(v)-vhat)
    def hesse(v):
        return v.Operator("hesse")
    def hessenn(v):
        return InnerProduct(n, hesse(v)*n)

    mean_dudn = 0.5*n * (grad(u)+grad(u.Other()))
    mean_dvdn = 0.5*n * (grad(v)+grad(v.Other()))


    alpha = 4
    h = specialcf.mesh_size
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u)*grad(v))
    a += SymbolicBFI(alpha*order**2/h*jump_u*jump_v, skeleton=True)
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


if __name__ == "__main__":
    biharmonic_DG()


