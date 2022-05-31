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

    mesh = Mesh (unit_square.GenerateMesh(maxh=1/6))
    order = 4
    fes = L2(mesh, order=order, dgjumps=True)
    wh,vh = fes.TnT()
    gfu = GridFunction(fes)  # solution

    g = 1.0
    f = 1.0
    alpha = 0.1 # has to be strictly above zero
    gamma = 0.1
    n = specialcf.normal(2)
    h = specialcf.mesh_size

    def hesse(v):
        return v.Operator("hesse")
    def hesse_nn(v):
        return InnerProduct(n,  hesse(v)*n)
    def mean_nn(v):
        return 0.5* (hesse_nn(v)+hesse_nn(v.Other()))
    def jump_n(v):
        return n*(grad(v)-grad(v.Other()))

    dS = dx(element_boundary=False)
    A = BilinearForm(fes, symmetric=True)


    # Alternative 1
    # A += ( alpha*wh*vh )*dx
    # A += InnerProduct(hesse(wh),hesse(vh))*dx

    # Alternative 2 (gives same results as Alternative 1)
    A += SymbolicBFI( alpha*wh*vh +  InnerProduct( hesse(wh), hesse(vh)))

    # Boundary terms
    A += SymbolicBFI( mean_nn(wh)*jump_n(vh),  VOL, skeleton=True )
    A += SymbolicBFI( mean_nn(vh)*jump_n(wh),  VOL, skeleton=True )
    A += SymbolicBFI( ( gamma/h )*jump_n(wh)*jump_n(vh),  VOL,  skeleton=True)

    A.Assemble()

    F = LinearForm(fes)
    F += SymbolicLFI(f*vh)
    F += SymbolicLFI(g*vh, BND, skeleton=True)
    F.Assemble()

    gfu = GridFunction(fes, name="uDG")
    gfu.vec.data = A.mat.Inverse() * F.vec
    Draw (gfu)


if __name__ == "__main__":
    biharmonic_DG()
    # poission_DG()


