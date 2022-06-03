#!/usr/bin/env python

from ngsolve import *
from netgen.geom2d import unit_square


def biharmonic_DG():
    mesh = Mesh (unit_square.GenerateMesh(maxh=1/10))
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


