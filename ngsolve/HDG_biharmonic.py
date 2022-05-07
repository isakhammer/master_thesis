#!/usr/bin/env python

from ngsolve import *
# from ngsolve.webgui import Draw
from netgen.geom2d import unit_square

def biharmonic_HDG():
    # 2.9 Fourth order equations

    # We consider the Kirchhoff plate equation: Find $w \in H^2$, such that
    # $$
    # \int \nabla^2 w : \nabla^2 v = \int f v
    # $$

    # A conforming method requires $C^1$ continuous finite elements. But there is no good option available, and thus there is no $H^2$ conforming finite element space in NGSolve.

    # We have the following two alternatives:

    ## Hybridized $C^0$-continuous interior penalty method:

    # A simple way out is to use continuous elements, and treat the missing $C^1$-continuity by a Discontinuous Galerkin method. A DG formulation is

    # $$
    # \sum_T \nabla^2 w : \nabla^2 v
    # - \int_{E} \{\nabla^2 w\}_{nn} \, [\partial_n v]
    # - \int_{E} \{\nabla^2 v\}_{nn} \, [\partial_n w] + \alpha \int_E  [\partial_n w]  [\partial_n v]
    # $$

    # [Baker 77, Brenner Gudi Sung, 2010]

    # We consider its hybrid DG version, where the normal derivative is a new, facet-based variable:

    # $$
    # \sum_T \nabla^2 w : \nabla^2 v
    # - \int_{\partial T} (\nabla^2 w)_{nn} \, (\partial_n v - \widehat{v_n})
    # - \int_{\partial T} (\nabla^2 v)_{nn} \, (\partial_n w - \widehat{w_n}) + \alpha \int_E (\partial_n v - \widehat{v_n}) (\partial_n w - \widehat{w_n})
    # $$

    # The facet variable is the normal derivative $n_E \cdot \nabla w$, what is oriented along the arbitrarily chosen edge normal-vector.
    # We cannot use the FacetSpace since it does not have the orientation, but we can use the normal traces of an HDiv space.
    # We don't need inner basis functions, so we set order inner to 0:

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

    Draw(u.components[0], mesh, name="disp_DG")
    Draw(grad (u.components[0]), mesh, "grad")
    Draw(hesse (u.components[0]), mesh, "hesse")

if __name__ == "__main__":
    biharmonic_HDG()
