from ngsolve import *
from ngsolve.webgui import Draw
from netgen.geom2d import SplineGeometry
import numpy as np
import pandas as pd
import os
from datetime import datetime
from matplotlib import pyplot as plt

L, m, r = (1,1,1)
u_ex = CoefficientFunction(100*cos(m*( 2*np.pi/L )*x)*cos(r*( 2*np.pi/L )*y))


def run(order, n):  # Mesh related parameters
    mesh = Mesh(unit_square.GenerateMesh(maxh=(1/n)))

    V = H1(mesh, order=order)
    u = V.TrialFunction()
    v = V.TestFunction()

    n = specialcf.normal(2)
    h = specialcf.mesh_size
    penalty = 3*order**2


    a = BilinearForm(V, symmetric=True)
    a += grad(u)*grad(v)*dx
    a += (-n*grad(u)*v - n*grad(v)*u + penalty/h*u*v)*ds(skeleton=True)
    a.Assemble()

    f = LinearForm(V)
    f += 1 * v * dx
    f += ( -n*grad(v)*u_ex + penalty/h*u_ex*v)*ds(skeleton=True)
    f.Assemble()

    u_h = GridFunction(V)
    u_h.vec.data = a.mat.Inverse() * f.vec

    def interpolate(u_h, u_ex):
        order_ex = 4
        Vh_ex = H1(mesh, order=order_ex)
        u_ex_inter = GridFunction(Vh_ex)
        u_ex_inter.Set(u_ex)

        u_h_inter = GridFunction(Vh_ex)
        u_h_inter.Set(u_h)
        return u_ex_inter, u_h_inter

    u_ex_inter, u_h_inter = interpolate(u_h, u_ex)
    e = u_h_inter - u_ex_inter
    de = u_h_inter.Deriv() - u_ex_inter.Deriv()

    el2 = sqrt(Integrate(e*e, mesh))
    eh1 = sqrt(Integrate(de*de, mesh))
    eh_energy = sqrt(Integrate(e*e + de*de, mesh))
    return el2, eh1, eh_energy


def print_results(ns, el2s, eh1s, ehs_energy, order):
    # Compute convergence rates

    ns = np.array(ns)
    hs = 1/ns
    el2s = np.array(eh1s)
    eh1s = np.array(eh1s)
    ehs_energy = np.array(ehs_energy)

    def compute_eoc(hs, errs):
        eoc = np.log(errs[:-1]/errs[1:])/np.log(hs[:-1]/hs[1:])
        return eoc

    eoc_el2 = compute_eoc(hs, el2s)
    eoc_eh1 = compute_eoc(hs, eh1s)
    eoc_eh_energy = compute_eoc(hs, eh1s)

    print("\n==============")
    print("SUMMARY")
    print("Order =", order," Mesh sizes = ", hs)
    print("L2 errors = ", el2s)
    print("H1 errors = ", eh1s)
    print("Energy errors  = ", eh1s)
    print("EOC L2 = ", eoc_el2)
    print("EOC H1 = ", eoc_eh1)
    print("EOC Energy = ", eoc_eh_energy)
    print("==============\n")



def convergence_analysis(orders, ns, dirname):

    for order in orders:
        el2s = []
        eh1s = []
        ehs_energy = []
        for n in ns:
            (el2, eh1, eh_energy) = run(order, n)
            el2s.append(el2)
            eh1s.append(eh1)
            ehs_energy.append(eh_energy)
        print_results(ns, el2s, eh1s, ehs_energy, order)

    return 0



if __name__ == "__main__":

    dirname = "figures/"+datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")
    print("figures in ", dirname)
    os.makedirs(dirname, exist_ok=True)
    orders = [2, 3, 4]
    ns = [2**2, 2**3, 2**4, 2**5, 2**6]

    convergence_analysis(orders, ns, dirname)

