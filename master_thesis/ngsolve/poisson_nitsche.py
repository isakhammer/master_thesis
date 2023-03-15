from ngsolve import *
from ngsolve.webgui import Draw
from netgen.geom2d import SplineGeometry
import numpy as np
import pandas as pd
import os
from datetime import datetime
from matplotlib import pyplot as plt

# -Δu = f in Ω, and u = g on Γ

# Example 1
L, m, r = (1,1,1)
u_ex = CoefficientFunction(100*cos(m*( 2*np.pi/L )*x)*cos(r*( 2*np.pi/L )*y))
f = 2*u_ex
g = u_ex

# Example 2
u_ex = CoefficientFunction(x**2*y**2)
f = -2*(y**2 + x**2)
g = u_ex

def run(order, n, vtk_dirname=None):  # Mesh related parameters
    mesh = Mesh(unit_square.GenerateMesh(maxh=(1/n)))

    V = H1(mesh, order=order)
    u = V.TrialFunction()
    v = V.TestFunction()

    n_Gamma = specialcf.normal(2)
    h = specialcf.mesh_size
    gamma = 5

    a = BilinearForm(V, symmetric=True)
    a += grad(u)*grad(v)*dx
    a += (-n_Gamma*grad(u)*v - n_Gamma*grad(v)*u + (gamma/h)*u*v)*ds(skeleton=True)
    a.Assemble()

    l = LinearForm(V)
    l += f * v * dx
    l += ( -n_Gamma*grad(v)*g + (gamma/h)*g*v)*ds(skeleton=True)
    l.Assemble()

    u_h = GridFunction(V)
    u_h.vec.data = a.mat.Inverse() * l.vec

    def interpolate(u_ex, u_h):
        order_ex = order
        Vh_ex = H1(mesh, order=order_ex)
        u_ex_inter = GridFunction(Vh_ex)
        u_ex_inter.Set(u_ex)

        u_h_inter = GridFunction(Vh_ex)
        u_h_inter.Set(u_h)
        return u_ex_inter, u_h_inter

    u_ex_inter, u_h_inter = interpolate(u_ex, u_h)
    e = u_h_inter - u_ex_inter
    de = u_h_inter.Deriv() - u_ex_inter.Deriv()

    el2 = sqrt(Integrate(e*e, mesh))
    eh1 = sqrt(Integrate(de*de, mesh))
    eh_energy = sqrt(Integrate(e*e + de*de, mesh))

    if vtk_dirname != None:
        filename=vtk_dirname+"/order_"+str(order)+"_n_"+str(n)
        print("vtk in ", filename)
        vtk = VTKOutput(mesh,
                        coefs=[u_ex, u_h, e, de],
                        names=["u_ex", "u_h", "e", "de"],
                        filename=filename,
                        subdivision=3)
        vtk.Do()

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
            (el2, eh1, eh_energy) = run(order, n, vtk_dirname=dirname)
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

