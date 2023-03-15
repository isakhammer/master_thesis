from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh

from netgen.geom2d import SplineGeometry
import numpy as np
import pandas as pd
import os
from datetime import datetime
from matplotlib import pyplot as plt
import sympy as sy


# -Δu = f in Ω, and u = g on Γ
def man_sol(u_sy, x_sy, y_sy):

    # Symbolic differentiation
    Delta_u = sy.diff(u_sy, x_sy,x_sy) + sy.diff(u_sy, y_sy, y_sy)
    f_sy = - Delta_u

    # Convert problem to ngsolve coefficients
    u_ex = CoefficientFunction(eval(str(u_sy)))
    f = CoefficientFunction(eval(str(f_sy)))
    g = u_ex
    return u_ex, f, g


# Define Manufactured solution
x_sy, y_sy = sy.symbols('x y')
(L,m,r) = (1,1,1)
u_sy = 100*sy.cos(x_sy * 2*sy.pi/L)*sy.sin(y_sy * 2*sy.pi/L)
# u_sy = sy.cos(x_sy)*sy.cos(y_sy)
# u_sy = x_sy**2 + y_sy**2

# Transform to manufactured solution
u_ex, f, g = man_sol(u_sy, x_sy, y_sy)




def run(order, n, vtk_dirname=None):  # Mesh related parameters

    mesh = MakeStructured2DMesh(quads=True, nx=n,ny=n)

    V = H1(mesh, order=order)
    u = V.TrialFunction()
    v = V.TestFunction()

    n_Gamma = specialcf.normal(2)
    h = specialcf.mesh_size
    gamma = order*(order+1)


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
        order_ex = order + 2
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
    el2s = np.array(el2s)
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

    dirname = "figures/"+datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    print("figures in ", dirname)
    os.makedirs(dirname, exist_ok=True)
    orders = [1, 2, 3, 4]
    ns = [2**2, 2**3, 2**4, 2**5, 2**6, 2**7]

    convergence_analysis(orders, ns, dirname)

