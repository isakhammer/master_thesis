from ngsolve import *
from ngsolve.meshes import MakeStructured2DMesh
import os
from datetime import datetime
import numpy as np
import sympy as sy


# -Δ^2*u + alpha*u = f in Ω,
# grad_u*n = g_1 and grad_Delta_u*n = g_2 on Γ
def man_sol(u_sy, x_sy, y_sy):

    # Symbolic differentiation
    Delta_u = sy.diff(u_sy, x_sy,x_sy) + sy.diff(u_sy, y_sy, y_sy)
    Delta2_u = sy.diff(Delta_u, x_sy,x_sy) + sy.diff(Delta_u, y_sy, y_sy)
    f_sy =  Delta2_u + alpha*u_sy
    grad_u_sy = (sy.diff(u_sy, x_sy), sy.diff(u_sy, y_sy))
    grad_Delta_u_sy = (sy.diff(Delta_u, x_sy), sy.diff(Delta_u, y_sy))

    # Convert problem to ngsolve coefficients
    u_ex = CoefficientFunction(eval(str(u_sy)))
    f = CoefficientFunction(eval(str(f_sy)))
    grad_u_ex = CoefficientFunction(eval(str(grad_u_sy)))
    grad_Delta_u_ex = CoefficientFunction(eval(str(grad_Delta_u_sy)))
    return u_ex, f, grad_u_ex, grad_Delta_u_ex


# Define Manufactured solution
x_sy, y_sy = sy.symbols('x y')

alpha = 1
(L,m,r) = (1,1,1)
u_sy = 100*sy.cos(x_sy * 2*sy.pi/L)*sy.cos(y_sy * 2*sy.pi/L)

# Transform to manufactured solution
u_ex, f, grad_u_ex, grad_Delta_u_ex = man_sol(u_sy, x_sy, y_sy)




def run(order, n_grid, vtk_dirname=None):  # Mesh related parameters

    mesh = MakeStructured2DMesh(quads=True, nx=n_grid,ny=n_grid)

    V = H1(mesh, order=order, dgjumps=True)

    u = V.TrialFunction()
    v = V.TestFunction()

    n = specialcf.normal(2)
    h = specialcf.mesh_size
    gamma = order*(order+1)

    grad_n      = lambda u: grad(u)*n

    hesse       = lambda u: u.Operator("hesse")
    hesse_nn    = lambda u: InnerProduct(n, hesse(u)*n)

    mean        = lambda u: 0.5*(u + u.Other())
    mean_n      = lambda u: 0.5*(grad(u) + grad(u.Other()))*n
    mean_nn     = lambda u: InnerProduct(n, 0.5*( hesse(u) )*n ) + InnerProduct(n, 0.5*( hesse(u.Other()) )*n )

    jump        = lambda u: u - u.Other()
    jump_n      = lambda u: (grad(u) - grad(u.Other()))*n
    jump_nn     = lambda u: InnerProduct(n,  (hesse(u) - hesse(u.Other) )*n)


    # Inner facets

    g_1 = grad_u_ex*n
    g_2 = grad_Delta_u_ex*n

    a = BilinearForm(V, symmetric=True)
    a += (InnerProduct(hesse(u), hesse(v)) + alpha*(u*v))*dx

    a += (-mean_nn(u)*jump_n(v) - mean_nn(v)*jump_n(v) )*dx(skeleton=True)
    a += (-hesse_nn(v)*grad_n(u) - hesse_nn(u)*grad_n(v) )*ds(skeleton=True)
    a += (gamma/h)*( jump_n(u)*jump_n(v) )*dx(skeleton=True)
    a += (gamma/h)*( grad_n(u)*grad_n(v) )*ds(skeleton=True)
    a.Assemble()

    l = LinearForm(V)
    l += (f*v) * dx
    l += ( -g_2*v + g_1*( -hesse_nn(v) + (gamma/h)*grad_n(v) ) )*ds(skeleton=True)
    l.Assemble()

    u_ex_h = GridFunction(H1(mesh, order=order+2))
    u_ex_h.Set(u_ex)

    u_h = GridFunction(V)
    u_h.vec.data = a.mat.Inverse() * l.vec

    e = u_h - u_ex
    de = grad(u_h) - grad(u_ex_h)
    el2 = sqrt(Integrate(e*e, mesh))
    eh1 = sqrt(Integrate(e*e + de*de, mesh))
    eh_energy = sqrt(Integrate( ( de*de )*dx, mesh) \
            + Integrate( ( ( h/gamma )*( de*n )*( de*n ) + gamma*h**(-1)*e*e )*ds(skeleton=True), mesh)
            # + Integrate( ( h*mean_n(e)*mean_n(e) +h**(-1)*jump(e)*jump(e) )*dx(skeleton=True), mesh) \
            )

    if vtk_dirname != None:
        filename=vtk_dirname+"/order_"+str(order)+"_n_"+str(n)
        vtk = VTKOutput(mesh,
                        coefs=[u_ex, u_h, e, de],
                        # coefs=[u_ex, u_h, e],
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


if __name__ == "__main__":

    dirname = "figures/biharmonic_CIP/"+datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    print("figures in ", dirname)
    os.makedirs(dirname, exist_ok=True)
    orders = [2, 3, 4]
    ns = [2**2, 2**3, 2**4, 2**5, 2**6, 2**7]

    convergence_analysis(orders, ns, dirname=dirname)

