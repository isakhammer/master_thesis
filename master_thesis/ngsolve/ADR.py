from ngsolve import *
from ngsolve.webgui import Draw
from netgen.geom2d import SplineGeometry
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


jump = lambda v : v-v.Other()
avg  = lambda v : 0.5*(v+v.Other())
flux = lambda v, n : 0.5*(Grad(v)+Grad(v).Other())*n

def ADRSolver(epsilon, b, c, eta, u_D, f, u_ex, errors, order_u=1,
                 init_maxh=1.0, num_ref=0):  # Mesh related parameters

    # Generating mesh
    mesh = Mesh(unit_square.GenerateMesh(maxh=init_maxh))
    for ref in range(num_ref):
          mesh.Refine()

    # Function space (discontinuous)
    V = L2(mesh, order=order_u, dgjumps=True)
    u,v = V.TnT() # Trial and test function

    # Parameters in bilinear and linear form
    gamma = Parameter(10)*order_u**2
    gamma_bnd = Parameter(10)*order_u**2
    n = specialcf.normal(2)
    h = specialcf.mesh_size

    ## Set up bilinear form:
    a = BilinearForm(V)
    # Symmetric interior point formulation for diffusion part
    aSIP = grad(u)*grad(v) * dx \
            - (n*grad(u)*v + n*grad(v)*u) * ds(skeleton=True) \
            - (flux(u,n)*jump(v) + flux(v,n)*jump(u)) * dx(skeleton=True) \
            + gamma * (1/h*u*v) * ds(skeleton=True) \
            + gamma * (1/h*jump(u)*jump(v)) * dx(skeleton=True)
    a += epsilon * aSIP
    # Central flux formulation for advection-reaction part
    aCF = (b*grad(u)*v) * dx + (c*u*v) * dx \
            - (b*n*jump(u)*avg(v)) * dx(skeleton=True) \
            - (b*n*u*v) * ds(definedon=mesh.Boundaries("bottom"))
    a += aCF
    # Upwind term for recovered convergence rate
    absbn = IfPos(b*n, b*n, -b*n)
    upwind = eta/2*(absbn*jump(u)*jump(v)) * dx(skeleton=True)
    a += upwind
    a.Assemble()

    ## Set up linear form / rhs:
    l = LinearForm(V)
    # SIP:
    l += f*v * dx
    l += - epsilon * (n*grad(v)*u_D) * ds(skeleton=True)
    l += epsilon*gamma_bnd * (1/h*u_D*v) * ds(skeleton=True)
    # CF:
    l += -(b*n*u_D*v) * ds(definedon=mesh.Boundaries("bottom"))
    l.Assemble()

    # Solve system
    u_h = GridFunction(V)
    u_h.vec.data = a.mat.Inverse(freedofs= V.FreeDofs())*l.vec

    ## Calculate error norms
    e = u_h - u_ex

    # Interpolate u_ex and u_h into higher order space so we can compute derivative
    order_ex = order_u+2
    Vh_ex = H1(mesh, order=order_ex)
    u_ex_h = GridFunction(Vh_ex)
    u_ex_h.Set(u_ex)
    u_h_inter = GridFunction(Vh_ex)
    u_h_inter.Set(u_h)
    de = u_h_inter.Deriv() - u_ex_h.Deriv()

    errors["eD"].append(epsilon * sqrt(Integrate(de*de, mesh)))
    errors["eR"].append(sqrt(Integrate(c**2*e*e, mesh)))

    # https://ngsolve.org/forum/ngspy-forum/1349-how-to-compute-the-sum-of-jump-errors-on-edges-in-dg
    errors["eF"].append(sqrt(Integrate(1/h*jump(e)*jump(e)*dx(element_boundary=True), mesh)))
    errors["eAd1"].append(sqrt(Integrate(absbn*jump(e)*jump(e)*dx(element_boundary=True), mesh)))

    bc = 1 # L infinity norm of b
    errors["eAd2"].append(sqrt(Integrate(h/bc*(b*de)*(b*de), mesh)))

    return u_h, errors

def getTestProblem(epsilon):
    # Set up test problem from given manufactured solution.
    b = CoefficientFunction((0, 1))
    c = CoefficientFunction(0.1)

    u = CoefficientFunction(cos(pi*x) * (1-exp((y-1)/epsilon)) / (1-exp(-2/epsilon)) \
                               + 0.5*cos(pi*x)*sin(pi*y))

    # b•grad(u) in our case is just the y derivative:
    dy_u = CoefficientFunction(-cos(pi*x)/epsilon * exp((y-1)/epsilon)/(1-exp(-2/epsilon)) \
                                   + 0.5*pi*cos(pi*x)*cos(pi*y))

    laplace_u = CoefficientFunction(-pi**2 * cos(pi*x) * (1-exp((y-1)/epsilon))/(1-exp(-2/epsilon))\
                                       - cos(pi*x)/epsilon**2 * exp((y-1)/epsilon)/(1-exp(-2/epsilon))\
                                       - pi**2 * cos(pi*x) * sin(pi*y))
    f = -epsilon*laplace_u + dy_u + c*u
    return b, c, u, f



# Compute convergence rates
compute_eoc = lambda errors : np.insert( (np.log(errors[:-1]) - np.log(errors[1:]))/np.log(2), 0, np.Inf)

def resultsTable(errors):
    eoc1 = compute_eoc(errors["eD"])
    eoc2 = compute_eoc(errors["eR"])
    eoc3 = compute_eoc(errors["eF"])
    eoc4 = compute_eoc(errors["eAd1"])
    eoc5 = compute_eoc(errors["eAd2"])

    table = pd.DataFrame({'Error 1': errors["eD"], 'EOC 1': eoc1,
                          'Error 2': errors["eR"], 'EOC 2': eoc2,
                          'Error 3': errors["eF"], 'EOC 3': eoc3,
                          'Error 4': errors["eAd1"], 'EOC 4': eoc4,
                          'Error 5': errors["eAd2"], 'EOC 5': eoc5})
    return table

def plotErrors(errors, epsilon):
    plt.loglog(errors["mesh_sizes"],errors["eD"], '.--', label=r"$\epsilon||\nabla e_k||_{\Omega}$")
    plt.loglog(errors["mesh_sizes"],errors["eR"], '.--', label=r"$c_o||e_k||_\Omega$")
    plt.loglog(errors["mesh_sizes"],errors["eF"], '.--', label=r"$||h^{-1/2}[e_k]||_{\mathcal{F}_h}$")
    plt.loglog(errors["mesh_sizes"],errors["eAd1"], '.--', label=r"$|| |b \cdot n_F|^{1/2}[e_k]||_{\mathcal{F}_h}$")
    plt.loglog(errors["mesh_sizes"],errors["eAd2"], '.--', label=r"$|| (h/b_c)^{1/2} b \cdot \nabla e_k||_\Omega$")
    plt.xlabel(r"Mesh size $h_k$")
    plt.ylabel(r"Contribution $E_k$")
    plt.title(r"Error norms for numerical solution when $\epsilon=%.1e$"%epsilon)
    plt.grid()
    plt.legend()
    plt.show()

def runConvergenceTest(epsilon, eta=1, max_num_ref=4, order_u=2, draw=False):
    # Set eta=0 to remove upwinding
    # To collect errors
    errors = {
        "eD" : [],
        "eR" : [],
        "eF" : [],
        "eAd1" : [],
        "eAd2" : [],
    }

    solver_prms = {
        "order_u": order_u,
        "init_maxh": 0.2,
        "num_ref" : 0,
    }

    for num_ref in range(max_num_ref+1):
        b, c, u_ex, f = getTestProblem(epsilon)
        u_h, errors = ADRSolver(epsilon, b, c, eta, u_ex, f, u_ex, errors, **solver_prms)

        # Plot solution
        if draw:
            mesh = u_h.space.mesh
            Draw(u_h, mesh, "u_h")

        # Increase number of refinements
        solver_prms["num_ref"] += 1

    errors["mesh_sizes"] = [solver_prms["init_maxh"] * 2**(-i) for i in range(0,max_num_ref+1)]

    print(resultsTable(errors))

    plotErrors(errors, epsilon)

if __name__ == "__main__":
    runConvergenceTest(epsilon=1.0, eta=1, max_num_ref=5)

