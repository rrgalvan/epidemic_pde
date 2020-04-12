from fenics import *
from mshr import *
import matplotlib.pyplot as plt

test_name = 'SIR_react_difuss'

# Deafault biological parameters [1]
default_parameters = {
    # Numerical method parameters
    "nx": None,
    "dt": None,
    "n_iter": None,
    "theta": 0.5, # Crank-Nicoloson constant
    # Biological method parameters
    "mu": 0.04,       # year^-1
    "nu": 24,         # year^-1
    "N" : 3e+5,       # total number of individuals
    "alpha" : 0.01,   # diffusion coeff.
    "beta": 4.4333e-4 # infection rate
}

def plot_postprocess(plot_options, filename="out.png"):
    print(f"  plotting {plot_options}, {filename}")
    if "show" in plot_options and plot_options["show"]:
        plt.show()
    if "save" in plot_options and plot_options["save"]:
        plt.savefig ( filename )
    plt.close()

def run_test (
    parameters = default_parameters,
    plot_options={}
    ):
    param = parameters
    print(param)
    nx = param["nx"]
    dt = param["dt"]
    n_iter = param["n_iter"]
    theta = param["theta"]

#
    # Mesh
    #

    mesh = UnitIntervalMesh ( nx )
    if plot_options:
        plot(mesh, title = f'{test_name} Mesh')
        plot_postprocess(plot_options, filename = f'{test_name}_mesh.png')
    #
    # Finite elements
    #
    Vh = FunctionSpace( mesh, "Lagrange", 1 )
    S = TrialFunction( Vh )
    I = TrialFunction( Vh )
    SS = TestFunction( Vh )
    II = TestFunction( Vh )

    #
    # Time discrtization
    #
    t = 0
    Tmax = n_iter*dt

    #
    # Biological parameters<
    #
    mu = Constant( param["mu"] )
    nu = Constant( param["nu"] )
    N  = Constant( param["N"] )
    alpha = Constant( param["alpha"] )
    beta = Constant( param["beta"] )
    print("---------------nu =", param["nu"])

    #
    # Initial data
    #
    S_init = Expression ( "C/2-C*abs(x[0]-0.5)", C=325000, degree = 1 )
    S0 = interpolate(S_init, Vh)
    I_init = Expression ( "C/2-C*abs(x[0]-0.5)", C=7500, degree = 1 )
    I0 = interpolate(I_init, Vh)
    yield S0.vector(), I0.vector()
    if plot_options:
        plot(S0, title = f'{test_name} initial S')
        plot_postprocess(plot_options, filename = f'{test_name}_S0.png')
        plot(I0, title = f'{test_name} initial I')
        plot_postprocess(plot_options, filename = f'{test_name}_I0.png')

    #
    # Variational formulation
    #

    # S: bilinear form
    a_S = (1 + dt*mu) * S*SS*dx +\
     dt*alpha*theta * dot( grad(S), grad(SS) ) * dx
    a_S += dt*beta * I0 * S*SS * dx # Nonlinear term
    # S: linear form
    b_S = S0*SS * dx + dt*mu*N*SS *dx +\
     dt*alpha*(1-theta) * dot( grad(S0), grad(SS) ) * dx

    # I: bilinear form
    a_I = (1+dt*(mu+nu)) * I*II * dx +\
     dt*alpha*theta * dot( grad(I), grad(II) ) * dx
    # I: linear form
    b_I = I0*II * dx + dt*beta*S0*I0*II * dx +\
     dt*alpha*(1-theta) * dot( grad(I0), grad(II) ) * dx

    #
    # Time loop
    #

    S1 = Function(Vh)
    I1 = Function(Vh)

    for iter in range(1,n_iter+1):
        t = t+dt
        print(f'Time setep {iter} (t={t:.3})')

        solve ( a_S == b_S, S1 )
        solve ( a_I == b_I, I1 )

        #
        # Print info
        #
        S_max, S_min = max(S1.vector()), min(S1.vector())
        I_max, I_min = max(I1.vector()), min(I1.vector())
        print(f"  máx(S)={S_max:.5e}, min(S)={S_min:.5e}")
        print(f"  máx(I)={I_max:.5e}, min(I)={I_min:.5e}")
        integral_S = assemble(S1*dx); print(f"  int(S)={integral_S:.5e}" )
        integral_I = assemble(I1*dx); print(f"  int(I)={integral_I:.5e}" )

        #
        # Plot
        #
        if plot_options and iter%1 == 0:
            plot(S1, title = f'{test_name} S, iter {iter}')
            plot_postprocess(plot_options, filename = f'{test_name}_S1_{iter}.png')
            plot(I1, title = f'{test_name} I, iter {iter}')
            plot_postprocess(plot_options, filename = f'{test_name}_I1_{iter}.png')

        #
        # Prepare next iteration
        #
        S0.assign(S1)
        I0.assign(I1)

        yield S1.vector(), I1.vector()
