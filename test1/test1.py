#
# Resolución numérica del modelo SI+reacción+difusión
#
#   S_t - alpha \Delta S = mu N - mu S - beta S I
#   I_t - alpha \Delta I = -(mu + nu) I  + beta S I
#   + homogeneous b.c.
# Here:
#   S, I: susceptible and infectious individuals,
#   N: total population (fixed positive constant)
# From this model, one can compute R (recovered individuals) as
#   R_t - alpha I = nu I + mu R
#
# Test basado en:
#   [1] Numerical modelling of an SIR epidemic model with diffusion
#   Settapat Chinviriyasit, Wirawan Chinviriyasit
#   Applied Mathematics and Computation 216 (2010) 395–409
#
import numpy as np
import matplotlib.pylab as plt

import sys
sys.path.append("../") # El módulo está en el directorio superior
import schemes.SIR_react_diffusion_1d as SIR_scheme
from schemes.SIR_react_diffusion_1d import test_name, run_test

#plt.style.use('seaborn-pastel')
plt.style.use('seaborn-dark')
plt.rcParams["figure.dpi"] = 300

def test_1 ( ):
    from copy import copy

    parameters = copy(SIR_scheme.default_parameters)
    parameters["nx"] = 10
    parameters["dt"] = 1.e-2
    parameters["n_iter"] = 30
    parameters["beta"] = 8.e-4 # Infection rate
    parameters["mu"] = 0 # Nacimientos

    import time
    print ( f'{test_name}:' )
    print ( '  FENICS/Python version' )
    print ( '  Implicit solver.' )
    print ( ' ', time.ctime ( time.time() ) )

    plot_all_iterations = False
    S_total = [] # Total number of susceptible at each time step
    I_total = [] # Total number of infected at each time step
    R_total = [] # Total number of recovered at each time step
    for S, I, R in run_test ( parameters, plot_options={ "save": False } ):
        S_array = S.get_local(); S_total.append(np.sum(S_array))
        I_array = I.get_local(); I_total.append(np.sum(I_array))
        R_array = R.get_local(); R_total.append(np.sum(R_array))
        print(f"Total Susceptibles: {S_total[-1]}")
        print(f"Total Infectados  : {I_total[-1]}")
        print(f"Total Recovered   : {R_total[-1]}")
        if plot_all_iterations:
            plt.plot(S, lw=2, label="Susceptibles")
            plt.plot(I, lw=2, label="Infectados")
            plt.plot(R, lw=1, label="Revovered")

            plt.legend()
            plt.grid()
            plt.show()

    dt=parameters["dt"]
    nt=parameters["n_iter"]
    t_grid = np.linspace(0, nt*dt, nt+1)
    plt.title("Evolution over time of S, I, R")
    plt.plot(t_grid, S_total, "-.", lw=2, label="(S)usceptibles")
    plt.plot(t_grid, I_total, "-", lw=3, label="(I)nfected")
    plt.plot(t_grid, R_total, "--", lw=2, label="(R)ecovered")
    plt.xlabel("Time")
    plt.ylabel("Individuals")
    plt.grid()
    plt.legend()
    fig_file = "test1.png"
    plt.savefig(fig_file)
    plt.show()

    print ( '' )
    print ( f"{test_name}:" )
    print ( '  Normal end of execution.' )
    print ( '' )
    print ( time.ctime ( time.time() ) )

if ( __name__ == '__main__' ):

    test_1()
