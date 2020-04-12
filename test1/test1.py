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
import schemes.SI_react_diffusion_1d as SI_scheme
from schemes.SI_react_diffusion_1d import test_name, run_test
#plt.style.use('seaborn-pastel')
plt.style.use('seaborn-dark')
plt.rcParams["figure.dpi"] = 300

def test_1 ( ):
    from copy import copy

    parameters = copy(SI_scheme.default_parameters)
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

    plot_each_iteration = False
    S_total = [] # Total number of susceptible at each time step
    I_total = [] # Total number of infected at each time step
    for S, I in run_test ( parameters, plot_options={ "save": False } ):
        s_array = S.get_local(); S_total.append(np.sum(s_array))
        i_array = I.get_local(); I_total.append(np.sum(i_array))
        print(f"Total Susceptibles: {S_total[-1]}")
        print(f"Total Infectados  : {I_total[-1]}")
        if plot_each_iteration:
            plt.plot(S, lw=2, label="Susceptibles")
            plt.plot(I, lw=2, label="Infectados")
            SI = S+I
            plt.plot(SI, lw=1, label="Suma I+S")

            plt.legend()
            plt.grid()
            plt.show()

    dt=parameters["dt"]
    nt=parameters["n_iter"]
    t_grid = np.linspace(0, nt*dt, nt+1)
    plt.title("Evolution over time of S, I")
    plt.plot(t_grid, S_total, lw=2, label="(S)usceptibles")
    plt.plot(t_grid, I_total, lw=2, label="(I)nfected")
    plt.xlabel("Time")
    plt.ylabel("Individuals")
    plt.grid()
    plt.legend()
    plt.show()

    print ( '' )
    print ( f"{test_name}:" )
    print ( '  Normal end of execution.' )
    print ( '' )
    print ( time.ctime ( time.time() ) )

if ( __name__ == '__main__' ):

    test_1()
