"""
N-STARK (Non-STAtionary loads Routing on networKs) -- https://github.com/aleable/N-STARK
Contributors:
    Alessandro Lonardi
    Caterina De Bacco
"""

import networkx as nx
import numpy as np

from initialization import *
from dynamics import *


class NSTARK:
    """Optimal Transport algorithms"""

    def __init__(self, exec_mode, topol, gseed, museed, pflux, verbose, relax_linsys,
                 tau_cost_dyn, time_step, time_step_f, time_meta, tot_time,
                 massseed, mode1, mode2, mode3, n_choice):

        # graph topology
        self.gseed = gseed                                  # Delaunay triangulation seed
        self.g = nx.Graph()                                 # graph topology
        self.length = np.zeros(self.g.number_of_edges())    # length of the edges

        # dynamical system parameters
        self.exec_mode = exec_mode                          # integrated fast dynamics / long run
        self.topol = topol                                # Bordeaux bus network / Delaunay triangulation

        self.pflux = pflux                                  # beta (gamma = 2 - beta)
        self.relax_linsys = relax_linsys                    # relaxation for Laplacian
        self.museed = museed                                # seed init conductivities
        self.masseed = massseed                             # random mass assignation seed
        self.tdens = np.zeros(self.g.number_of_edges())     # conductivities

        self.time_step = time_step                          # time step slow dynamics
        self.time_meta = time_meta                          # stopping time for dynamics with time dependent forcings
        self.time_step_f = time_step_f                      # time step for slow dynamics with Fourier coefficients
        self.tot_time = tot_time                            # stopping time (safety variable)

        self.mode1 = mode1
        self.mode2 = mode2
        self.mode3 = mode3
        self.n_choice = n_choice                            # number of resonant points

        self.C = np.zeros((self.g.number_of_nodes(), self.g.number_of_nodes()))     # matrix of Fourier coefficients
        self.forcing_f = None                               # time-varying forcing function

        # convergence parameters
        self.tau_cost_dyn = tau_cost_dyn                    # threshold convergence conductivities slow dynamics

        # misc
        self.verbose = verbose                              # verbose

        # variables to serialize
        self.opttdens_dyn = np.zeros(self.g.number_of_edges())
        self.tdens_stack = []
        self.flux_stack = []
        self.cost_stack = []

    def ot_setup(self):
        """Building graph topology

        Returns:
            self.g: nx.Graph(), network topology
            self.length: np.array, lengths of edges
            self.C: sparse.array, Fourier coefficient matrix
            self.forcing_f: function, forcing functions
        """

        print("* graph topology construction")

        if self.topol == "delaunay":
            self.g, self.length = delaunay_topology(self)
            self.C, self.forcing_f = rhs(self)

        if self.topol == "real":
            self.g, self.length = real_topology(self)
            self.C, self.forcing_f = rhs_real(self)

        return self.g, self.length, self.C, self.forcing_f

    def dyn_exec(self):
        """Run dynamics

        Returns:
            if self.exec_mode == "evolution":
                self.tdens_stack: np.array, conductivities over time
                self.flux_stack: np.array, fluxes over time
                self.cost_stack: np.array, cost over time
            if self.exec_mode == "fourier":
                self.tdens_stack: np.array, conductivities over time
                self.tdens: np.array, conductivities at convergence
                self.cost_stack: np.array, cost over time

        """

        # conductivities initialization
        self.tdens = tdensinit(self)
        if self.exec_mode == "evolution":       # dynamics with time-varying forcings
            self.tdens_stack, self.flux_stack, self.cost_stack = dyn(self)
        if self.exec_mode == "fourier":         # dynamics with Fourier coefficients
            self.tdens_stack, self.opttdens_dyn, self.cost_stack = dyn(self)

    def export_evolution(self):
        """Export variables dynamics with time-varying forcings

        Returns:
            self.tdens_stack: np.array, conductivities over time
            np.flux_stack: np.array, fluxes over time
            self.cost_stack: np.array, cost over time
            """

        print("* serialize time dependent loads dynamics")

        return np.array(self.tdens_stack), np.array(self.flux_stack), np.array(self.cost_stack)

    def export_fourier(self):
        """Export variables dynamics with Fourier expansion

         Returns:
            self.tdens_stack: np.array, conductivities over time
            self.tdens: np.array, conductivities at convergence
            self.cost_stack: np.array, cost over time
            """

        print("* serialize Fourier dynamics")

        return np.array(self.tdens_stack), np.array(self.opttdens_dyn), np.array(self.cost_stack)
