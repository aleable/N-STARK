"""
N-STARK (Non-STAtionary loads Routing on networKs) -- https://github.com/aleable/N-STARK
Contributors:
    Alessandro Lonardi
    Caterina De Bacco
"""

import numpy as np
import networkx as nx
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse import identity


def tdensinit(self):
    """Initialization of the conductivities: mu_e ~ U(0,1)

    Parameters:
        seed: int, seed random initialization conductivities

    Returns:
        self.tdens: np.array, initialized conductivities
    """

    prng = np.random.RandomState(seed=self.museed)
    self.tdens = np.array([prng.uniform(0, 1) for i in range(self.g.number_of_edges())])

    return self.tdens


def dyn(self):
    """Execute dynamics

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

    if self.exec_mode in "evolution":
        print("* running [long run dynamics]")
    if self.exec_mode in "fourier":
        print("* running [fourier dynamics]")

    time_iteration = 0
    prng = np.random.RandomState(seed=self.museed)

    # initialization quantities of dynamics
    nnode = self.g.number_of_nodes()
    inc_mat = csr_matrix(nx.incidence_matrix(self.g, nodelist=list(range(nnode)), oriented=True))
    inc_transpose = csr_matrix(inc_mat.transpose())
    inv_len_mat = diags(1/self.length, 0)
    td_mat = diags(self.tdens, 0)

    stiff = inc_mat * td_mat * inv_len_mat * inc_transpose      # weighted Laplacian
    stiff_relax = stiff + self.relax_linsys*identity(nnode)     # avoid empty kernel

    if self.exec_mode in "evolution":
        forcing = np.array(self.forcing_f(0))                   # forcing at time step t = 0
        pot = spsolve(stiff_relax, forcing.transpose())         # potential

    if self.exec_mode in "fourier":
        a = td_mat * inv_len_mat * inc_transpose * np.linalg.pinv(csr_matrix.todense(stiff_relax))
        flux_mat_squared = np.diagonal(a * self.C * a.transpose())       # squared fluxes

    # executing dynamics
    convergence_achieved = False
    cost = 0
    while not convergence_achieved:

        ### UPDATE TDENS-POT SYSTEM
        tdens_old = self.tdens

        ### TIME VARYING FORCINGS
        if self.exec_mode in "evolution":
            pot_old = pot

            ### DYNAMICS - ONE STEP
            self.tdens, pot, info = update(self, pot, inc_mat, inc_transpose, inv_len_mat, time_iteration)
            self.tdens_stack.append(self.tdens)

            # cost
            td_mat = diags(self.tdens, 0)  # diagonal matrix conductivities
            flux = td_mat * inv_len_mat * inc_transpose * pot  # right hand side ode
            cost = np.sum((self.length/self.tdens) * flux**2) + 1/(2-self.pflux)*np.sum(self.length*self.tdens**(2-self.pflux))
            self.cost_stack.append(cost)
            self.flux_stack.append(flux)

        ### FOURIER COEFFICIENTS
        if self.exec_mode in "fourier":
            flux_mat_squared_old = flux_mat_squared

            ### DYNAMICS - ONE STEP
            self.tdens, flux_mat_squared, info = update_fourier(self, flux_mat_squared, inc_mat, inc_transpose, inv_len_mat)
            self.tdens_stack.append(self.tdens)

        if info != 0:
            self.tdens = tdens_old + prng.rand(*self.tdens.shape) * np.mean(tdens_old) / 1000.
            if self.exec_mode in "evolution":
                pot = pot_old + prng.rand(*pot.shape) * np.mean(pot_old) / 1000.
            if self.exec_mode in "fourier":
                flux_mat_squared = flux_mat_squared_old + prng.rand(*flux_mat_squared.shape) * np.mean(flux_mat_squared_old) / 1000.

        ### OUT OF RUNNING TIME
        if time_iteration*self.time_step >= self.tot_time:
            convergence_achieved = True
            self.tdens = tdens_old
            print("\tERROR: dynamics did NOT converge [iteration > maxit]")

        ### OSCILLATING STATE CONVERGENCE
        if time_iteration*self.time_step >= self.time_meta and self.exec_mode in "evolution":
            convergence_achieved = True
            self.tdens = tdens_old
            print("\tdynamics time dependent loads converged [oscillating state]")

        ### CONVERGENCE FOURIER
        if self.exec_mode in "fourier":
            convergence_achieved, cost, abs_diff_cost = cost_convergence_fourier(self, flux_mat_squared, cost, convergence_achieved)
            self.cost_stack.append(cost)

        ### VERBOSE
        if self.exec_mode in "evolution" and self.verbose:
            if time_iteration % 5 == 0:
                print('\r', 'beta=%1.2f, it=%3d, time=%1.2f, J=%8.2e' % (self.pflux, time_iteration,
                                                                             time_iteration*self.time_step, cost))
        if self.exec_mode in "fourier" and self.verbose:
            if time_iteration % 5 == 0:
                print('\r', 'beta=%1.2f, it=%3d, time=%1.2f, J_err=%8.2e' % (self.pflux, time_iteration,
                                                                             time_iteration*self.time_step, abs_diff_cost))

        time_iteration += 1

    if self.exec_mode == "evolution":
        return self.tdens_stack, self.flux_stack, self.cost_stack
    if self.exec_mode == "fourier":
        return self.tdens_stack, self.tdens, self.cost_stack


def update(self, pot, inc_mat, inc_transpose, inv_len_mat, time_iteration):
    """One step update dynamics for time-varying forcings

    Parameters:
        pot: np.array, potential matrix on nodes
        inc_mat: sparse.matrix, oriented incidence matrix
        inc_transpose: sparse.matrix, oriented incidence matrix transposed
        inv_len_mat: sparse.matrix, diagonal matrix 1/l_e

    Returns:
        self.tdens: np.array, updated conductivities
        pot: np.array, updated potential matrix on nodes
        info: bool, sanity check flag spsolve
    """

    grad = inv_len_mat * inc_transpose * pot
    rhs_ode = (self.tdens**self.pflux)*(grad**2) - self.tdens

    self.tdens = self.tdens + self.time_step*rhs_ode
    td_mat = diags(self.tdens, 0)

    stiff = inc_mat*td_mat*inv_len_mat*inc_transpose                              # update stiffness matrix
    stiff_relax = stiff + self.relax_linsys * identity(self.g.number_of_nodes())  # avoid zero kernel

    forcing = np.array(self.forcing_f(time_iteration*self.time_step))             # time evaluation forcing
    pot = spsolve(stiff_relax, forcing.transpose())                               # update potentials

    if np.any(np.isnan(pot)):
        info = -1
        pass
    else:
        info = 0

    return self.tdens, pot, info


def update_fourier(self, flux_mat_squared, inc_mat, inc_transpose, inv_len_mat):
    """One step update dynamics with Fourier coefficients

    Parameters:
            flux_mat_squared: np.array, old flux squared
        inc_mat: sparse.matrix, oriented incidence matrix
        inc_transpose: sparse.matrix, oriented incidence matrix transposed
        inv_len_mat: sparse.matrix, diagonal matrix 1/l_e

    Returns:
        self.tdens: np.array, updated conductivities
        flux_mat_squared: np.array, updated flux squared
        info: bool, sanity check flag spsolve
    """

    rhs_ode = (self.tdens ** (self.pflux - 2)) * flux_mat_squared - self.tdens      # right hand side ode
    self.tdens = self.tdens + self.time_step_f * rhs_ode                            # update conductivities
    td_mat = diags(self.tdens, 0)

    stiff = inc_mat * td_mat * inv_len_mat * inc_transpose          # update stiffness matrix
    stiff_relax = stiff + self.relax_linsys * identity(self.g.number_of_nodes())  # avoid zero kernel

    a = td_mat * inv_len_mat * inc_transpose * np.linalg.pinv(csr_matrix.todense(stiff_relax))
    flux_mat_squared = np.diagonal(a * self.C * a.transpose())           # update fluxes

    if np.any(np.isnan(flux_mat_squared)):
        info = -1
        pass
    else:
        info = 0

    return self.tdens, flux_mat_squared, info


def cost_convergence_fourier(self, flux_mat_squared, cost, convergence_achieved):
    """Evaluating convergence

    Parameters:
        flux_mat_squared: np.array, flux squared
        cost: float, cost
        convergence_achieved: bool, convergence flag

    Returns:
        convergence_achieved: bool, updated convergence flag
        cost_update: float, updated cost
        abs_diff_cost: float, difference cost
    """

    cost_update = np.sum((self.length/self.tdens)*flux_mat_squared) + 1/(2-self.pflux)*np.sum(self.length*self.tdens**(2-self.pflux))
    abs_diff_cost = abs(cost_update - cost)/(self.time_step_f*cost_update)

    if self.pflux >= 0.0:    # sanity check
        if abs_diff_cost < self.tau_cost_dyn:
            convergence_achieved = True
            print("\tdynamics fourier converged [stopping criteria]")

    return convergence_achieved, cost_update, abs_diff_cost