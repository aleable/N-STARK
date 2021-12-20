![alt text](https://github.com/aleable/N-STARK/blob/main/misc/logo.svg)

## Implementation details

### Table of Contents  

- [What's included](#whats-included)  
- [Requirements](#requirements)  
- [How to use](#how-to-use)  
    - [Parameters](#parameters)  
- [I/O format](#io-format)  
    - [Input](#input)
    - [Output](#output)
- [Additional implementation details](#additional-implementation-details)
    - [Fourier series computation](#fourier-series-computation)
    - [Convergence criteria](#convergence-criteria)
- [Usage examples](#usage-examples)  
- [Contacts](#contacts)
- [License](#license)


## What's included

- ```dashboard.ipynb```: Jupyter notebook containing an easy-to-use interface with N-STARK
- ```dashboard_misc.ipynb```: Complementary functions needed by ```dashboard.ipynb```
- ```main.py```: main function containing the N-STARK class
- ```initialization.py```: initialization of the routing problem, i.e. construction of the graph topology and on the forcings
- ```fourier.py```: script computing the matrix of Fourier coefficients to run our new dynamics
- ```dynamics.py```: finite difference scheme of the dynamical systems

## Requirements

All the dependencies needed to run the algorithms can be installed using ```setup.py```.
This script can be executed with the command:

```bash
python setup.py
```

Now, you are ready to use the code! To do so, you can simply use the notebook ```dashboard.ipynb```, from which you can access our solver. <br/>

## How to use

### Parameters

The parameters you can pass to N-STARK are described in detail in ```dashboard.ipynb```. They are:

- *Problem construction*

    - ```topol``` = ```"delaunay"``` or ```"real"```: choose which graph topology you want to build, use ```"synth"``` to generate a [Delaunay triangulation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html) of nodes placed randomly in the square [0,1] x [0,1], ```"real"``` to use the Bordeaux bus network
    - ```exec_mode``` = ```"evolution"``` or ```"fourier"```: choose which dynamics you want to run, use ```"evolution"``` to consider time dependent loads, and ```"fourier"``` to use the Fourier coefficients as forcings
    - ```gseed``` (```type=int```): seed for the nodes' coordinates random initialization: x<sub>v</sub> ~ Uniform(0,1), y<sub>v</sub> ~ Uniform(0,1)
    - ```museed``` (```type=int```): seed for the conductivities random initialization μ<sub>e</sub>(0) ~ Uniform(0,1)
    - ```massseed``` (```type=int```): seed for the random initialization of the mass loads (more later)
    - ```pflux``` (```type=float```): 0 < β < 2 (β = 2 - γ), regulatory parameter for traffic congestion

- *Forcing construction*

    The construction of the right hand side is different for synthetic and real networks. In particular, using our code you can test the the dynamics in these two setups:

    1. If ```topol``` = ```"delaunay"```, then the input loads are constructed injecting S<sub>v</sub> (t) = 20 A<sub>v</sub> cos(2 π m<sub>1</sub> t) +  10 A<sub>v</sub> cos(2 π m<sub>2</sub> t) +  5 A<sub>v</sub> cos(2 π m<sub>3</sub> t), for all v ∈ V, with A ~ Dirichlet(**k** = **1**<sub>|V|</sub>) - 1/|V|

    2. If ```topol``` = ```"real"```, the input loads are constructed as S<sub>v</sub>(t) = ∑<sub>v ∈ {m<sub>1</sub>, m<sub>2</sub>, m<sub>3</sub>} </sub> A<sub>v</sub><sup>i</sup> cos (2 π i t), with A<sup>i</sup> = 100/```n_choice``` for ```n_choice``` randomly extracted nodes, and A<sup>i</sup> = -100/(|V|- ```n_choice```) for the remaining |V|-```n_choices``` ones. If m<sub>1</sub>, m<sub>2</sub> or m<sub>3</sub> is equal to -1, the amplitudes are set to 0

    The parameters are:

    - ```mode1``` (```type=int```): m<sub>1</sub>
    - ```mode2``` (```type=int```): m<sub>2</sub>
    - ```mode3``` (```type=int```): m<sub>3</sub>
    - ```n_choice``` (```type=int```): number of nodes with positive amplitudes forcings

- *Dynamics parameters*

    - ```time_step``` (```type=float```): time step for the dynamics with time dependent forcings (it should be << ```time_step_f```)
    - ```time_step_f```  (```type=float```):  time step dynamics with Fourier coefficients
    - ```time_meta``` (```type=float```): upped bound time dynamics with time dependent forcings
    - ```tot_time``` (```type=int```): upper bound on the number of time steps (safety variable)
    - ```relax_linsys``` (```type=float```): relaxation for the weighted Laplacian  
    - ```tau_cond_dyn``` (```type=float```): threshold for convergence, using the conductivities


- *Misc*

    - ```VERBOSE``` (```type=bool```)

## I/O format

### Input

If you want to test the code on the Bordeaux bus network you need these input files:

- ```bus_edges.pkl``` (```type=list of tuples```): list of the edges of the network, each element is of the form ```u v```, with ```u``` and ```v``` nodes connected by an edge ```e```
- ```bus_nodes.pkl``` (```type=list of int```): list of nodes of the network
- ```bus_nodes_coord.pkl``` (```type=dict```): dictionary containing the nodes' coordinates. The formatting is ```{v : (x[v], y[v])}```, with ```v``` node index, and ```x[v]```,```y[v]``` cartesian coordinates of ```v```

### Output

The outputs returned by our schemes are:

- ```self.tdens_stack```: evolution of the conductivities for the two dynamics. This can be exported in the setting of time dependent loads using ```export_evolution()```, and when taking the Fourier coefficients using ```export_fourier()```
- ```self.flux_stack```: evolution of the fluxes in the dynamics with time dependent loads. It can be serialized with ```export_evolution()```
- ```self.opttdens_dyn```: stationary conductivities of the dynamics with Fourier coefficients. It can be serialized with ```export_fourier()```
- ```self.cost_stack```: evolution of the conductivities for the two dynamics. This can be exported in the setting of time dependent loads using ```export_evolution()```, and when taking the Fourier coefficients using ```export_fourier()```

Additionally, using the module ```ot_setup()```, one can serialize:

- ```self.g```: graph topology
- ```self.length```: length of the edges
- ```self.C```: forcing matrix constructed with Fourier coefficients
- ```self.forcing_f```: time dependent forcing function

## Additional implementation details

### Fourier series computation

To calculate the Fourier series expansion of the input loads, we use def ```fs_coeff()``` and ```coefficient_matrix()```, both defined in the ```fourier.py``` model. <br>
The first computes the complex coefficient of a periodic function as:

```python
def fs_coeff(f, N, return_complex):
    """Calculates the first 2*N+1 Fourier series coeff. of a periodic function.
    Taken from: https://stackoverflow.com/questions/4258106/how-to-calculate-a-fourier-series-in-numpy

    f(t) ~= a0/2+ sum_{k=1}^{N} ( a_k*cos(2*pi*k*t/T) + b_k*sin(2*pi*k*t/T) )

    If return_complex is set to True, it returns instead the coefficients
    {c0,c1,c2,...}
    such that:

    f(t) ~= sum_{k=-N}^{N} c_k * exp(i*2*pi*k*t/T)

    where we define c_{-n} = complex_conjugate(c_{n})

    Parameters
    ----------
    f : the periodic function, a callable like f(t)
    N : the function will return the first N_max + 1 Fourier coeff.
    return_complex : bool parameter, set to True in our implementation

    Returns
    -------
    if return_complex == False, the function returns:
        a0 : float
        a,b : numpy float arrays describing respectively the cosine and sine coeff.
    if return_complex == True, the function returns:
        c : numpy 1-dimensional complex-valued array of size N+1
    """

    T = 1                                   # period assumed to be 1
    f_sample = 2 * N                        # sample frequency
    t, dt = np.linspace(0, T, f_sample + 2, endpoint=False, retstep=True)

    y = np.fft.rfft(f(t)) / t.size          # here the DFT is computed, see:
                                            # https://numpy.org/doc/stable/reference/generated/numpy.fft.rfft.html
    if return_complex:
        return y
    else:
        y *= 2
        return y[0].real, y[1:-1].real, -y[1:-1].imag
```

The second function combines the coefficients as follows:
```python
def coefficient_matrix(forcings, N_chosen, return_complex):
    """Compute the Fourier coefficients matrix

    Parameters:
        forcings: function, time vaying forcings
        T: float, period
        N_chosen: int, truncation of Fourier series
        return_complex: bool, complex or real values in coeff matrix

    Returns:
        coeff_matrix: numpy.array(complex), matrix C
    """

    # coefficients calculation
    coeff_ = fs_coeff(forcings, N_chosen, return_complex)

    # coefficients flattening, we omit the code to clarify the important steps
    [...]
    coeff_extended = ... # array with complex fourier coefficients

    # building the Fourier coefficient matrix, C
    nnode = coeff_.shape[0]
    coeff_matrix = np.zeros((nnode, nnode), dtype=float)
    for i in range(nnode):
        for j in range(nnode):
            coeff_matrix[i][j] = np.dot(coeff_extended[i], coeff_extended[j].conj())

    coeff_matrix = csr_matrix(coeff_matrix)

    return coeff_matrix
```

### Convergence criteria

The convergence criteria we choose for the dynamics with Fourier coefficients is the following: we stop the scheme if the cost variation Δℒ<sub>γ</sub> := (ℒ<sub>γ</sub>(t+1) - ℒ<sub>γ</sub>(t)) / (ℒ<sub>γ</sub>(t+1)```time_step_f```) is below a certain threshold. In detail:
```python
# abs_diff_cost is the difference of the cost in two consecutive time steps
if self.pflux >= 0.0:                       # sanity check
    if abs_diff_cost < self.tau_cost_dyn:   # threshold check
        convergence_achieved = True
        print("dynamics fourier converged [stopping criteria]")
```

## Usage examples

For a basic usage example of the code you can simply take a look look at ```dashboard.ipynb```. <br/>
The execution of N-STARK is performed in three steps:
- *Initialization*: first, you need to pass the necessary parameters to the N-STARK class. Similarly to what done in ```dashboard.ipynb```, you can simply run:

```python
from main import *  # import N-STARK class
[...]   # you may want to load other modules

# problem construction
topol = "real"
exec_mode = "fourier"
gseed = 0
museed = 0
massseed = 0
pflux = 1.0

# forcing construction
mode1 = 1
mode2 = -1
mode3 = -1
n_choice = 1

# time dependent loads dynamucs
relax_linsys = 1e-5
time_step = 0.01         # << time_step_f to capture the behavior of the forcings
time_meta = 10
tot_time = 1000

# dynamics with Fourier coefficients
time_step_f = 0.1        # time step Fourier dynamics
tau_cost_dyn = 1e-5      # threshold convergence cost dynamics

VERBOSE = False

# create an object
nstark_f = NSTARK(exec_mode, topol, gseed, museed, pflux, VERBOSE, relax_linsys,
               tau_cost_dyn, time_step, time_step_f, time_meta, tot_time,
               massseed, mode1, mode2, mode3, n_choice)

# construct the topology
nstark_f.ot_setup()     # step 1
```
- *Execution*: now, you can execute our scheme. You just need to run:
```python
nstark_f.dyn_exec()     # step 2
```

- *Serialization*: once the algorithm reached convergence, you can export the results with:

```python
tdens_stack_f, tdens_opt_f, cost_stack_f = nstark_f.export_fourier()    # step 3
```

Similarly, if you want to run the dynamics with time dependent loads, you can simply execute:
```python
[...]       # same initialization of the parameters
# IMPORTANT: we need to change
exec_mode = "evolution"

nstark = NSTARK(exec_mode, topol, gseed, museed, pflux, VERBOSE, relax_linsys,
               tau_cost_dyn, time_step, time_step_f, time_meta, tot_time,
               massseed, mode1, mode2, mode3, n_choice)

# we also export the inputs so we can plot the results, see how in dashboard.ipynb
graph, length, C, forcing = nstark.ot_setup()                   # step 1
nstark.dyn_exec()                                               # step 2
tdens_stack, flux_stack, cost_stack = nstark.export_evolution() # step 3
```

## Contacts

For any issues or questions, feel free to contact us sending an email to <a href="alessandro.lonardi@tuebingen.mpg.de">alessandro.lonardi@tuebingen.mpg.de</a>.

## License

Copyright (c) 2021 <a href="https://aleable.github.io/">Alessandro Lonardi</a> and <a href="https://www.cdebacco.com/">Caterina De Bacco</a>

<sub>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:</sub>

<sub>The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.</sub>

<sub>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.</sub>
