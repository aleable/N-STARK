"""
N-STARK (Non-STAtionary loads Routing on networKs) -- https://github.com/aleable/N-STARK
Contributors:
    Alessandro Lonardi
    Caterina De Bacco
"""

from fourier import coefficient_matrix
import numpy as np
from scipy.spatial import distance, Delaunay
import networkx as nx
import pickle5 as pkl


def delaunay_topology(self):
    """Generation of graph from Delaunay triangulation

    Returns:
        self.g: nx.Graph(), Graph topology
        self.length: np.array, lengths of edges
    """

    # delaunay triangulation
    nnode = 20
    gs = np.random.RandomState(seed=self.gseed)
    x = np.array([gs.uniform(0, 1) for i in range(nnode)])
    y = np.array([gs.uniform(0, 1) for i in range(nnode)])
    points = np.transpose([x, y])
    tri = Delaunay(points)
    paths = np.append(tri.simplices, np.transpose([tri.simplices[:, 0]]), axis=1)
    for path in paths:
        nx.add_path(self.g, path)

    # nodes coordinates
    for inode, xy in enumerate(points):
        self.g.nodes[inode]["pos"] = tuple(xy)

    # length edges
    self.length = np.zeros(self.g.number_of_edges())
    for i, edge in enumerate(self.g.edges()):
        self.length[i] = distance.euclidean(x[edge[0]], y[edge[1]])


    return self.g, self.length


def real_topology(self):
    """Generation of Bordeaux bus graph

    Returns:
        self.g: nx.Graph(), Graph topology
        self.length: np.array, lengths of edges
    """

    nodes_original = np.array(pkl.load(open("../data/input/bus_nodes.pkl", "rb")))
    edges_original = pkl.load(open("../data/input/bus_edges.pkl", "rb"))
    coord_original = pkl.load(open("../data/input/bus_nodes_coord.pkl", "rb"))

    coord_x = np.array(list(coord_original.values()))[:,0]
    coord_y = np.array(list(coord_original.values()))[:,1]
    min_x = min(coord_x)
    max_x = max(coord_x)
    min_y = min(coord_y)
    max_y = max(coord_y)
    # trimming the graph to a central area of the city
    nodes_original = np.array([node for node in nodes_original if coord_original[node][0] > min_x + abs(max_x - min_x)/3 and coord_original[node][0] < min_x + 2*abs(max_x - min_x)/3 and coord_original[node][1] > min_y + abs(max_y - min_y)/3 and coord_original[node][1] < min_y + 2*abs(max_y - min_y)/3])

    nodes = [int(np.where(nodes_original == node)[0]) for node in nodes_original]
    edges = []

    for edge in edges_original:
        try:
            index_1 = int(np.where(nodes_original == edge[0])[0])
            index_2 = int(np.where(nodes_original == edge[1])[0])
            edges.append((index_1, index_2))
        except:
            pass

    coord = {}
    for i in coord_original.items():
        try:
            coord[int(np.where(np.array(nodes_original) == i[0])[0])] = i[1]
        except:
            pass

    self.g.add_nodes_from(nodes)
    self.g.add_edges_from(edges)

    # nodes coordinates
    for it in coord.items():
        self.g.nodes[it[0]]["pos"] = tuple(it[1])

    # length edges
    self.length = np.zeros(self.g.number_of_edges())
    for i, edge in enumerate(self.g.edges()):
        self.length[i] = distance.euclidean(coord[edge[0]], coord[edge[1]])

    return self.g, self.length


def rhs(self):
    """Building the Fourier coefficients matrix C

    Returns:
        self.C: Fourier coefficients matrix for the integrated dynamics
    """

    N = self.g.number_of_nodes()

    def forcings(t):
        """Define the time-varying input loads acting on the nodes of the network

        Parameters:
            t: float, time

        Returns:
            forcing: function, loads
        """

        omega = 2 * np.pi  # period T is assumed to be 1

        np.random.seed(self.masseed)

        index_1 = list(set(np.random.randint(low=0, high=N, size=20000)))[:int(N)]  # index A1
        index_2 = list(set(np.random.randint(low=0, high=N, size=20000)))[:int(N)]  # index A2
        index_3 = list(set(np.random.randint(low=0, high=N, size=20000)))[:int(N)]  # index A3

        A1 = np.zeros(N)
        A2 = np.zeros(N)
        A3 = np.zeros(N)
        A1[index_1] = (np.random.dirichlet(np.ones(N), size=1)[0]-1/N)*20
        A2[index_2] = (np.random.dirichlet(np.ones(N), size=1)[0]-1/N)*10
        A3[index_3] = (np.random.dirichlet(np.ones(N), size=1)[0]-1/N)*5

        mode1 = self.mode1               # oscillation mode 1
        mode2 = self.mode2               # oscillation mode 2
        mode3 = self.mode3               # oscillation mode 2

        return tuple(A1[i] * np.cos(omega * mode1 * t) +
                     A2[i] * np.cos(omega * mode2 * t) +
                     A3[i] * np.cos(omega * mode3 * t) for i in range(N))

    N_chosen = 10
    return_complex = True
    self.C = coefficient_matrix(forcings, N_chosen, return_complex)
    self.forcings = forcings

    return self.C, self.forcings


def rhs_real(self):
    """Building the Fourier coefficients matrix C

    Returns:
        self.C: Fourier coefficients matrix for the integrated dynamics
    """

    N = self.g.number_of_nodes()

    def forcings(t):
        """Define the time-varying input loads acting on the nodes of the network

        Parameters:
            t: float, time

        Returns:
            forcing: function, loads
        """

        np.random.seed(self.masseed)
        index_1 = list(set(np.random.randint(low=0, high=N, size=20000)))[:int(N)]
        omega = 2 * np.pi  # period T is assumed to be 1

        np.random.seed(self.masseed)
        node1 = np.random.choice(np.array(index_1), size=self.n_choice)
        np.random.seed(self.masseed+1000)
        node2 = np.random.choice(np.array(index_1), size=self.n_choice)
        np.random.seed(self.masseed+2000)
        node3 = np.random.choice(np.array(index_1), size=self.n_choice)

        A1 = -np.array([100 / (N - self.n_choice)]*N)
        A1[node1] = np.array(100 / self.n_choice)
        A2 = -np.array([100 / (N - self.n_choice)] * N)
        A2[node2] = np.array(100 / self.n_choice)
        A3 = -np.array([100 / (N - self.n_choice)] * N)
        A3[node3] = np.array(100 / self.n_choice)

        mode1 = self.mode1              # oscillation mode 1
        mode2 = self.mode2              # oscillation mode 2
        mode3 = self.mode3              # oscillation mode 3


        if mode1 == -1:
            A1 = np.zeros(N)
        if mode2 == -1:
            A2 = np.zeros(N)
        if mode3 == -1:
            A3 = np.zeros(N)

        return tuple(A1[i] * np.cos(omega * mode1 * t) +
                     A2[i] * np.cos(omega * mode2 * t) +
                     A3[i] * np.cos(omega * mode3 * t) for i in range(N))

    N_chosen = 10
    return_complex = True
    self.C = coefficient_matrix(forcings, N_chosen, return_complex)
    self.forcings = forcings

    return self.C, self.forcings
