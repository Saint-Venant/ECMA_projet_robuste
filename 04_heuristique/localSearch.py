import numpy as np
from transformation_fichiers import load_data
from classes import make_yki, make_x, make_edges
from SP2 import SP2k, SP1
import time



class Data:
    def __init__(self, n, L, W, K, B, w_v, W_v, lh, coordonnees, lij):
        self.n = n
        self.L = L
        self.W = W
        self.K = K
        self.B = B
        self.w_v = w_v
        self.W_v = W_v
        self.lh = lh
        self.coordonnees = coordonnees
        self.lij = lij
        self.edges = make_edges(n)

def init_yki(n, K):
    yki = [n*[0] for k in range(K)]
    for i in range(n):
        k = np.random.randint(K)
        yki[k][i] = 1
    return yki

def compute_SP2k(yki, data):
    values_SP2k = []
    valMax = 0
    indexMax = -1
    for k in range(data.K):
        val = SP2k(yki, k, data.w_v, data.W_v, data.W, data.n)
        values_SP2k.append(val)
        if val > valMax:
            valMax = val
            indexMax = k
    return values_SP2k, valMax, indexMax

def getFeasibleSolution(data, verbose=False):
    n = data.n
    K = data.K
    B = data.B
    yki = init_yki(n, K)

    # compute the values of all SP2-k sub-problems and get the max
    values_SP2k, valMax, kMax = compute_SP2k(yki, data)
    if verbose:
        print(values_SP2k)
    feasible = (valMax <= B)

    maxIter = 10000
    iteration = 0
    while not(feasible) and (iteration < maxIter):
        # change one vertex to a random cluster
        vertices = [i for i in range(n) if yki[kMax][i] == 1]
        assert(len(vertices) > 0)
        np.random.shuffle(vertices)

        acceptSwitch = False
        for i_change in vertices:
            for new_k in range(K):
                yki[kMax][i_change] = 0
                yki[new_k][i_change] = 1

                _, newValMax, _ = compute_SP2k(yki, data)

                if (newValMax < valMax):
                    acceptSwitch = True
                    break
                else:
                    yki[kMax][i_change] = 1
                    yki[new_k][i_change] = 0
            if acceptSwitch:
                break

        if not(acceptSwitch):
            # random change
            i_change = vertices[0]
            new_k = np.random.randint(K)
            yki[kMax][i_change] = 0
            yki[new_k][i_change] = 1

        values_SP2k, valMax, kMax = compute_SP2k(yki, data)
        if verbose:
            print(values_SP2k)
        feasible = (valMax <= B)
        iteration += 1

    if feasible:
        print("Found feasible solution (in {} iterations)\n\n".format(iteration))
    else:
        print("No feasible solution found (in {} iterations)\n\n".format(iteration))
    return yki

def V1(yki, n, K):
    yki_bis = [x.copy() for x in yki]
    i = np.random.randint(n)
    for k in range(K):
        yki_bis[k][i] = 0
    k = np.random.randint(K)
    yki_bis[k][i] = 1
    return yki_bis

def V2(yki, n, K):
    yki_bis = [x.copy() for x in yki]
    i1 = np.random.randint(n)
    i2 = np.random.randint(n)
    while (i2 == i1):
        i2 = np.random.randint(n)
    cluster1 = 0
    cluster2 = 0
    for k in range(K):
        if yki_bis[k][i1] == 1:
            cluster1 = k
            yki_bis[k][i1] = 0
        if yki_bis[k][i2] == 2:
            cluster2 = k
            yki_bis[k][i2] = 0
    yki_bis[cluster1][i2] = 1
    yki_bis[cluster2][i1] = 1
    return yki_bis

def isFeasible(yki, data):
    values_SP2k, valMax, kMax = compute_SP2k(yki, data)
    return valMax <= data.B

def draw_V():
    alpha = np.random.random()
    if alpha < 0.5:
        return V1
    else:
        return V2

def descent(yki, data, dt=30, verbose=False):
    n = data.n
    K = data.K
    edges = data.edges
    t_begin = time.time()
    x=make_x(yki, edges, K)
    value = SP1(x, edges, data.lij, data.L, n, data.lh)
    if verbose:
        print(' > value : ', value)
    while (time.time() - t_begin) < dt:
        V = draw_V()
        yki_bis = V(yki, n, K)
        x_bis = make_x(yki_bis, edges, K)
        if isFeasible(yki_bis, data):
            valueBis = SP1(x_bis, edges, data.lij, data.L, n, data.lh)
            if valueBis <= value:
                yki = yki_bis
                value = valueBis
                if verbose:
                    print(' > value : ', value)
    return yki, value

def heuristique(instanceFile, dtSearch=30, nbSearch=4):
    # load data
    instanceFile = '..\\modified_data\\52_berlin_9.dat'
    n,L,W,K,B,w_v,W_v,lh,coordonnees,lij=load_data(path=instanceFile)
    data = Data(n, L, W, K, B, w_v, W_v, lh, coordonnees, lij)

    values = []
    for p in range(nbSearch):
        yki = getFeasibleSolution(data)
        yki, value = descent(yki, data, dt=dtSearch)
        values.append(value)
    print(values)
    return np.min(values)

if __name__ == '__main__':
    #np.random.seed()
    # load data
    instanceFile = '..\\modified_data\\52_berlin_9.dat'
    '''n,L,W,K,B,w_v,W_v,lh,coordonnees,lij=load_data(path=instanceFile)
    data = Data(n, L, W, K, B, w_v, W_v, lh, coordonnees, lij)

    yki = getFeasibleSolution(data)
    yki, _ = descent(yki, data, verbose=True)'''
    heuristique(instanceFile, dtSearch=45, nbSearch=3)
