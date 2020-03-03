import numpy as np

from transformation_fichiers import load_data
from classes import make_yki, make_x, make_edges
from SP2 import SP2k, SP1


def init_yki(n, K):
    yki = [n*[1]]
    for k in range(1, K):
        yki.append(n*[0])
    return yki

def compute_SP2k(yki, w_v, W_v, W, n):
    values_SP2k = []
    valMax = 0
    indexMax = -1
    for k in range(K):
        val = SP2k(yki,k,w_v,W_v,W,n)
        values_SP2k.append(val)
        if val > valMax:
            valMax = val
            indexMax = k
    return values_SP2k, valMax, indexMax


if __name__ == '__main__':
    # load data
    n,L,W,K,B,w_v,W_v,lh,coordonnees,lij=load_data(path='..\\modified_data\\52_berlin_3.dat')
    edges=make_edges(n)
    
    yki = init_yki(n, K)
    x=make_x(yki,edges,K)

    # compute the values of all SP2-k sub-problems and get the max
    values_SP2k, valMax, kMax = compute_SP2k(yki, w_v, W_v, W, n)
    print(values_SP2k)
    feasible = (valMax <= B)

    maxIter = 1000
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

                _, newValMax, _ = compute_SP2k(yki, w_v, W_v, W, n)

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

        values_SP2k, valMax, kMax = compute_SP2k(yki, w_v, W_v, W, n)
        print(values_SP2k)
        feasible = (valMax <= B)
        iteration += 1

    if feasible:
        print("Found feasible solution")
    else:
        print("No feasible solution found")
