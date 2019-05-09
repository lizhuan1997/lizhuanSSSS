import numpy as np
from buckyballgenerator import generator
import time
import matplotlib.pyplot as plt

def energy_per_config(structure, vector):
    energy = 0
    for edge in edges:
        i = edge[0]
        j = edge[1]
        energy += J * vector[i] * vector[j]
    return energy

def gapp(structure, statevector, J, flipnode):
        oldEnergy = 0
        newEnergy = 0
        flip = -statevector[flipnode]
        neigh = structure[flipnode]
        for i in range(3):
            oldEnergy += J * statevector[flipnode] * statevector[neigh[i]]
            newEnergy += J * flip * statevector[neigh[i]]
        gap = oldEnergy - newEnergy
        return gap, flip

structure, edges = generator()
E_min = -66
E_max = 90
J = 1
beta = 1
logG = {}
Hist = {}
f = 2.718281828

for E in range(E_min, E_max+1):
    logG[E] = 1
    Hist[E] = 0
statevector = np.random.rand(60)
# statevector = np.where(statevector > 0.5, 1, -1)

for i in range(60):
    if statevector[i] > 0.5:
        statevector[i] = 1
    else:
        statevector[i] = -1
start = time.time()

statevector = np.ones(60)
E = energy_per_config(structure, statevector)

while(f > 1+10e-8):

    for u in range(E_min, E_max+1):
        Hist[u] = 0
    for iter in range(2000):
        for j in range(60):
        # j = np.random.randint(0, 60)
            gap, flip = gapp(structure, statevector, J, j)
            E_new = E - gap
            if np.exp(logG[E] - logG[E_new]) > 1:
                E = E_new
                logG[E] += np.log(f)
                Hist[E] += 1
                statevector[j] = flip
            else:
                p = np.random.rand()
                if p < np.exp(logG[E] - logG[E_new]):
                    E = E_new
                    logG[E] += np.log(f)
                    Hist[E] += 1
                    statevector[j] = flip
                else:
                    logG[E] += np.log(f)
                    Hist[E] += 1
   
        Hist_list = [Hist[i] for i in np.arange(E_min, E_max+1, 2)]
        Hist_list = Hist_list[0:len(Hist_list)-3]
        error = (np.max(Hist_list) - np.mean(Hist_list)) / np.mean(Hist_list)
        # if iter % 1000 == 0:
        #     print(error)

            # plt.ion()
            # X = range(E_min, E_max)
            # Y = np.array([logG[i] for i in X])

            # index = np.nonzero(Y==1)
            # # print(index)

            # Y = np.delete(Y, index)
            # X = np.delete(X, index)
            # plt.plot(X, Y)
            # plt.show()
            # plt.pause(0.01)
            # plt.clf()
    f = f ** 0.5
           

X = list(range(E_min, E_max + 1))
Y = np.zeros(E_max - E_min + 1)
for i in range(Y.shape[0]):
    Y[i] = logG[X[i]]

Y_remove = [Y[i] for i in np.arange(0, E_max - E_min, 2)]
Y = [Y[i] for i in np.arange(0, E_max - E_min, 2)]

while 1 in Y_remove:
    Y_remove.remove(1)

Y = np.exp(Y - np.min(Y_remove))
print(np.min(Y_remove))
print(Y[0] / np.sum(Y) * 2**60)
E = (np.arange(E_min, E_max, 2))
Z = np.dot(Y, np.exp(-E))/sum(Y)*2**60
stop = time.time()

print(np.log(Z)/60, stop-start)
