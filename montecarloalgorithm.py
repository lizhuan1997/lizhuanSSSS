import numpy as np
import networkx as nx
import time
import matplotlib.pyplot as plt
from buckyballgenerator import generator


def energy_per_config(structure, vector):
    energy = 0
    
    for edge in edges:
        i = edge[0]
        j = edge[1]
        energy += vector[i] * vector[j]
    return energy


structure, edges = generator()

RANGE = 5000000
beta = 1.5
start = time.time()
statematrix = []
# statevector = np.ones(60)
statevector = np.random.rand(60)
for i in range(60):
    if statevector[i] > 0.5:
        statevector[i] = 1
    else:
        statevector[i] = -1
energy = 0
for k in range(RANGE):
    for j in range(60):
        oldEnergy = 0
        newEnergy = 0
        flip = -statevector[j]
        neigh = structure[j]
        for i in range(3):
            oldEnergy += statevector[j] * statevector[neigh[i]]
            newEnergy += flip * statevector[neigh[i]]
            gap = oldEnergy - newEnergy
        if gap > 0:
            statevector[j] = flip
        else:
            p = np.random.rand()
            if p < np.exp(beta * gap):
                statevector[j] = flip
    energy += energy_per_config(structure, statevector)

    statematrix.append(list(statevector))


stop = time.time()
energy = energy / RANGE
print(energy, stop - start)



# Ground State

state = []

for i in range(RANGE):
    statevec = statematrix[i]
    energy = energy_per_config(structure, statevec)
    if energy == -66:
        state.append(statevec)
state = np.unique(state, axis=0)
print(state.shape)
exit()
# entropy

statematrixunique = list(np.unique(statematrix, axis=0))

entropy = 0
for i in range(len(statematrixunique)):
    vector = list(statematrixunique[i])
    freq = statematrix.count(vector) / RANGE
    entropy += freq * np.log(freq)
print("entropy :" , entropy)
