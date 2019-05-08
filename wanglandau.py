import numpy as np
from buckyballgenerator import generator
import time

def energy_per_config(structure, vector):
    energy = 0
    for edge in edges:
        i = edge[0]
        j = edge[1]
        energy += J * vector[i] * vector[j]
    return energy


structure, edges = generator()

iteration = 100000
E_min = -66
E_max = 90
J = 1
beta = 1
logG = {}
Hist = {}
f = 1
for E in range(E_min, E_max+1):
    logG[E] = 1
    Hist[E] = 0
statevector = np.random.rand(60)
for i in range(60):
    if statevector[i] > 0.5:
        statevector[i] = 1
    else:
        statevector[i] = -1
energy = 0


start = time.time()
while(True):
    f = f/2
    print(f)
    if f < 10e-6:
        break
    for E in range(E_min, E_max+1):
        Hist[E] = 0
    while(True):
        j = np.random.randint(0, 60)
        vector1 = list(statevector)
        vector2 = list(statevector)
        vector2[j] = -vector1[j]
        energy1 = energy_per_config(structure, vector1)
        energy2 = energy_per_config(structure, vector2)
        gap = np.exp(logG[energy1] - logG[energy2])
        # print(gap)
        if gap > 1:
            logG[energy2] += f
            Hist[energy2] += 1
            statevector = vector2
        else:
            p = np.random.rand()
            if p < gap:
                logG[energy2] += f
                Hist[energy2] += 1
                statevector = vector2
            else:
                logG[energy1] += f
                Hist[energy1] += 1
                
        error = np.std(np.array(list(Hist.values()))) / np.mean(np.array(list(Hist.values()))) - 1
        # print("ERROR", error)
        if error < 0.05:
            break
print(logG)

X = list(range(-66, 91))
Y = np.zeros(91+66)
for i in range(91+66):
    Y[i] = logG[X[i]]

Y = np.exp(Y - 479.88743591308594)
Y = [Y[i] for i in np.arange(0, 66+90, 2)]
# print(Y[0]/sum(Y)* 2**60)
E = (np.arange(-66, 90, 2))
Z = np.dot(Y, np.exp(-E))/sum(Y)*2**60
stop = time.time()

print(np.log(Z)/60, stop-start)