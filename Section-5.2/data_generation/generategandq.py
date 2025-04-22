import numpy as np

Tempfile = "T.dat"
Qfile = "Q.dat"

gfile = "g.dat"
qfile = "q.dat"

# geometrical parameters
L = [0.04, 0.14, 0.22, 0.16, 11.44]  # domain widths in m
NN = [10, 20, 30, 20, 150]  # number of elements per domain

LE = np.zeros(len(L))

for i in range(len(L)):
    LE[i] = L[i] / NN[i]  # get length of one element
    
N = 0
NNN = np.zeros(len(NN))  # assignment which element numbers to which domain
count = 0
for i in NN:
    N += i
    NNN[count] = N
    count += 1

Ts = np.loadtxt(Tempfile)
Qs = np.loadtxt(Qfile)

steps = np.shape(Ts)[0]

g = np.zeros((steps, len(L) + 1))
q = np.zeros((steps, len(L) + 1))

g[:, 0] = Ts[:, 0]
q[:, 0] = Qs[:, 0]

for i in range(len(L)):
    if i == 0:  # get q sensor in middle
        idq = int(NN[i] / 2)
    else:
        idq = int((NN[i] + NN[i - 1]) / 2)
    idg1 = idq - 1  # sensor left of q sensor
    idg2 = idq + 1  # sensor right of q sensor
    g[:, i + 1] = (Ts[:, idg2] - Ts[:, idg1]) / (2.0 * LE[i])
    q[:, i + 1] = Qs[:,idq]

np.savetxt(fname=gfile, X=g)
np.savetxt(fname=qfile, X=q)
