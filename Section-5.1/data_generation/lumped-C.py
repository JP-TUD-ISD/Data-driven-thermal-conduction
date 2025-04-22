import numpy as np

L = [0.04, 0.14, 0.22, 11.44 + 0.16]  # domain widths in m 
# important: In the last layer, 0.16 m are from FSS not subsoil
# this is a simplyfiyng assumption, that the 0.16 m can be neglected for this data generation
# the data set for layer 3 needs to be also applied for layer 4 in the simulation later on
C = [310.0 * 2400, 300.0 * 2400, 810.0 * 2400, 970.0 * 2400]  # heat capacity per domain 

K = 0.8  # thermal conductivity for one domain ... W / (m K)
# here it is assumed that just for one domain the conductivity is known
idK = 0  # identifier to highlight which domain is the one with known K
# remark: Python counting, so first material -> idk = 0

tbottom = 10.0  # temperature at 12m, see reference

sensorfile = "sensor.dat"  # file with real world data

a = np.loadtxt(sensorfile)

times = len(a[:, 0])  # get length of data set divided by 2

b = np.zeros((times * 2, len(L) + 1))  # plus one for so that the time can be saved as well

b[0 : times,       0] = a[:, 0]  # array storing the temperature gradient
b[times:2 * times, 0] = a[:, 0]

for i in range(len(L)):
    print(i)
    if i == (len(L) - 1):
        b[0:times        , i + 1] = (tbottom - a[:, i + 1]) / L[i]
        b[times:2 * times, i + 1] = -b[0:times        , i + 1]  # negative gradient is also valid, so double the data set size
    else:
        b[0:times, i + 1] = (a[:, i + 2] - a[:, i + 1]) / L[i]
        b[times:2 * times, i + 1] = -b[0:times        , i + 1]  # negative gradient is also valid, so double the data set size

d = np.zeros((times * 2, len(L) + 2))  # plus two for time, and layer 4, which now is a seperate data set
d[:, 0] = b[:, 0]
d[:, 1] = b[:, 1]
d[:, 2] = b[:, 2]
d[:, 3] = b[:, 3]
d[:, 4] = b[:, 3]
d[:, 5] = b[:, 4]

np.savetxt(fname="g.dat",X=d)

if idK > -1:
    c = np.zeros((times * 2, len(L) + 1))  # array storing the heat flux, plus one, since first entry is time
    c[:, 0] = b[:, 0]  # store time
    c[:, idK + 1] = K * b[:, idK + 1]
    for i in range(1, len(L) + 1):
        if i  == idK + 1:
            continue
        CH = (C[i - 1] * L[i - 1] + C[i - 2] * L[i - 2]) / 2.0
        for j in range(times):
            if j == 0:
                c[j, i] = c[j, idK + 1]
                c[j + times, i] = c[j + times, idK + 1]
            else:
                dt = a[j, 0] - a[j - 1, 0]  # first entry is time
                Tn1 = a[j, i]
                Tn  = a[j - 1, i]
                c[j, i] = c[j, i - 1] - CH * (Tn1 - Tn) / dt
                c[j + times, i] = c[j + times, i - 1] - CH * (Tn1 - Tn) / dt
    
    e = np.zeros((times * 2, len(L) + 2))  # array storing the heat flux, plus one, since first entry is time , plus one for layer 4
    e[:, 0] = c[:, 0]
    e[:, 1] = c[:, 1]
    e[:, 2] = c[:, 2]
    e[:, 3] = c[:, 3]
    e[:, 4] = c[:, 3]
    e[:, 5] = c[:, 4]
    np.savetxt(fname="q.dat",X=e)
