import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
import matplotlib as mpl

'''
plotting script that plots the distance for each element over time and cooridante
x -> contains coordinate of the 1-D problem
y -> is time at which a simulation is performed
z -> distance
'''

mpl.style.use("./WileyNJDv5.cls")
mpl.rc("text", usetex=True)  # make plots with latex fonts
plt.rc('text.latex', preamble=r"\usepackage{siunitx}" #r"\documentclass[APA,Times1COL]{WileyNJDv5}"
                              r"\sisetup{per-mode=fraction}"
                              r"\usepackage{amsmath}")
mpl.rc("font", family="serif")
mpl.rc("axes", grid="true")
mpl.rc("grid", color="0.8", linestyle='dashed', linewidth=0.5)
mpl.rc("lines", linewidth=2)
mpl.rc("font", size="12")


# geometrical parameters
L = [0.04, 0.14, 0.22, 0.16, 11.44]  # domain widths in m
NN = [3, 12, 14, 15, 16]  # number of elements per domain
labels = ["SMA", "AC", "FSS$_1$", "FSS$_2$", "subsoil"]
dt = 353.14634458378345
fnam = "00_distances"
usee = True
uselog = True
usesurf = True

cmi = 1/2.54  # centimeters in inches

N = 0  # total number of elements used, is calculated
NNN = np.zeros(len(NN))  # assignment which element numbers to which domain
count = 0
for i in NN:
    N += i
    NNN[count] = N
    count += 1

LL = 0.0  # total length
for i in L:
    LL += i

x = np.zeros(N + 1)
xe = np.zeros(N + 1)
count = 1  # zeroth entry of x is zero, so start at one
for i in range(len(NN)):
    Le = L[i] / NN[i]
    for j in range(int(NN[i])):
        x[count] = Le + x[count - 1]  # just add the element length to previous coordinate
        xe[count] = count
        count += 1
        
Z = np.loadtxt(fnam)
ZZ = np.log10(Z[1:,:])
maZZ = np.max(ZZ)
miZZ = np.min(ZZ)
a1 = np.shape(Z)[0]
y = np.zeros(a1 - 1)
for i in range(0, a1 - 1):
    y[i] = ((i + 1) * dt)*1.e-6

if usee:
    X, Y = np.meshgrid(xe[:-1], y)
else:
    X, Y = np.meshgrid(xe[:-1], y)

fig = plt.figure(figsize=(15 * cmi, 18.5 * cmi))
ax = plt.axes(projection='3d')

if usesurf:
    if uselog:
        surf = ax.plot_surface(X,Y,ZZ, cmap=cm.coolwarm, rcount=100, ccount=100)
    else:
        surf = ax.plot_surface(X,Y,Z[1:,:], cmap=cm.coolwarm)
else:
    if uselog:
        if usee:
            for i in range(len(y)):
                ax.scatter(xe[1:], ZZ[i, :], zs=y[i], zdir="y", c=cm.coolwarm( (ZZ[i, :]-miZZ)/(maZZ - miZZ) ), edgecolor="none")
        else:
            for i in range(len(y)):
                ax.scatter(x[1:], ZZ[i, :], zs=y[i], zdir="y", cmap=cm.coolwarm( (ZZ[i, :]-miZZ)/(maZZ - miZZ) ), edgecolor="none")
    else:
        print("fuck")
        exit()

if usee:
    ax.set_xlim(0, xe[-1])
else:
    ax.set_xlim(0, LL)

ax.set_ylim(0, 6.04727800e+0)
#ax.set_ylim(0, 6.0)
if uselog:
    ax.set_zlim(-12, 0.5)
else:
    ax.set_zlim(0, 8)

if usee:
    ax.set_xlabel("material")
else:
    ax.set_xlabel("depth")
ax.set_ylabel("time in 10$^6$ s")

if uselog:
    ax.set_zlabel("logarithm of distance")
else:
    ax.set_zlabel("distance")

# ax.set_xticks(list(ax.get_xticks()) + NN)
ax.set_xticks(NNN, labels=labels)

if uselog:
    ax.view_init(elev=6., azim=-90.)
else:
    ax.view_init(elev=13., azim=-70.)

# ax.set_zscale('log')

# fig.colorbar(surf, shrink=0.5, aspect=15)

plt.savefig("test.png", dpi=600)
plt.savefig("distance-over-geometry-100.svg")

plt.show()
