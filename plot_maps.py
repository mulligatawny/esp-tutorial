###############################################################################
# Testbed for Eigenspace Subgrid Modeling - Barycentric Map Generator         #
# (run this script after DCP.cpp)                                             #
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# load initial and perturbed sgs anisotropy tensor
A_sgs = np.loadtxt('A_sgs.dat').reshape(3,3)
A_sgsp = np.loadtxt('A_sgsp.dat').reshape(3,3)
A_res = np.loadtxt('A_res.dat').reshape(3,3)

def get_bcc(A):
    [w, V] = np.linalg.eig(A)
    # sort eigenvalues from largest to smallest
    idx = w.argsort()[::-1]  
    w = w[idx]
    V = V[idx]
    # compute barycentric coordinates (17)
    x1c = np.matrix([ 0.0 , 0.0 ])
    x2c = np.matrix([ 1.0 , 0.0 ])
    x3c = np.matrix([ 0.5 , np.sqrt(3.0)/2.0 ])
    xbc = x1c*(w[0]-w[1]) + 2*x2c*(w[1]-w[2]) + x3c*(3*w[2]+1)
    return xbc

xbc = get_bcc(A_sgs)
xbcp = get_bcc(A_sgsp)
xbcr = get_bcc(A_res)

# plot anisotropy on barycentric map (Figure 1)
fig = plt.figure()
plt.scatter(xbc[0,0], xbc[0,1], label='original sgs')
plt.scatter(xbcp[0,0], xbcp[0,1], label='perturbed sgs')
plt.scatter(xbcr[0,0], xbcr[0,1], label='resolved')
plt.plot([xbc[0,0], xbcr[0,0]],[xbc[0,1], xbcr[0,1]], 'k', linewidth=0.5)
plt.legend()
plt.plot([1, 0.5], [0, 0.866], 'k-')
plt.plot([0, 0.5], [0, 0.866], 'k-')
plt.plot([0, 1], [0, 0], 'k-')
plt.plot([0.5, 1.01, 0.1], [-0.01, 0.5, 0.88], 'w.')
plt.text(0-.05,0-.05,'1c')
plt.text(1,0-.05,'2c')
plt.text(0.5-.025,0.866+.025,'3c')
plt.show()
