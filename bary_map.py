###############################################################################
# Barycentric map for turbulence anisotropy (Jofre et al. 2017)               #
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

state = 5
if state == 1:
    T = np.eye(3)
    T[0,0] = 0.0
    T[1,1] = 0.0
elif state == 2:
    T = np.eye(3)
    T[0,0] = 0.0
elif state == 3:
    T = np.eye(3)
else:
    T = np.random.rand(3,3)
# generate 3x3 positive semi-definite matrix
A = np.dot(np.transpose(T),T)
# make trace elements equal and three orders of magnitude larger
A[0,0] = A[0,0]*1e3
A[1,1] = A[0,0]
A[2,2] = A[0,0]
print('The matrix is\n', A)
# compute normalized anisotropy tensor (13)
aij = A/np.trace(A) - 1/3*np.eye(3)
# eigendecomposition
[w, V] = np.linalg.eig(aij)
# sort eigenvalues from largest to smallest
idx = w.argsort()[::-1]  
w = w[idx]
V = V[idx]
print('The eigenvalues are',w)
print('The sum of the eigenvalues is {:.2f}'.format(np.sum(w)))

# compute barycentric coordinates (17)
x1c = np.matrix([ 0.0 , 0.0 ])
x2c = np.matrix([ 1.0 , 0.0 ])
x3c = np.matrix([ 0.5 , np.sqrt(3.0)/2.0 ])
xbc = x1c*(w[0]-w[1]) + 2*x2c*(w[1]-w[2]) + x3c*(3*w[2]+1)

# plot anisotropy on barycentric map (Figure 1)
fig = plt.figure()
plt.scatter(xbc[0,0], xbc[0,1])
plt.plot([1, 0.5], [0, 0.866], 'k-')
plt.plot([0, 0.5], [0, 0.866], 'k-')
plt.plot([0, 1], [0, 0], 'k-')
plt.plot([0.5, 1.01, 0.1], [-0.01, 0.5, 0.88], 'w.') 
plt.text(0-.05,0-.05,'1c')
plt.text(1,0-.05,'2c')
plt.text(0.5-.025,0.866+.025,'3c')
plt.show()
