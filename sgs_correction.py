import numpy as np
import matplotlib.pyplot as plt

np.random.seed(31)

s = np.random.rand(3,3)
s = np.dot(np.transpose(s), s) # sgs tensor
#print('SGS stress tensor is\n', s)
r = np.random.rand(3,3)
r = np.dot(np.transpose(r), r)
#print('Resolved stress tensor is\n', r)

# get eigenvalues of resolved stress tensor
rij = r/np.trace(r) - 1/3*np.eye(3)
[wr, Vr] = np.linalg.eig(rij)
idr = wr.argsort()[::-1]
wr = wr[idr]
Vr = Vr[idr]
#print(wr)

x1c = np.matrix([ 0.0 , 0.0 ])
x2c = np.matrix([ 1.0 , 0.0 ])
x3c = np.matrix([ 0.5 , np.sqrt(3.0)/2.0 ])

# barycentric coordinates of resolved stress tensor
xbcr = x1c*(wr[0]-wr[1]) + 2*x2c*(wr[1]-wr[2]) + x3c*(3*wr[2]+1)

# compute new eigenvalues of sgs tensor

sij = s/np.trace(s) - 1/3*np.eye(3)
[ws, Vs] = np.linalg.eig(sij)
ids = ws.argsort()[::-1]
ws = ws[ids]
Vs = Vs[ids]
print(sij)

xbcs = x1c*(ws[0]-ws[1]) + 2*x2c*(ws[1]-ws[2]) + x3c*(3*ws[2]+1)

deltaB = 0.5
xnew = xbcs[0,0] + deltaB*(xbcr[0,0] - xbcs[0,0])
ynew = xbcs[0,1] + deltaB*(xbcr[0,1] - xbcs[0,1])

l1 = (1 - deltaB)*ws[0] + deltaB*ws[0]
l2 = (1 - deltaB)*ws[1] + deltaB*ws[1]
l3 = (1 - deltaB)*ws[2] + deltaB*ws[2]

# assemble new sgs tensor
wst = np.eye(3)
wst[0,0] = l1
wst[1,1] = l2
wst[2,2] = l3

rstar = np.dot(Vs, wst)
rstar = np.dot(rstar, np.transpose(Vs))
print(rstar)



plt.scatter(xbcr[0,0], xbcr[0,1])
plt.scatter(xbcs[0,0], xbcs[0,1])
plt.scatter(xnew, ynew)
plt.plot([1, 0.5], [0, 0.866], 'k-')
plt.plot([0, 0.5], [0, 0.866], 'k-')
plt.plot([0, 1], [0, 0], 'k-')
plt.plot([0.5, 1.01, 0.1], [-0.01, 0.5, 0.88], 'w.') 
plt.text(0-.05,0-.05,'1c')
plt.text(1,0-.05,'2c')
plt.text(0.5-.025,0.866+.025,'3c')
#plt.show()
