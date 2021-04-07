###############################################################################
# Barycentric map for turbulence anisotropy (Jofre et al. 2017)               #
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

state = 4
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
    S_sgs = np.array([[2,-1,0],[-1,2,-1],[0,-1,2]])
    S_res = np.array([[24, 14, 4],[14, 13, 17],[4, 17, 51]])

# compute normalized anisotropy tensor (13)
A_sgs = S_sgs/np.trace(S_sgs) - 1/3*np.eye(3)
A_res = S_res/np.trace(S_res) - 1/3*np.eye(3)

print('Original sgs stress\n', A_sgs)

def get_eigs(aij):
    [w, V] = np.linalg.eig(aij)
    # sort eigenvalues from largest to smallest
    idx = w.argsort()[::-1]  
    w = w[idx]
    V = V[idx]
    return w, V
    
w_sgs, V_sgs = get_eigs(A_sgs)
w_res, V_res = get_eigs(A_res)

#print(V_sgs)

# compute barycentric coordinates (17)
def get_bcxy(w):
    x1c = np.matrix([ 0.0 , 0.0 ])
    x2c = np.matrix([ 1.0 , 0.0 ])
    x3c = np.matrix([ 0.5 , np.sqrt(3.0)/2.0 ])
    xbc = x1c*(w[0]-w[1]) + 2*x2c*(w[1]-w[2]) + x3c*(3*w[2]+1)
    return xbc

x_sgs = get_bcxy(w_sgs)
x_res = get_bcxy(w_res)



# perturb w_sgs toward w_res
w_sgsp = np.zeros_like(w_sgs)

dB = 1.0

w_sgsp[0] = (1 - dB)*w_sgs[0] + dB*w_res[0]
w_sgsp[1] = (1 - dB)*w_sgs[1] + dB*w_res[1]
w_sgsp[2] = (1 - dB)*w_sgs[2] + dB*w_res[2]

# form perturbed stress

#S_sgsp = V_sgs @ w_sgsp*np.eye(3) @ np.transpose(V_sgs)
S_sgsp = - V_sgs.dot(w_sgsp*np.eye(3)).dot(np.transpose(V_sgs))

print('Perturbed stress is\n', S_sgsp)

# plot anisotropy on barycentric map (Figure 1)
fig = plt.figure()
plt.scatter(x_sgs[0,0], x_sgs[0,1])
plt.scatter(x_res[0,0], x_res[0,1])
plt.plot([1, 0.5], [0, 0.866], 'k-')
plt.plot([0, 0.5], [0, 0.866], 'k-')
plt.plot([0, 1], [0, 0], 'k-')
plt.plot([0.5, 1.01, 0.1], [-0.01, 0.5, 0.88], 'w.') 
plt.text(0-.05,0-.05,'1c')
plt.text(1,0-.05,'2c')
plt.text(0.5-.025,0.866+.025,'3c')
#plt.show()
