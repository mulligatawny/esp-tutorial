###############################################################################
# Barycentric map for turbulence anisotropy (Jofre et al. 2017)               #
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

S_sgs = np.array([[2,-1,0],[-1,2,-1],[0,-1,2]])
S_res = np.array([[24, 14, 4],[14, 13, 17],[4, 17, 51]])

# compute normalized anisotropy tensor (13)
A_sgs = S_sgs/np.trace(S_sgs) - 1/3*np.eye(3)
A_res = S_res/np.trace(S_res) - 1/3*np.eye(3)

#print('-----------------------------------------------------------------------')
#print('The resolved anisotropy tensor is:\n', A_res)
print('The sgs anisotropy tensor is:\n', A_sgs)

def get_eigs(aij):
    [w, V] = np.linalg.eig(aij)
    # sort eigenvalues from largest to smallest
    idx = w.argsort()[::-1]  
    w = w[idx]
    V = V[idx]
    return w, V
    
w_sgs, V_sgs = get_eigs(A_sgs)
w_res, V_res = get_eigs(A_res)

#print('The resolved eigenvalues are\n', w_res)
#print('The resolved eigenvectors are\n', V_res)
#print('The sgs eigenvalues are\n', w_sgs)
#print('The sgs eigenvectors are\n', V_sgs)

# perturb w_sgs toward w_res
w_sgsp = np.zeros_like(w_sgs)

dB = 0.789

w_sgsp[0] = (1.0 - dB)*w_sgs[0] + dB*w_res[0]
w_sgsp[1] = (1.0 - dB)*w_sgs[1] + dB*w_res[1]
w_sgsp[2] = (1.0 - dB)*w_sgs[2] + dB*w_res[2]

# form perturbed stress
#S_sgsp = - V_sgs.dot(w_sgsp*np.eye(3)).dot(np.transpose(V_sgs))
S_sgsp = np.transpose(V_sgs).dot(w_sgsp*np.eye(3)).dot(V_sgs)

print('The perturbed sgs anisotropy tensor is\n', S_sgsp)
