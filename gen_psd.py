import numpy as np

A = np.array([[4,2,1],[-2,0,5],[2,3,5]])
A = np.dot(np.transpose(A), A)
print(A)
