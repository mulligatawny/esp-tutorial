import numpy as np

a = np.array([[1,2],[2,3]])

v, w = np.linalg.eigh(a)

print(v)
print(w)

print(a)
print(w @ v*np.eye(2) @ np.transpose(w))
