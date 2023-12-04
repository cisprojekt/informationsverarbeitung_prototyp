import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

def mds(distance_matrix, dimensions=2):
    """
    Does multidimensional scaling (also for non-euclidian distances)
    Will only use positive eigenvalues according to
    https://mediaspace.nottingham.ac.uk/media/MDSA+Non-Euclidean+distance+matrices/1_nsh2vdgr
    """
    
    n = distance_matrix.shape[0]

    # Create centering matrix H
    C = np.eye(n) - np.ones((n, n))/n

    # Double centering
    DC = -C.dot(distance_matrix**2).dot(C)/2

    # Get eigenvalues and eigenvectors
    evals, evecs = eigh(DC)

    # Sort eigenvalues and -vectors
    idx   = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]

    # Only use eigenvalues larger than 0
    w, = np.where(evals > 0)
    
    # Create diagonal Matrix from square roots of positive eigenvalues
    L  = np.diag(np.sqrt(evals[w]))
    V  = evecs[:,w]
    M  = V.dot(L)
    
    return M[:,:dimensions]

def is_euclidean(distance_matrix):
    """
    Checks if matrix is euclidean by testing for
    triangle inequality
    """
    n = distance_matrix.shape[0]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if distance_matrix[i, j] > distance_matrix[i, k] + distance_matrix[k, j]:
                    return False
    return True

# Example distance matrix
D = np.array([
    [0, 3, 8, 2],
    [3, 0, 2, 1],
    [8, 2, 0, 3],
    [2, 1, 3, 0]
])

print("Matrix is euclidean: {}".format(is_euclidean(D)))

M = mds(D)

# Plot the points
plt.figure(figsize=(8, 6))
plt.scatter(M[:, 0], M[:, 1], color="blue")
for i in range(M.shape[0]):
    plt.text(M[i, 0], M[i, 1], str(i+1), color="red")
plt.title("MDS of the given distance matrix")
plt.grid(True)
plt.show()
