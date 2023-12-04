#!/usr/bin/env python
# coding: utf-8

"""
Implementation of the SMACOF-Algorithm for MDS as it is described in
[1] Modern Multidimensional Scaling. (2005). In Springer Series in Statistics.
Springer New York. https://doi.org/10.1007/0-387-28981-x
"""

from sklearn.manifold import MDS
from matplotlib import pyplot as plt
import numpy as np
import random


def create_rand_3d_data(n):
    """
    Create n random points with x-, y- and z-coordinates between 0 and 100
    """
    points = []
    for i in range(n):
        points.append(
            [random.randint(0, 100), random.randint(0, 100), random.randint(0, 100)]
        )

    points = np.array(points)
    return points


def distance(point_a, point_b):
    """
    Calculate Euclidean distance between two points
    """
    return np.linalg.norm(point_a - point_b)


def distance_matrix(points):
    """
    Calculate distance matrix of size n x n for given list of points
    """
    n = len(points)
    dist_mat = np.ndarray((n, n))
    for i, point_a in enumerate(points):
        for j, point_b in enumerate(points):
            dist_mat[i, j] = distance(point_a, point_b)

    return dist_mat


def create_plot(points):
    """
    Create a scatter-plot for a list of points where each
    point has an x- and y-coordinate
    """

    # Create list of x- and y-coordinates
    x = [point[0] for point in points]
    y = [point[1] for point in points]

    plt.scatter(x, y)

    # Set limit for consistent viewframe
    plt.xlim([-100, 100])
    plt.ylim([-100, 100])

    # Add numbering to points
    for i in range(len(x)):
        plt.annotate(i, (x[i] + 2, y[i]))

    plt.show()


def calculate_V(weights):
    """
    Calculate the matrix V according to eq. (8. 18)
    """
    V = np.ndarray(weights.shape)
    for iy, ix in np.ndindex(V.shape):
        if iy == ix:
            V[iy, ix] = sum(
                w if (j != ix) else 0 for j, w in enumerate(weights[..., ix])
            )
        else:
            V[iy, ix] = -weights[iy, ix]

    return V


def calculate_random_Z(dim_n, dim_m):
    """
    Calculate random data for Z as n x m Matrix
    """
    Z = np.ndarray((dim_n, dim_m))
    for iy, ix in np.ndindex(Z.shape):
        Z[iy, ix] = random.uniform(-1, 1)

    return Z


def calculate_weights(dist_mat):
    """
    Calculate weights, matrix entry is 1 if (dis)similarity is present
    in distance matrix, otherwise set 0
    """
    weights = np.zeros(dist_mat.shape)
    for iy, ix in np.ndindex(weights.shape):
        if dist_mat[iy, ix]:
            weights[iy, ix] = 1

    return weights


def calculate_B(Z, weights, dist_mat):
    """
    Calculate B according to e.q. (8.24)
    """
    D = distance_matrix(Z)
    B = np.ndarray(D.shape)

    # Calculate non-diagonal elemnts first
    for iy, ix in np.ndindex(B.shape):
        if iy != ix:
            if D[iy, ix] != 0:
                B[iy, ix] = -(weights[iy, ix] * dist_mat[iy, ix]) / (D[iy, ix])
            else:
                B[iy, ix] = 0

    # Now calculate diagonal elements
    for iy, ix in np.ndindex(B.shape):
        if iy == ix:
            B[iy, ix] = -sum(b if (j != ix) else 0 for j, b in enumerate(B[..., ix]))

    return B


def calculate_const(weights, dist_mat):
    """
    Calculate the first constant term according to e.q. (8.15)
    """
    const = 0
    for iy, ix in np.ndindex(weights.shape):
        if iy < ix:
            const += weights[iy, ix] * dist_mat[iy, ix] ** 2

    return const


def stress_function(X, V, Z, B, weights, dist_mat, verbose=False):
    """
    Calculate the stress according to e.q. (8.27)
    """
    term_1 = calculate_const(weights, dist_mat)
    term_2 = np.trace(X.T @ V @ X)
    term_3 = -2 * np.trace(X.T @ B @ X)

    stress = term_1 + term_2 + term_3

    return stress


def guttman_transform(n, B, Z, weights):
    """
    Calculate Guttman transformation according to e.q. (8.28, 8.29)
    """
    weights_one = True

    # Check all weights to select the correct formula
    for iy, ix in np.ndindex(weights.shape):
        if weights[iy, ix] != 1:
            weights_one = False
            break

    # If all weights are 1 use e.q. (8.29)
    # TODO: implement Moore-Penrose inverse if there exists weight unequal to 1
    X_u = (1 / n) * B @ Z

    return X_u


def calculate_mds(
    dist_mat, max_it=50, eps=10e-6, dim=2, book_example=False, verbose=False
):
    """
    Calculate a lower-dimensional representation of a distance-matrix
    using the SMACOF-Algorithm
    
    Args:
        dist_mat (np.ndarray): Input distance matrix (dissimilarities)
        max_it (int): Maximum number of iterations (default is 50)
        eps (float): Stop condition for stress difference (default is 10e-6)
    
    Returns:
        X_u (np.ndarray): List of x- and y-coordinates in low-dimensional space
    """

    # Counter for the iterations
    k = 0

    # Create weight matrix
    weights = calculate_weights(dist_mat)

    # Step 1: Set X_0 = Z, with random start configuration
    Z = calculate_random_Z(dist_mat.shape[0], dim)
    X = Z

    # Use example from book if book_example is set
    if book_example:
        Z = np.array(
            [[-0.266, -0.539], [0.451, 0.252], [0.016, -0.238], [-0.200, 0.524]]
        )
        X = Z

    # Step 2: Compute initial stress
    V = calculate_V(weights)
    B = calculate_B(Z, weights, dist_mat)
    stress = stress_function(X, V, X, B, weights, dist_mat)

    if verbose:
        print("initial stress: {}".format(stress))

    # Step 3: While-loop, stop if maximal iteration is reached
    while k < max_it:

        # Step 4: Increase iteration counter
        k += 1

        # Step 5: Compute Guttman Transformation
        B = calculate_B(Z, weights, dist_mat)
        X_u = guttman_transform(B.shape[0], B, Z, weights)

        # Step 6: Compute stress of updated X
        new_stress = stress_function(X_u, V, Z, B, weights, dist_mat)
        if verbose:
            print("iteration {} stress: {}".format(k, new_stress))

        # Step 7: Update Z
        Z = X_u

        # TODO: Step 8: End while

    return X_u


# (1) Use random distance matrix
dist_mat = distance_matrix(create_rand_3d_data(5))

# Using Scikit
embedding = MDS(n_components=2, dissimilarity="precomputed", normalized_stress=False)
conf_1 = embedding.fit_transform(dist_mat)

# Using own implementation
conf_2 = calculate_mds(dist_mat)

# Compare both results
create_plot(conf_1)
create_plot(conf_2)

# (2) Use example in book
dist_mat = [[0, 5, 3, 4], [5, 0, 2, 2], [3, 2, 0, 1], [4, 2, 1, 0]]
dist_mat = np.array(dist_mat)
conf_3 = calculate_mds(dist_mat, book_example=True, verbose=True)
create_plot(conf_3)
