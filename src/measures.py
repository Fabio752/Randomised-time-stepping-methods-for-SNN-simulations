import numpy as np


def mse(coherences, ground_truth):
    return np.sum(np.square(coherences - ground_truth))


def sync_measure(V):
    # Initialise lists.
    v_i_squared = np.square(V)
    avg_potentials = np.mean(V, axis = 0)

    # Compute delta_n
    squared_potentials = np.square(avg_potentials)
    delta_n = np.mean(squared_potentials) - np.square(np.mean(avg_potentials))
    
    # Compute delta
    mean_v_i_squared = np.square(np.mean(V, axis = 1))
    delta = np.mean(np.mean(v_i_squared, axis = 1) - mean_v_i_squared)
    
    # Compute sigma_n
    sigma_n = delta_n / delta 
    return sigma_n
