import numpy as np


def eucl_dist(sim, obs):
    total = 0
    for key in sim:
        if key in 'loc':
            continue
        x = np.array(sim[key])
        y = np.array(obs[key])
        # z = np.array(obs["SEM_" + key])
        if x.size != y.size:
            print("size not match")
            return np.inf

        total += np.sum((x - y) ** 2)
    return total

