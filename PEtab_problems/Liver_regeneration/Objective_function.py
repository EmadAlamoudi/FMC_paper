import numpy as np


def eucl_dist(sim, obs):
    if sim == -15:
        print("timeout")
        return np.inf

    total = 0
    for key in sim:
        if key in ('loc', "condition1__time", "condition1__cell.id", "condition1__Tension"):
            continue

        x = np.array(sim[key])
        y = np.array(obs[key])
        z = np.array(obs[key + "_SEM"])
        # simulation does not finish successfuly, only partial part of the
        # result wrtten. In such case, ignore the parameter vector
        if x.size != y.size:
            # size does not match
            return np.inf
        # if np.max(y) != 0:
        #     x = x/np.max(y)
        #     y = y/np.max(y)
        total += np.sum(((x - y)/z) ** 2)
    print("total: ", total)
    return total
