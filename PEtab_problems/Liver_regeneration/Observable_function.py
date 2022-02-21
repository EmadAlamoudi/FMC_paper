import sys
import numpy as np

def additive_noise(sumstat: dict = {}):
    noisy_sumstat = {}
    sumstat_edit = prepare_data(sumstat)
    sigma = {'condition1__YAP_nuclear_observable': 0.061763933333333,
             'condition1__YAP_total_observable': 0.050105066666667}
    for key, val in sumstat_edit.items():
        if key == 'loc': continue
        if key == 'condition1__YAP_nuclear_observable':
            noisy_sumstat[key] = val + sigma[key] * np.random.randn(len(val))
        else:
            noisy_sumstat[key] = val + sigma[key] * np.random.randn(len(val))

    return noisy_sumstat


def main(sumstat):
    noisy_sumstat = additive_noise(sumstat)
    new_sumstat = {}
    for key, val in noisy_sumstat.items():
        if key == "loc":
            continue
        new_sumstat[key] = noisy_sumstat[key]
    return new_sumstat

def prepare_data(sim):
    unique_obs = [0,8,15,20,30,50]
    step = 10
    itemindex = []
    new_dict = dict()
    uniqe_key_list = []
    for key in [*sim]:
        if key in ['loc', 'condition1__time', 'condition1__Tension', 'condition1__cell.id']:
            continue
        new_dict[key]=[]
        for index in unique_obs:
            new_dict[key].extend(sim[key][index*step:(index*step)+step])
    for key, val in new_dict.items():
        new_dict[key] = np.array(val)
    return new_dict

if __name__ == "__main__":
    sumstat = sys.argv[1]
    main(sumstat)
