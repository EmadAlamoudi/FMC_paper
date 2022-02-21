import os
from pyabc.external.r import R
import pyabc


def obj_func_old(sumstat: dict = {}):
    cur_path = os.getcwd()
    r = R(os.path.join(cur_path, 'abc_r_blueprint.R'))

    # Target distances
    dist_speed_tar = r.distance("distance_speed_tar")
    dist_trn_tar = r.distance("distance_trn_tar")
    dist_arr_tar = r.distance("distance_arr_tar")
    dist_str_tar = r.distance("distance_str_tar")
    dist_msd_tar = r.distance("distance_msd_tar")

    # Infected distances
    dist_speed_inf = r.distance("distance_speed_inf")
    dist_trn_inf = r.distance("distance_trn_inf")
    dist_arr_inf = r.distance("distance_arr_inf")
    dist_str_inf = r.distance("distance_str_inf")
    dist_msd_inf = r.distance("distance_msd_inf")

    # Aggregate distances
    distance = pyabc.AggregatedDistance([dist_speed_tar, dist_trn_tar,
                                         dist_arr_tar,
                                         dist_str_tar,
                                         dist_msd_tar,
                                         dist_speed_inf,
                                         dist_trn_inf,
                                         dist_arr_inf,
                                         dist_str_inf,
                                         dist_msd_inf])
    return distance


def obj_func(x, x_0):
    cur_path = os.getcwd()
    r = R(os.path.join(
        os.path.join(os.path.dirname(os.path.dirname(cur_path)), 'HIV'),
        'abc_r_blueprint.R'))
    # Target distances
    dist_speed_tar = r.distance("distance_speed_tar")
    dist_trn_tar = r.distance("distance_trn_tar")
    dist_arr_tar = r.distance("distance_arr_tar")
    dist_str_tar = r.distance("distance_str_tar")
    dist_msd_tar = r.distance("distance_msd_tar")

    # Infected distances
    dist_speed_inf = r.distance("distance_speed_inf")
    dist_trn_inf = r.distance("distance_trn_inf")
    dist_arr_inf = r.distance("distance_arr_inf")
    dist_str_inf = r.distance("distance_str_inf")
    dist_msd_inf = r.distance("distance_msd_inf")

    # Aggregate distances
    distance = pyabc.AggregatedDistance([dist_speed_tar, dist_trn_tar,
                                         dist_arr_tar,
                                         dist_str_tar,
                                         dist_msd_tar,
                                         dist_speed_inf,
                                         dist_trn_inf,
                                         dist_arr_inf,
                                         dist_str_inf,
                                         dist_msd_inf])
    return distance(x,x_0)
