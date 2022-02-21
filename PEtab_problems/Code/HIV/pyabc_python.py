# -*- coding: utf-8 -*-
import pyabc
import numpy as np
from pyabc.external.r import R
from pyabc import Distribution, RV, ABCSMC
from pyabc import populationstrategy
from pyabc.distance.scale import median
import fitmulticell as fmc
import os
from pyabc.sampler import RedisEvalParallelSampler
import argparse
import glob
import pandas as pd
from fitmulticell.sumstat import SummaryStatistics as ss

# Password for the redis_server

# Morpheus Model
cur_path = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/HIV'
# cur_path = "/home/emad/Documents/HIV"
filename = "HIV_2D_model.xml"
file_ = os.path.join(cur_path, filename)
par_map = {'PS_tar': './Global/Constant[@symbol="PS_tar"]',
           'DT_tar': './Global/Constant[@symbol="DT_tar"]',
           'PS_inf': './Global/Constant[@symbol="PS_inf"]',
           'DT_inf': './Global/Constant[@symbol="DT_inf"]',
           'JTT': './Global/Constant[@symbol="JTT"]',
           'JTI': './Global/Constant[@symbol="JTI"]',
           'JIM': './Global/Constant[@symbol="JIM"]',
           'JTM': './Global/Constant[@symbol="JTM"]',
           'JTC': './Global/Constant[@symbol="JTC"]',
           'JIC': './Global/Constant[@symbol="JIC"]',
           'JII': './Global/Constant[@symbol="JII"]',
           'C_scs': './Global/Constant[@symbol="C_scs"]',
           'C_vcs': './Global/Constant[@symbol="C_vcs"]'}

r = R(os.path.join(cur_path, 'abc_r_blueprint.R'))
sum_stat = r.summary_statistics("mySummaryStatistics", is_py_model=True) # If pyabc ver 0.9.22 add is_py_model=True


def sum_stat_wrapper(loc):
    print("======= sumstat =======")
    print("loc: ", loc)
    return sum_stat({"loc": loc})


sumstat = sum_stat_wrapper
model = fmc.model.MorpheusModel(model_file = file_,
                                sumstat=sumstat,
                                executable="/home/emad/morpheus-2.2.5",
                                par_map=par_map,
                                show_stdout=False,
                                show_stderr=False,
                                raise_on_error=False)
# R - Distances and sumstats
par_obs = {'PS_tar': 30,
           'DT_tar': 15,
           'PS_inf': 30,
           'DT_inf': 15,
           'JTT': 400.74,
           'JTI': 279.90,
           'JIM': 368.05,
           'JTM': 311.32,
           'JTC': 364.56,
           'JIC': 442.31,
           'JII': 66.38,
           'C_scs': 145.86,
           'C_vcs': 201.05}

model.sample(par_obs)
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
distance = pyabc.AggregatedDistance([dist_speed_tar,
                                     dist_trn_tar,
                                     dist_arr_tar,
                                     dist_str_tar,
                                     dist_msd_tar,
                                     dist_speed_inf,
                                     dist_trn_inf,
                                     dist_arr_inf,
                                     dist_str_inf,
                                     dist_msd_inf])

# Sumstats


def distance_func():
    print("======= distance =======")
    return distance


dt_domain = np.arange(60)

prior = Distribution(PS_tar=RV("uniform", 0, 100),  # low,low+interval
                     DT_tar=RV("rv_discrete", values=(dt_domain, [1/60]*60)),
                     PS_inf=RV("uniform", 0, 100),
                     DT_inf=RV("rv_discrete", values=(dt_domain, [1/60]*60)),
                     JTT=RV("uniform", 0, 500),  
                     JTI=RV("uniform", 0, 500),  
                     JIM=RV("uniform", 0, 500),
                     JTM=RV("uniform", 0, 500),
                     JTC=RV("uniform", 0, 500), 
                     JIC=RV("uniform", 0, 500),
                     JII=RV("uniform", 0, 500),
                     C_scs=RV("uniform", 0, 100),
                     C_vcs=RV("uniform", 0, 100)) 


transition = pyabc.transition.AggregatedTransition(mapping={
    'PS_tar': pyabc.MultivariateNormalTransition(),
    'DT_tar': pyabc.DiscreteJumpTransition(domain=dt_domain, p_stay=0.7),
    'PS_inf': pyabc.MultivariateNormalTransition(),
    'DT_inf': pyabc.DiscreteJumpTransition(domain=dt_domain, p_stay=0.7),
    'JTT': pyabc.MultivariateNormalTransition(),
    'JTI': pyabc.MultivariateNormalTransition(),
    'JIM': pyabc.MultivariateNormalTransition(),
    'JTM': pyabc.MultivariateNormalTransition(),
    'JTC': pyabc.MultivariateNormalTransition(),
    'JIC': pyabc.MultivariateNormalTransition(),
    'JII': pyabc.MultivariateNormalTransition(),
    'C_scs': pyabc.MultivariateNormalTransition(),
    'C_vcs': pyabc.MultivariateNormalTransition()})



# ABCSMC
pop_strategy = pyabc.ConstantPopulationSize(nr_particles=2)
eps = pyabc.QuantileEpsilon(65, alpha=0.7)  # Defines which distances is needed in next populatio
abc = ABCSMC(model, prior, distance,
             transitions=transition, population_size=pop_strategy, eps=eps)

obs = r.observation("Exp")

# tryjectory = model.sample(par_map_sample)
# distance(tryjectory, r.observation("Exp"))
# Database
db_name = str.split(filename, sep = '.')[0]
db = "sqlite:///" + os.path.join(cur_path, db_name + '.db') # Define name of the database

abc.new(db, r.observation("Exp"))
#abc.load(db, 1, r.observation("Exp"))

# Abbruchkriterium
history = abc.run(max_nr_populations=2)
