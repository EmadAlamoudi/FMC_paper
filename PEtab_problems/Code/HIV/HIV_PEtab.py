import petab_MS
from fitmulticell.PEtab.base import PetabImporter
from fitmulticell.model import MorpheusModel as morpheus_model
from fitmulticell.model import MorpheusModels as morpheus_models
from pyabc.sampler import RedisEvalParallelSampler
from pyabc.sampler import MulticoreEvalParallelSampler
from pyabc import QuantileEpsilon
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import pyabc
import matplotlib.pylab as plt
from pathlib import Path
import os
import tempfile
from pyabc.external.r import R
from pyabc import Distribution, RV, ABCSMC

# import pathlib
#
# os.chdir('/home/emad/PycharmProjects/PEtab_FMC_extension')
# print(pathlib.Path().resolve())
# petab_problem_path = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/Liver_regeneration' + '/Meyer_MolSystBiol_2020.yaml'
# petab_problem = petab_MS.Problem.from_yaml(petab_problem_path)
# importer = PetabImporter(petab_problem)
# model = importer.import_petab_problem()

petab_problem_path = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/HIV' + '/HIV.yaml'
petab_problem = petab_MS.Problem.from_yaml(petab_problem_path)
importer = PetabImporter(petab_problem)
PEtab_prior = importer.create_prior()
par_map_imported = importer.get_par_map()
obs_pars_imported = petab_problem.get_x_nominal_dict(scaled=True)
PEtab_par_scale = petab_problem.get_optimization_parameter_scales()
dict_data_imported = petab_problem.get_measurement_dict()
PEtab_model = importer.create_model(executable="/home/emad/morpheus-2.2.5")

cur_path = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/HIV'
r = R(os.path.join(cur_path, 'abc_r_blueprint.R'))
sum_stat = r.summary_statistics("mySummaryStatistics", is_py_model=True) # If pyabc ver 0.9.22 add is_py_model=True


def sum_stat_wrapper(loc):
    return sum_stat({"loc": loc})
sumstat = sum_stat_wrapper

PEtab_model.sumstat = sumstat
# PEtab_tryjectory = PEtab_model.sample(obs_pars_imported)
obj_function = importer.get_objective_function()


dt_domain = np.arange(60)
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


pop_strategy = pyabc.ConstantPopulationSize(nr_particles=2)
eps = pyabc.QuantileEpsilon(65, alpha=0.7)  # Defines which distances is needed in next populatio
abc = ABCSMC(PEtab_model, PEtab_prior, obj_function,
             transitions=transition, population_size=pop_strategy, eps=eps)


filename = "HIV_2D_model.xml"
db_name = str.split(filename, sep='.')[0]
cur_path = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/HIV'

db = "sqlite:///" + os.path.join(cur_path, db_name + '.db')

abc.new(db, dict_data_imported["condition1"])   
#abc.load(db, 1, r.observation("Exp"))

# Abbruchkriterium
history = abc.run(max_nr_populations=2)

print("Done!")
