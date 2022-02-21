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




# import pathlib
#
# os.chdir('/home/emad/PycharmProjects/PEtab_FMC_extension')
# print(pathlib.Path().resolve())
# petab_problem_path = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/Liver_regeneration' + '/Meyer_MolSystBiol_2020.yaml'
# petab_problem = petab_MS.Problem.from_yaml(petab_problem_path)
# importer = PetabImporter(petab_problem)
# model = importer.import_petab_problem()

# petab_problem_path = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/Liver_regeneration' + '/Meyer_MolSystBiol_2020.yaml'
# petab_problem = petab_MS.Problem.from_yaml(petab_problem_path)
# importer = PetabImporter(petab_problem)
# PEtab_prior = importer.create_prior()
# par_map_imported = importer.get_par_map()
# obs_pars_imported = petab_problem.get_x_nominal_dict(scaled=True)
# PEtab_par_scale = petab_problem.get_optimization_parameter_scales()
# dict_data_imported = petab_problem.get_measurement_dict()
# PEtab_model = importer.create_model()
# PEtab_model.timeout = 900
# PEtab_model.ignore_list = ["cell.id", "Tension", "time"]
#
# PEtab_tryjectory = PEtab_model.sample(obs_pars_imported)
# model_dir = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/Liver_regeneration' + '/YAP_Signaling_Liver_Regeneration_Model_reparametrized_further.xml'
#
# abc = pyabc.ABCSMC(PEtab_model, PEtab_prior, eucl_dist, population_size=2,
#                    eps=QuantileEpsilon(alpha=0.3), all_accepted=False)
#
# db_path = ("sqlite:///" +
#            os.path.join(tempfile.gettempdir(), "test.db"))
# history = abc.new(db_path, dict_data_imported)
# abc.run(max_nr_populations=2)

########################### NEW IMPLEMINTATION ################################

petab_problem_path = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/Liver_regeneration' + '/Meyer_MolSystBiol_2020.yaml'
petab_problem = petab_MS.Problem.from_yaml(petab_problem_path)

importer = PetabImporter(petab_problem)
PEtab_prior = importer.create_prior()
par_map_imported = importer.get_par_map()
obs_pars_imported = petab_problem.get_x_nominal_dict(scaled=True)
PEtab_par_scale = petab_problem.get_optimization_parameter_scales()
dict_data_imported = petab_problem.get_measurement_dict()
PEtab_model = importer.create_model()
PEtab_model.timeout = 900
PEtab_model.ignore_list = ["cell.id", "Tension", "time"]
Petab_obj_func = importer.get_objective_function()
PEtab_tryjectory = PEtab_model.sample(obs_pars_imported)
model_dir = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/Liver_regeneration' + '/YAP_Signaling_Liver_Regeneration_Model_reparametrized_further.xml'

abc = pyabc.ABCSMC(PEtab_model, PEtab_prior, Petab_obj_func, population_size=2,
                   eps=QuantileEpsilon(alpha=0.3))

db_path = ("sqlite:///" +
           os.path.join(tempfile.gettempdir(), "test.db"))
history = abc.new(db_path, dict_data_imported)
abc.run(max_nr_populations=2)

