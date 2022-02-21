
from time import time
from fitmulticell.sumstat import SummaryStatistics as ss
import matplotlib.pyplot as plt
from fitmulticell.PEtab.base import PetabImporter
import pyabc
from fitmulticell.model import MorpheusModel
import numpy as np
import petab_MS
import tempfile
import os


pop_size = 2
min_eps = 750
min_eps_ori = min_eps
max_nr_pop = 2

logfilepath = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper/PEtab_problems/Code/Tumor_2d/TumorStats.txt"
problempath = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper/PEtab_problems/Tumor_2D/Tumour_Spheroid_ScenI_1e.xml"

observation_par = {"k_div_max": 4.17e-2,
                   "L_init": 1.2e1,
                   "q_init": 7.5e-1,
                   "L_div": 100,
                   "ke_pro": 5e-3,
                   "ke_deg": 8e-4,
                   "e_div": 1e-2}


petab_problem_path = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/Tumor_2D' + '/Tumor_2D.yaml'
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

PEtab_obj_function = petab_problem.get_objective_function()
observation_morpheus = PEtab_model.sample(observation_par)
PEtab_model.par_scale = "log10"

abc = pyabc.ABCSMC(models=PEtab_model,
                   parameter_priors=PEtab_prior,
                   distance_function=PEtab_obj_function,
                   population_size=pop_size)

db_path = "sqlite:///" + "/tmp/" + "test_14param_Felipe.db"

abc.new(db_path, observation_morpheus)
history_f = abc.run(max_nr_populations=max_nr_pop, minimum_epsilon=min_eps_ori)

