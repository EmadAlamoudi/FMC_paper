
from time import time
import tumor2d
from fitmulticell.sumstat import SummaryStatistics as ss
import matplotlib.pyplot as plt
from string import capwords
import os
import pyabc
from fitmulticell.model import MorpheusModel
import numpy as np
import scipy


def eucl_dist(sim, obs):
    total = 0
    for key in sim:
        if key in 'loc':
            continue
        total += scipy.stats.ks_2samp(sim[key], obs[key]).statistic
    return total


pop_size = 2
min_eps = 750
min_eps_ori = min_eps
max_nr_pop = 2

# logfilepath = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper/PEtab_problems/Code/Tumor_2d/TumorStats.txt"
problempath = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper/PEtab_problems/Tumor_2D/Tumour_Spheroid_ScenI_1e.xml"

par_map = {'k_div_max': './Global/Constant[@symbol="k_div_max"]',
           'L_init': './Global/Constant[@symbol="L_init"]',
           'q_init': './Global/Constant[@symbol="q_init"]',
           'L_div': './Global/Constant[@symbol="L_div"]',
           'ke_pro': './Global/Constant[@symbol="ke_pro"]',
           'ke_deg': './Global/Constant[@symbol="ke_deg"]',
           'e_div': './Global/Constant[@symbol="e_div"]',
           }


start_time = time()
observation_par = {"k_div_max": 4.17e-2,
                   "L_init": 1.2e1,
                   "q_init": 7.5e-1,
                   "L_div": 100,
                   "ke_pro": 5e-3,
                   "ke_deg": 8e-4,
                   "e_div": 1e-2}


sumstat = ss(output_file="logger_1.csv", ignore=["cell.id", "time"])
model = MorpheusModel(
    model_file=problempath,
    par_map=par_map,
    executable="/home/emad/morpheus-2.2.5",
    sumstat=sumstat,
)

observation_morpheus = model.sample(observation_par)
model.par_scale = "log10"
# observation_origin = tumor2d.simulate(division_rate=4.17e-2,
#                                       initial_spheroid_radius=1.2e1,
#                                       initial_quiescent_cell_fraction=7.5e-1,
#                                       division_depth=100,
#                                       ecm_production_rate=5e-3,
#                                       ecm_degradation_rate=8e-4,
#                                       ecm_division_threshold=1e-2)

limits = dict(k_div_max=(-3, -1),
              L_init=(1, 3),
              q_init=(0, 1.2),
              L_div=(-5, 0),
              ke_pro=(-5, 0),
              ke_deg=(-5, 0),
              e_div=(-5, 0))
#
prior = pyabc.Distribution(**{key: pyabc.RV("uniform", a, b - a)
                              for key, (a, b) in limits.items()})

# data_mean = tumor2d.load_default()[1]  # (raw, mean, var)

# In[6]:


# redis_sampler = pyabc.sampler.RedisEvalParallelSampler(host=host, port=port, look_ahead = False)

abc = pyabc.ABCSMC(models=model,
                   parameter_priors=prior,
                   distance_function=eucl_dist,
                   population_size=pop_size)

db_path = "sqlite:///" + "/tmp/" + "test_14param_Felipe.db"

abc.new(db_path, observation_morpheus)
history_f = abc.run(max_nr_populations=max_nr_pop, minimum_epsilon=min_eps_ori)





# petab_problem_path = "/home/emad/Insync/blackhand.3@gmail.com/Google_Drive/Bonn/Github/FMC_paper" + '/PEtab_problems' + '/Tumor_2D' + '/Tumor_2D.yaml'
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
