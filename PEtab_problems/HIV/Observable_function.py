import os
from pyabc.external.r import R


def sumstat(loc):
    cur_path = os.getcwd()
    r = R(os.path.join(
        os.path.join(os.path.dirname(os.path.dirname(cur_path)), 'HIV'),
        'abc_r_blueprint.R'))
    sum_stat = r.summary_statistics("mySummaryStatistics", is_py_model=True) # If pyabc ver 0.9.22 add is_py_model=True
    print("loc: ", loc)
    return sum_stat({"loc": loc})

