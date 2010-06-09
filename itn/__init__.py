# from mcmc import *
from model import *
from generic_mbg import stukel_invlogit, FieldStepper
import pymc as pm
import numpy as np

nugget_labels = {'sp_sub': 'V'}
obs_labels = {'sp_sub': 'eps_p_f'}
non_cov_columns = {'pos': 'float', 'neg': 'float'}

def check_data(input):
    if np.any(input.pos+input.neg)==0:
        raise ValueError, 'Some sample sizes are zero.'
    if np.any(np.isnan(input.pos)) or np.any(np.isnan(input.neg)):
        raise ValueError, 'Some NaNs in input'
    if np.any(input.pos<0) or np.any(input.neg<0):
        raise ValueError, 'Some negative values in pos and neg'

def itn_map(sp_sub, a1, a2):
    itn = sp_sub.copy('F')
    itn = stukel_invlogit(itn,a1,a2)
    return itn
    
map_postproc = [itn_map]

def itn_val(data):
    obs = data.pos
    n = data.pos + data.neg
    def f(sp_sub, a1, a2):
        return pm.stukel_invlogit(sp_sub,a1,a2)
    return obs, n, f

validate_postproc = [itn_val]

def mcmc_init(M):
    M.use_step_method(pm.gp.GPEvaluationGibbs, M.sp_sub, M.V, M.eps_p_f, ti=M.ti)
    M.use_step_method(pm.gp.GPParentAdaptiveMetropolis, [M.amp, M.scale, M.diff_degree, M.V], scales={M.amp: .001, M.scale: .001, M.diff_degree: .001, M.V: .001})
    # M.use_step_method(pm.AdaptiveMetropolis, [M.a1, M.a2, M.V], scales={M.a1: .001, M.a2: .001, M.V: .001})

metadata_keys = ['ti','fi','ui']