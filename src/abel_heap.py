import os
import numpy as np
import pandas as pd
import sys
from sklearn.metrics import auc
from tqdm import tqdm
from joblib import Parallel, delayed
import subprocess

def calculate_conditional_probability(L, ETA, A, model, K_BINS, use_cache):
    inp = os.path.join("data", f"train_{L}_{model}.out")
    out = os.path.join("data", "conditional_probability", f"train_{L}_{model}_{ETA}_{A}_{K_BINS}.out")
    calc = os.path.join(".", "exec_files", "conditional_probability_L_Eta_A_Kbins")

    if not use_cache or not os.path.exists(out):
        subprocess.run(
            [calc, inp, out, str(L), str(ETA), str(A), str(K_BINS)],
            check=True
        )

    with open(out) as cond_prob:
        decision_variables = np.array(list(map(float, cond_prob.readline().split())))
        cond_prob_X1y = np.array(list(map(float, cond_prob.readline().split())))

    return decision_variables, cond_prob_X1y

def calculate_conditional_probability_parallel(args, K_BINS, use_cache, njobs):
    results = list(Parallel(n_jobs=njobs)(delayed(calculate_conditional_probability)(L, ETA, A, model, K_BINS, use_cache) 
        for L, ETA, A, model in args))
    
    return pd.DataFrame({"L": list(map(lambda x: x[0], args)),
                         "ETA": list(map(lambda x: x[1], args)),
                         "A": list(map(lambda x: x[2], args)),
                         "model": list(map(lambda x: x[3], args)),
                         "results": results})

def calculate_roc_curve(L, ETA, A, model, K_BINS, use_cache):
    inp = os.path.join("data", f"test_{L}_{model}.out")
    inp_cond_prob = os.path.join("data", "conditional_probability", f"train_{L}_{model}_{ETA}_{A}_{K_BINS}.out")
    out = os.path.join("data", "AUC", f"test_{L}_{ETA}_{A}_{K_BINS}.out")
    calc = os.path.join(".", "exec_files", "roc_curve_L_Eta_A")

    if not use_cache or not os.path.exists(out):
        subprocess.run(
            [calc, inp, inp_cond_prob, out, str(L), str(ETA), str(A)],
            check=True
        )

    with open(out) as auc:
        fpr = np.array(list(map(float, auc.readline().split())))
        tpr = np.array(list(map(float, auc.readline().split())))
        thresholds = np.array(list(map(float, auc.readline().split())))

    return fpr, tpr, thresholds

def calculate_roc_curve_parallel(args, K_BINS, use_cache, njobs):
    results = list(Parallel(n_jobs=njobs)(delayed(calculate_roc_curve)(L, ETA, A, model, K_BINS, use_cache) 
        for L, ETA, A, model in args))
    
    return pd.DataFrame({"L": list(map(lambda x: x[0], args)),
                         "ETA": list(map(lambda x: x[1], args)),
                         "A": list(map(lambda x: x[2], args)),
                         "model": list(map(lambda x: x[3], args)),
                         "results": results})

def calculate_events(L, ETAs, model):
    inp = os.path.join("data", f"test_{L}_{model}.out")
    calc = os.path.join(".", "exec_files", "event_rate")
    stdout = subprocess.run([calc, inp],
                            check=True,
                            text=True,
                            input=' '.join(list(map(str, ETAs))),
                            capture_output=True).stdout

    return list(map(int, stdout.split()))
