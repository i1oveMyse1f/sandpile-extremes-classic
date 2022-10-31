import os
os.chdir("../")
import sys
sys.path.insert(1, 'src/')

import abel_heap
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
#from utils import *

def plot_global_auc_curves(auc_curve, gamma, ax=None, legend=True):
    sns.lineplot(data=auc_curve, x="T", y="eps", hue="L", style="p", markers=True, palette="tab10", markersize=10, ax=ax)
    str_gamma = "{" + f"{gamma:.3}" + "}"
    if ax is None:
        #plt.title("Quality of predictions, $\gamma=$" + str(gamma))
        plt.xlabel("T")
        plt.ylabel("$\epsilon$")
        plt.xscale("log")
        if legend:
            plt.legend(bbox_to_anchor=(1.02,1.02), loc="upper left", markerscale=2)
    else:
        #ax.set_title("Quality of predictions, $\gamma=$" + str(gamma))
        ax.set_xlabel("T")
        ax.set_ylabel("$\epsilon$")
        ax.set_xscale("log")
        if legend:
            ax.legend(bbox_to_anchor=(1.02,1.02), loc="upper left", markerscale=2)

def calculate_roc_curve_gamma(gamma, PERCENTS_ETA, LS, TS, model, njobs=-1, grid_lt=True):
    if not grid:
        assert len(LS) == len(TS)
        
    # build args
    get_eta = lambda L, percent: int(L**gamma * percent)
    get_a = lambda T: np.exp(-1. / T)
    
    args = []
    add_t = []
    add_p = []
    add_gamma = []
    
    if grid:
        for L in LS:
            for T in TS:
                A = get_a(T)
                for percent in PERCENTS_ETA:
                    ETA = get_eta(L, percent)

                    add_t.append(T)
                    add_p.append(percent)
                    add_gamma.append(gamma)
                    args.append((L, ETA, A, model))
    else:
        for L, T in zip(LS, TS):
            A = get_a(T)
            for percent in PERCENTS_ETA:
                ETA = get_eta(L, percent)

                add_t.append(T)
                add_p.append(percent)
                add_gamma.append(gamma)
                args.append((L, ETA, A, model))
    
    # run subprograms
    abel_heap.calculate_conditional_probability_parallel(args, K_BINS=200, use_cache=True, njobs=njobs)
    auc_curve = abel_heap.calculate_roc_curve_parallel(args, K_BINS=200, use_cache=True, njobs=njobs)

    # calculate some useful columns
    auc_curve["eps"] = auc_curve["results"].apply(lambda x: np.min(x[0] + (1 - x[1]))) # fpr + (1 - tpr)
    auc_curve["fpr"] = auc_curve["results"].apply(lambda x: list(x[0]))
    auc_curve["tpr"] = auc_curve["results"].apply(lambda x: list(x[1]))
    auc_curve["T"] = add_t
    auc_curve["p"] = add_p
    auc_curve["gamma"] = add_gamma
    
    return auc_curve

def plot_auc_curve(auc_curve, ax=None):
    exploded = auc_curve.explode("fpr", ignore_index=True).drop(columns=["results", "tpr"])
    exploded["tpr"] = auc_curve["tpr"].explode(ignore_index=True)
    exploded["fpr"] = exploded["fpr"].astype(float)
    exploded["tpr"] = exploded["tpr"].astype(float)
    
    sns.lineplot(data=exploded[::5], x="fpr", y="tpr", hue="L", style="p", markers=True, palette="tab10", markersize=10, legend=False, ax=ax)
    if ax is not None:
        #ax.set_title(r"ROC curves, $\gamma=$" + str(gamma))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], [0, 0.2, 0.4, 0.6, 0.8, 1])
        plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], [0, 0.2, 0.4, 0.6, 0.8, 1])
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        #plt.legend(markerscale=2)
    else:
        #plt.title(r"ROC curves, $\gamma=$" + str(gamma))
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], [0, 0.2, 0.4, 0.6, 0.8, 1])
        plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1], [0, 0.2, 0.4, 0.6, 0.8, 1])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        #plt.legend(markerscale=2)
        
def plot_local_global(gamma, t_scaling, LS, PERCENTS_ETA, TS_local, TS_global, model):
    local_aucs = []

    for i, L in enumerate(LS):
        TS_local_scaled = [int(t_scaling**i * t) for t in TS_local]
        auc_curve = calculate_roc_curve_gamma(gamma, PERCENTS_ETA, [L], TS_local_scaled, model, njobs=7)
        local_aucs.append(auc_curve)

    local_auc_curves = pd.concat(local_aucs)

    global_auc_curve = calculate_roc_curve_gamma(gamma, PERCENTS_ETA, LS, TS_global, model, njobs=7)
    best_aucs = global_auc_curve.iloc[global_auc_curve.groupby(["L", "p"])["eps"].idxmin()]
    
    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    #fig.subplots_adjust(wspace=0.3)
    #plot_auc_curve(best_aucs, ax1)
    #plot_global_auc_curves(global_auc_curve, gamma, ax2)
    plot_auc_curve(best_aucs)
    plt.show()
    plot_global_auc_curves(global_auc_curve, gamma)
    plt.show()

    ax = sns.relplot(data=local_auc_curves, x="T", y="eps", style="p", hue="L", markers=True, palette="tab10", markersize=10, col="L", 
                     kind="line", facet_kws={'sharey': True, 'sharex': False})
    plt.suptitle("T(L) near to local minimum", y=1.1)
    ax.set_axis_labels("T(L)", "$\epsilon$")
    ax._legend.set_bbox_to_anchor([1.02, 0.61])
    ax._legend.set_frame_on(True)
    plt.show()

def calculate_event_rate(gamma, PERCENTS_ETA, LS, model):    
    event_rates = {
        "event_rate": [],
        "L": [],
        "p": [],
        "model": []
    }
    
    for L in LS:
        ETAs = []
        for percent in PERCENTS_ETA:
            ETAs.append(int(L**gamma * percent))
        count_events = abel_heap.calculate_events(L, ETAs, model)
        total_events = count_events[-1]
        count_events = count_events[:-1]
        
        count_events = list(map(lambda x: x / total_events, count_events))
        event_rates["event_rate"].extend(count_events)
        event_rates["L"].extend([L] * len(PERCENTS_ETA))
        event_rates["model"].extend([model] * len(PERCENTS_ETA))
        event_rates["p"].extend(PERCENTS_ETA)
    
    return pd.DataFrame(event_rates)

from shapely.geometry import Polygon
from IPython.display import display
from copy import deepcopy

def calculate_IoU_fixed_gamma(data):
    IoU = {}
    for p, by_p in data.groupby(["p"]):
        unions, intersects = None, None
        for L, by_p_l in by_p.groupby(["L"]):
            fpr = by_p_l["fpr"].iloc[0]
            tpr = by_p_l["tpr"].iloc[0]
            polygon = Polygon([(x, y) for x, y in zip(fpr, tpr)] + [(1, 0)])
            if unions is None:
                unions = deepcopy(polygon)
                intersects = deepcopy(polygon)
            else:
                intersects = intersects.intersection(polygon)
                unions = unions.union(polygon)
        
        IoU[p] = intersects.area / unions.area
        
    return IoU

def calculate_IoU(data):
    results = defaultdict(list)
    for gamma, by_gamma in data.groupby(["gamma"]):
        IoU_gamma = calculate_IoU_fixed_gamma(by_gamma)
        results["p"] += IoU_gamma.keys()
        results["IoU"] += IoU_gamma.values()
        results["gamma"] += [gamma] * len(IoU_gamma)
    return pd.DataFrame(data=results)

def calculate_roc_curve_gamma_min_event_rate(gamma, PERCENTS_ETA, LS, TS, model, njobs=-1, min_event_rate=1e-9, grid=True):
    if not grid:
        assert len(LS) == len(TS)
        
    # build args
    get_eta = lambda L, percent: int(L**gamma * percent)
    get_a = lambda T: np.exp(-1. / T)
    
    args = []
    add_t = []
    add_p = []
    add_gamma = []
    
    event_rate = calculate_event_rate(gamma, PERCENTS_ETA, LS, "determ")
    
    if grid:
        for L in LS:
            for T in TS:
                A = get_a(T)
                for percent in PERCENTS_ETA:
                    ETA = get_eta(L, percent)

                    if (event_rate[(event_rate["p"] == percent) & (event_rate["L"] == L)]["event_rate"] > min_event_rate).all():
                        add_t.append(T)
                        add_p.append(percent)
                        add_gamma.append(gamma)
                        args.append((L, ETA, A, model))
    else:
        for L, T in zip(LS, TS):
            A = get_a(T)
            for percent in PERCENTS_ETA:
                ETA = get_eta(L, percent)

                if (event_rate[(event_rate["p"] == percent) & (event_rate["L"] == L)]["event_rate"] > min_event_rate).all():
                    add_t.append(T)
                    add_p.append(percent)
                    add_gamma.append(gamma)
                    args.append((L, ETA, A, model))
                    
    # run subprograms
    abel_heap.calculate_conditional_probability_parallel(args, K_BINS=200, use_cache=True, njobs=njobs)
    auc_curve = abel_heap.calculate_roc_curve_parallel(args, K_BINS=200, use_cache=True, njobs=njobs)

    # calculate some useful columns
    auc_curve["eps"] = auc_curve["results"].apply(lambda x: np.NaN if len(x[0]) == 0 else np.min(x[0] + (1 - x[1]))) # fpr + (1 - tpr)
    auc_curve["fpr"] = auc_curve["results"].apply(lambda x: list(x[0]))
    auc_curve["tpr"] = auc_curve["results"].apply(lambda x: list(x[1]))
    auc_curve["T"] = add_t
    auc_curve["p"] = add_p
    auc_curve["gamma"] = add_gamma
    auc_curve = pd.merge(auc_curve, event_rate, on=["L", "p", "model"], how="left")
    
    return auc_curve