import numpy as np
import matplotlib.pyplot as plt
import cirq as cq
import stim
import pymatching as pm
import stimcirq as sc
import pickle
from datetime import date

from src.noisify import *
from src.circuit_generation import generate_split_circuit
from src.cirq_glue import *
from src.qutrit_gates import *
from src.tableu_sim import *
from joblib import Parallel, delayed

#stim infra
with open("noise_model.pkl", "rb") as f:
    noise_model = pickle.load(f)

with open("coherence_times.pkl", 'rb') as f:
    coherence_times = pickle.load(f)

t1, t2e ,t2s = [],[],[]
for qb in noise_model:
    if f"qb{qb}" in coherence_times:
        noise_model[qb]['t1'] = coherence_times[f'qb{qb}']['t1']
        noise_model[qb]['t2_star'] = coherence_times[f'qb{qb}']['t2_star']
        noise_model[qb]['t2_echo'] = coherence_times[f'qb{qb}']['t2_echo']
        
        t1.append(coherence_times[f'qb{qb}']['t1'])
        t2s.append(coherence_times[f'qb{qb}']['t2_star'])
        t2e.append(coherence_times[f'qb{qb}']['t2_echo'])

noise_model['average']['t1'] = np.average(t1)
noise_model['average']['t2_echo'] = np.average(t2e)
noise_model['average']['t2_star'] = np.average(t2s)


# readout='xy'
# for coherent_error_qubit_id in [2]:
for coherent_error_qubit_id in [2, 3, 7, 11, 15, 16]:
    for readout in ['xy', 'yx']:
        kwargs = dict(readout=readout,init='z', shots=5000, noise_model=noise_model,average=False,
                   p_eff=None, rotate=False, factor=1, simulation='cirq', arb_init=False, cirq_method='stim',
                    return_mmts=True, coherent_error_qubit_id=coherent_error_qubit_id)
        
        
        numbers = np.arange(20)
        results_list = Parallel(n_jobs=-1, verbose=1)(delayed(get_expectation)(**kwargs) for number in numbers)
        
        full_results = defaultdict(list)
        for r in results_list:
            for k, v in r.items():
                if hasattr(v, '__iter__'):
                    full_results[k].extend(v)
                else:
                    full_results[k].append(v)
        for k,v in results_list[0].items():
            if hasattr(v, '__iter__'):
                continue
            else:
                if 'shots_' in k:
                    full_results[k] = np.sum(full_results[k])
                else:
                    full_results[k] = np.mean(full_results[k])
        
        today = date.today().strftime('%Y%m%d')
        with open(f"./{today}_results_{readout}_cirq_{coherent_error_qubit_id}_100k.pkl", 'wb') as f:
            pickle.dump(full_results, f)