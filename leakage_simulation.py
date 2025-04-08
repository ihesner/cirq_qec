import stim
import numpy as np
import os
import sys
import pickle
cwd = os.getcwd()
sys.path.append(cwd+"/src/")

from tableu_sim import*
from circuit_generation import*
from noisify import*



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


# Number of times to run the simulation
#N = 10000
N = 57456 #experiment


num_qubits = 17

p_L = 0.0018 #leakage rate from RB experiment

T_D = 1e3 * 23.24 #ns #Estimated leakage decay time


# Loop to run the simulation N times


split_circuit, meas_tracker = generate_split_circuit()

#print(split_circuit)


noisy_circuit = noisify(split_circuit, noise_model=noise_model,
                        pipelined=True, virtual_z=True,
                        average=False)


#print(noisy_circuit)

leaked_tot = 0

non_leaked_tot = 0

for x in range(N):

    # Convert the circuit to a TableauSimulator
    #detectors, observables, leaked_bits = circuit_to_tableau_simulator(circuit,p_L,T_D,num_qubits)
    
    detectors, observables, leaked_bits = circuit_to_tableau_simulator(noisy_circuit,p_L,T_D,num_qubits)

    #print("Detection events ")
    #print(detectors)
    
    #print("Observables ")
    #print(observables)
    
    #print("Leaked qubits ? ")
    #print(np.sum(leaked_bits)%2)
    
    leaked = np.sum(leaked_bits)
    
    if leaked == 0: non_leaked_tot += 1
    else : leaked_tot += 1
    
   
    
    #print("")
    
print("Leakage discard rate ", leaked_tot/N)

print("Leakage non-discard rate ", non_leaked_tot/N)

print("Experimental leakage non-discard rate = ", 0.768)
    
    # Perform measurements
