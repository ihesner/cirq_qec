import stim
import numpy as np
from tableu_sim import*

# Number of times to run the simulation
N = 5000

# Number of CZ repetitions (must be odd)
M = 11

num_qubits = 2

p_L = 0.0018 #leakage rate from RB experiment

p_D = p_L*0.5 #leakage decay rate

p1 = 1e-3 #single qubit depolarizing error

p2 = 1e-2 #two qubit depolarizing error


# Loop to run the simulation N times

#Prepares a bell state and measures Z1 and Z2



        

zz_values = []
zz_values_pos = []

for i in range(N):

    # Create a Tableau simulator instance
    simulator = stim.TableauSimulator()
    
    bits = np.zeros(num_qubits) #classical bits that denote leakage: 0 is not leaked, 1 is leaked
    
    simulator.h(0,1)
    
    #Repeat the CZ gate M times
    for x in range(M):
    
        #apply leakage aware CZ gate
        #simulator.cz(0, 1)
        apply_cz_gate(simulator,bits,[0,1],p_L,p1,p2)
    
        #In this example unleaked after every CZ gate too
    
        unleak = np.random.binomial(1, p_D, 2)
    
        for idx,b in enumerate(bits):
        
            if unleak[idx] == 1 and bits[idx]==1 :
                bits[idx] = 0
                
                #the qubit is reset to the one state
                simulator.reset(idx)
                simulator.x(idx)
            
    
    #print(bits)
    
    simulator.h(1)
    
    #simulator.measure_many(0,1)
    simulator.measure(0)
    simulator.measure(1)
    
    outcomes = simulator.current_measurement_record()
    obs = (outcomes[0]+outcomes[1])%2

    
    #print(simulator.current_measurement_record())
    #print(obs)
    
    #postselect on leaked qubits
    
    zz_values.append(obs)
    
    if np.sum(bits) == 0: zz_values_pos.append(obs)
    
    
zz_values = np.array(zz_values)
zz_values_pos = np.array(zz_values_pos)

print("ZZ expectation value")
print(1-2*np.mean(zz_values))

print("ZZ expectation value: leaked post-selected")
print(1-2*np.mean(zz_values_pos))
