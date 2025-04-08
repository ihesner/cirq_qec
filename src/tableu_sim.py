import numpy as np


all_gates = ["C_XYZ", "C_ZYX", "H", "H_YZ", "I","Y","SQRT_Y","SQRT_Y_DAG", "CX", "CY","CZ","XCX", "XCY", "XCZ", "YCX", "YCY", "YCZ", "R", "RX", "RY","M", "MX", "MY"]

single_qubit_gates = ["H","I","Y","SQRT_Y","SQRT_Y_DAG", "R", "RX", "RY","M", "MX", "MY"]

all_errors = ["DEPOLARIZE1", "DEPOLARIZE2", "X_ERROR", "Z_ERROR","PAULI_CHANNEL_1"]

ANNOTATION_OPS = {"OBSERVABLE_INCLUDE", "DETECTOR", "SHIFT_COORDS", "QUBIT_COORDS"}

import stim

def apply_cz_gate(simulator,bits,qubit_indices,p_L):

    """
    Apply a CZ gate including incoherent leakage 
    bits : list of classical bits (all of them) that say when a qubit has leaked
    qubit_indices : indices of control(data) and target(ancilla)
    p_L : leakge probability
    p1 : probability of single qubit error when the partner is leaked
    p2 : gate error probability
    """

    
    #depolarizing probabilities due to leakage: 0 for now
    p1 = 0
    p2 = 0
    #FIXME: check a value for depolarizing channel
    
    
    qubit1 = qubit_indices[0]
    qubit2 = qubit_indices[1]
    
    bit1 = bits[qubit1]
    bit2 = bits[qubit2]

    if bit1+bit2 == 0 : #no qubit is leaked
    
        simulator.cz(qubit1,qubit2)
        #simulator.depolarize2(qubit1,qubit2,p=p2)
        
    elif bit1 == 0: #first qubit leaked, depolarizing error on second qubit
    
        simulator.depolarize1(qubit2,p=p1) #extra error due to leakage
        
    elif bit2 == 0: #second qubit leaked, depolarizing error on first qubit
    
        simulator.depolarize1(qubit1,p=p1) #extra error due to leakage
        
    
    #update classical bits
    
    leaked = np.random.binomial(1, p_L, 2) #leak a qubit with probability p_L
    
    for idx,b in enumerate(leaked):
        
        if b == 1: bits[qubit_indices[idx]] = 1
        else: continue
        

def apply_gates(name, targets,simulator,leak_info):

    pl = leak_info[0]
    bits0 = leak_info[1]
    
    if name == 'H':
        simulator.h(*targets)
    elif name == 'S':
        simulator.s(*targets)
    elif name == 'X':
        simulator.x(*targets)
    elif name == 'Y':
        simulator.y(*targets)
    elif name == 'Z':
        simulator.z(*targets)
    elif name == 'SQRT_Y':
        simulator.sqrt_y(*targets)
    elif name == 'SQRT_Y_DAG':
        simulator.sqrt_y_dag(*targets)
    elif name == 'CNOT':
        simulator.cnot(*targets)
    elif name == 'CX':
        simulator.cnot(*targets)
    elif name == 'CZ':
        
        #simulator.cz(*targets) #standard gate
        
        qubit_indices = []
        for target in targets:
            qubit_indices.append(target.value)
            
            if len(qubit_indices) ==  2:
                apply_cz_gate(simulator,bits0,qubit_indices,pl)
                qubit_indices = []
        
    elif name == 'SWAP':
        simulator.swap(*targets)
    elif name == 'M':
        simulator.measure_many(*targets)
    elif name == 'TICK':
        pass
    else:
        raise ValueError(f"Unsupported gate: {gate}")
        
        
def apply_error(name,probabilities,targets,simulator):

    if name == 'DEPOLARIZE1':
        simulator.depolarize1(*targets,p=probabilities[0])
    elif name == 'DEPOLARIZE2':
        simulator.depolarize2(*targets,p=probabilities[0])
    elif name == 'X_ERROR':
        simulator.x_error(*targets,p=probabilities[0])
    elif name == 'Z_ERROR':
        simulator.z_error(*targets,p=probabilities[0])
    elif name == 'PAULI_CHANNEL_1':
        the_circuit = stim.Circuit()
        the_circuit.append("PAULI_CHANNEL_1", *targets,probabilities)
        simulator.do_circuit(the_circuit)
    elif name == 'PAULI_CHANNEL_2':
        the_circuit = stim.Circuit()
        the_circuit.append("PAULI_CHANNEL_2", *targets,probabilities)
        simulator.do_circuit(the_circuit)
    else:
        raise ValueError("Unsupported gate: {}".format(gate))
        
def apply_annotation(name,targets,simulator):

    measurements = simulator.current_measurement_record()
    
    if name == 'DETECTOR':
        detector = 0
        for x in targets:
            detector += measurements[x.value]
    
        return detector%2
        
    elif name == 'OBSERVABLE_INCLUDE':
        obs = 0
        for x in targets:
            obs +=  measurements[x.value]
    
        return obs%2
        
    elif name == 'QUBIT_COORDS':
        pass
        
    elif name == 'SHIFT_COORDS':
        pass
            
    else:
        print(name)
        raise ValueError(f"Unsupported gate: {gate}")


def circuit_to_tableau_simulator(circuit: stim.Circuit, p_L, T_D, num_qubits):

    """Run a noisy stim circuit as a tableu simulator
       return list of detection events and observables 
    
    """
    
    simulator = stim.TableauSimulator()
    
    detectors = []
    observables = []
    
    bits = np.zeros(num_qubits+1) #classical bits for leak/unleak state one extra bit because of indexing
    
    leak_info = [p_L,bits]
    
    for operation in circuit:
        gate = operation.name  # Access the gate type
        targets = operation.targets_copy()  # Access the target qubits
        arguments = operation.gate_args_copy() # Access arguments of the operation
        
    
        
        if gate in all_gates: #gate
            apply_gates(gate,targets,simulator,leak_info)
            
            previous_gate = gate #store the last gate/measurement operation
            
            no_idling_qubits = []
            for tg in targets:
                no_idling_qubits.append(tg.value) #store the qubits that were touche
            
        elif gate in all_errors: #error
            apply_error(gate,arguments,targets,simulator)
            #FIXME: make leakage-aware errors
            
        elif gate in ANNOTATION_OPS:
            detect_event = apply_annotation(gate,targets,simulator)
            if gate == 'DETECTOR': detectors.append(detect_event)
            if gate == 'OBSERVABLE_INCLUDE': observables.append(detect_event)
            
            #FIXME: coordinates of the stabilizers are ignored
            
        elif gate == 'TICK': #Apply leakage decay
            
            if previous_gate == 'CZ': gate_time = 100 #ns
            elif previous_gate == 'M': gate_time = 400 #ns
            elif previous_gate in single_qubit_gates : gate_time = 24 #ns
    
            #T_D = 1e3 * 23.24 #ns #Estimated leakage decay time
            p_D = 1 - np.exp(-gate_time/T_D)
            
            num_no_idling_qubits = len(no_idling_qubits)
            
            num_idling = num_qubits - num_no_idling_qubits
            
            unleak = np.random.binomial(1, p_D, num_idling)
            
            #print(unleak)
            
            indx_unleak = 0
            
            #print(no_idling_qubits)
            #print(unleak)
    
            #for idx,b in enumerate(bits):
            for idx in range(1,num_qubits+1): #Circuit qubits start from 1 and end in N
                
                if idx in no_idling_qubits:
                    pass
                
                else:
        
                    if unleak[indx_unleak] == 1 and bits[idx]==1 :
                        bits[idx] = 0
                
                        #the qubit is reset to the one state
                        simulator.reset(idx)
                        simulator.x(idx)
                    
                    indx_unleak +=1
            
        
    return detectors, observables, bits
