import stim




def noisify(old_circuit:stim.Circuit,
            noise_model:dict=None,
            average=False,
            p_eff=None,
            factor=1,
            pipelined=True,
            virtual_z=False) -> stim.Circuit:
    """
    Takes stim circuit and adds noise given input noise model. 
    Noise model structure is:
    noise_model = {qb_id: {'single': single qubit gate depolarization noise,
                        'two_qubit': two qubit depolarization noise,
                        'idle': {gate: (x, y, z) Pauli noise model
                                for gate in ['CZ', 'M', 'else'] where else 
                                is for single qubit gates}
                        'ro': readout noise}}
    Example:
    {'Average': {'single': 0.0009030342058923018,
                'two_qubit': 0.02231052272840749,
                'idle': {'CZ': (0.00067977325, 0.00067977325,  0.005121881),
                        'M': (0.00271632172, 0.00271632172, 0.0203306602),
                        'else': (0.00016318774, 0.00016318774, 0.001231651)},
                'ro': 0.010660643794117634}}

    If average, then will use the "average" entry in noise_model instead 
    of using qubit specific noise.
    factor input will decrease the noise taken from the noise model by the 
    given factor.
    p_eff treats all noise weights as the same p value.

    Args:
        old_circuit (stim.Circuit): input stim circuit to add noise to
        noise_model (dict, optional): Device noise model. Defaults to None.
        average (bool, optional): If true, uses the average entry 
            for the noise model (average entry must be provided). 
            Defaults to False.
        p_eff (_type_, optional): p factor for flat noise model. 
            Defaults to None.
        factor (int, optional): Factor of improvement compared to noise model. 
            Defaults to 1.
        pipelined (bool, optional): Determines if idle noise added to other 
            qubits after a measurement. If True, idle noise not added as in 
            a pipelined approach. Default True.

    Returns:
        stim.Circuit: Noisified stim.Circuit
    """
    new_circuit = stim.Circuit()
    for line in old_circuit:
        gate = line.name
        targets = [t.value for t in line.targets_copy()]
        
        # Non-noisy lines of stim circuit, no noise to add
        if gate == "TICK":
            new_circuit.append("TICK")
        elif gate in ["DETECTOR", "OBSERVABLE_INCLUDE"]:
            args = line.gate_args_copy()
            new_circuit.append(gate, [stim.target_rec(t) for t in targets], args)
        elif gate in ["SHIFT_COORDS", "QUBIT_COORDS"]:
            args = line.gate_args_copy()
            new_circuit.append(gate, targets, args)

        # Everything else, should only be gates to add noise to.
        else:
            # Ro noise before measurement (apply before gate)
            if gate == "M":
                if average:
                    new_circuit.append("X_ERROR", targets, 
                                    noise_model['average']['ro'] / factor)
                elif p_eff is not None:
                    new_circuit.append("X_ERROR", targets, p_eff)
                else:
                    for t in targets:
                        new_circuit.append("X_ERROR", t, 
                                           noise_model[t]['ro'] / factor)

            # Add gate
            new_circuit.append(gate, targets)

            # Adding gate noise
            if gate == "M":
                pass # Already taken care of, nothing to add (should there be?)
                # Maybe state assignment error?
            elif gate == "CZ":
                # Gate noise
                if average:
                    new_circuit.append("DEPOLARIZE2", targets, 
                                noise_model['average']['two_qubit'] / factor)
                elif p_eff is not None:
                    new_circuit.append("DEPOLARIZE2", targets, p_eff)
                else:
                    for t in range(0, len(targets), 2):
                        new_circuit.append("DEPOLARIZE2", 
                            [targets[t], targets[t+1]], 
                            noise_model[targets[t]]['two_qubit'][targets[t+1]] / factor)
            
            # Only single qubit gates  

            elif gate == "Z":
                if virtual_z:
                    pass # Perfect instantaneous virtual Z gates
                else:
                    # Gate noise
                    if average:
                        new_circuit.append("DEPOLARIZE1", targets, 
                                        noise_model['average']['single'] / factor)
                    elif p_eff is not None:
                        new_circuit.append("DEPOLARIZE1", targets, p_eff)
                    else:
                        for t in targets:
                            new_circuit.append("DEPOLARIZE1", t, 
                                            noise_model[t]['single'] / factor)
            else:
                # Gate noise
                if average:
                    new_circuit.append("DEPOLARIZE1", targets, 
                                    noise_model['average']['single'] / factor)
                elif p_eff is not None:
                    new_circuit.append("DEPOLARIZE1", targets, p_eff)
                else:
                    for t in targets:
                        new_circuit.append("DEPOLARIZE1", t, 
                                           noise_model[t]['single'] / factor)

            # Idle noise for other qubits
            idle_qubits = [qb for qb in range(1, 18) 
                        if qb not in targets + ["average"]]
            idle_type = gate if gate in ['M', 'CZ'] else 'else'
            if (len(idle_qubits) > 0 or p_eff is not None) \
            and (gate != 'Z' or not virtual_z):
                if average:
                    new_circuit.append("PAULI_CHANNEL_1", idle_qubits,
                        [noise / factor
                        for noise in noise_model['average']['idle'][idle_type]])
                # For flat p just depolarize once at end of cycle
                elif p_eff is not None:
                    if gate == "M" and not pipelined:
                        new_circuit.append("DEPOLARIZE1", range(1, 18), p_eff)
                        #FIXME: Lazy hardcoded qbs 1-17
                else:
                    for qb in idle_qubits:
                        if gate != "M" or not pipelined:
                            new_circuit.append("PAULI_CHANNEL_1", qb,
                                [noise / factor
                                for noise in noise_model[qb]['idle'][idle_type]])
                        
    return new_circuit
# End noisify
