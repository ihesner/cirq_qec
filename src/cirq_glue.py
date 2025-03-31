import cirq as cq
import numpy as np
import pymatching as pm
import stimcirq as sc
import stim

def cirq_noisify_stim(clean_circuit:cq.Circuit, noise_model, 
                 use_average=True, pipelined=True):

    moments = []
    system_qubits = clean_circuit.all_qubits()

    # Likely need to initialize qutrits here.

    for moment in clean_circuit:
        pre_noise, post_noise = [], []
        operated_qubits = []

        # Adding Gate Noise
        for op in moment:

            for qubit in op.qubits:
                operated_qubits.append(qubit)
            qb = "average" if use_average else op.qubits[0].x

            # Measurement Noise
            if isinstance(op.gate, cq.MeasurementGate):
                pre_noise.append(
                    cq.bit_flip(p=noise_model[qb]['ro']).on_each(*op.qubits))

            # Two Qubit Gate Noise
            elif type(op.gate) in (cq.CZPowGate, cq.CXPowGate):
                q1, q2 = op.qubits[0].x, op.qubits[1].x
                p = noise_model['average']['two_qubit'] if use_average\
                    else noise_model[q1]['two_qubit'][q2]
                post_noise.append(
                    cq.depolarize(p, n_qubits=2).on(*op.qubits)
                )
                
            # Single Qubit Noise
            elif type(op.gate) in (cq.YPowGate, cq.XPowGate, 
                                cq.ops.pauli_gates._PauliY, 
                                cq.ops.pauli_gates._PauliX,
                               cq.ops.common_gates.HPowGate):
                post_noise.append(
                    cq.depolarize(noise_model[qb]['single']).on_each(*op.qubits))

            # Virtual Z Noise
            elif type(op.gate) in (cq.ZPowGate, 
                                cq.ops.pauli_gates._PauliZ):
                pass # print("Virtual Z gate")
            else:
                print("Whats Happening Here?", op) # Hopefully nothing

        # Add Idle Noise
        for qubit in system_qubits:
            if qubit in operated_qubits:
                pass # No idle noise added
            else:
                qb = 'average' if use_average else qubit.x
                # Z Gates
                if type(op.gate) in (cq.ZPowGate, 
                        cq.ops.pauli_gates._PauliZ):
                    pass # Virtual Z gate, instantaneous
                # Measurement 
                elif isinstance(op.gate, cq.MeasurementGate):
                    if pipelined:
                        pass # No idling for pipelined approach
                    else:
                        post_noise.append(cq.asymmetric_depolarize(
                                        p_x=noise_model[qb]['idle']['M'][0],
                                        p_y=noise_model[qb]['idle']['M'][1],
                                        p_z=noise_model[qb]['idle']['M'][2],
                                        ).on(qubit)
                        )
                # Two-qubit Gate
                elif isinstance(op.gate, cq.CZPowGate):
                    post_noise.append(cq.asymmetric_depolarize(
                                    p_x=noise_model[qb]['idle']['CZ'][0],
                                    p_y=noise_model[qb]['idle']['CZ'][1],
                                    p_z=noise_model[qb]['idle']['CZ'][2],
                                    ).on(qubit)
                    )
                # Single qubit gates
                elif type(op.gate) in (cq.YPowGate, cq.XPowGate, 
                               cq.ops.pauli_gates._PauliY, 
                               cq.ops.pauli_gates._PauliX,
                               cq.ops.common_gates.HPowGate):
                    post_noise.append(cq.asymmetric_depolarize(
                                    p_x=noise_model[qb]['idle']['else'][0],
                                    p_y=noise_model[qb]['idle']['else'][1],
                                    p_z=noise_model[qb]['idle']['else'][2],
                                    ).on(qubit)
                    )
                else:
                    raise NameError(f"Is noise configured for {op.gate} type gates?")

        moments.append(cq.Moment(pre_noise))
        moments.append(moment)
        moments.append(cq.Moment(post_noise))

    return cq.Circuit(moments)
# End cirq_noisify


def get_operator_values(mmts,
                        dets,
                        match:pm.Matching,
                        mtrack,
                        readout='zz',
                        ):

    shots = len(mmts)
    raw, decoded, post_selected = [], [], []
    for s in range(shots):
        # Need to define logical op from mtrack and raw qubit mmts,
        # otherwise stim defaults all mmts to zero (without this yy wont be -1)

        # Center data qubit from each rep code
        log_op_targets = mtrack.get_targets([3, 15], -1)
        # Center qubit from split mmt
        if readout in ["zz", "yy"]:
            log_op_targets += mtrack.get_targets([9], -2)
        if readout in ['xx', "yy"]:
            # Other rep code qubits
            log_op_targets += mtrack.get_targets([2, 7, 11, 16], -1)
            # X-stabilizer comparisons
            log_op_targets += mtrack.get_targets([4, 6, 12, 14], [-2, -3])
        
        log_op_value = (-1)**sum([mmts[s][t] for t in log_op_targets])

        # Raw value
        raw.append(log_op_value)

        # Post-selected value
        if sum(dets[s]) == 0:
            post_selected.append(log_op_value)

        correction = match.decode(dets[s])[0]
        decoded.append((-1)**correction * log_op_value)

    return {'raw': np.average(raw), 
            'decoded': np.average(decoded), 
            'post_selected': np.average(post_selected),
            'shots_raw': len(raw), 
            'shots_decoded': len(decoded), 
            'shots_ps': len(post_selected), 
           }
# End get_operator_values


def clean_stim_conversion(converted_circuit):
    
    clean_moments = []
    for moment in converted_circuit:
        clean_gates = []
        for gate in moment:
            if not isinstance(gate, sc.ShiftCoordsAnnotation)\
                and not isinstance(gate, sc.DetAnnotation)\
                and not isinstance(gate, sc.CumulativeObservableAnnotation):
                clean_gates.append(gate)
        clean_moments.append(cq.Moment(clean_gates))

    return cq.Circuit(clean_moments)
# End clean_stim_conversion


def convert_to_cirq(circuit:stim.Circuit, noise_model=None):

    cirq_circuit = sc.stim_circuit_to_cirq_circuit(circuit)
    clean_cirq_circuit = clean_stim_conversion(cirq_circuit)
    if noise_model is not None:
        clean_cirq_circuit = cirq_noisify_stim(clean_cirq_circuit, 
                                          noise_model=noise_model)

    return clean_cirq_circuit
# End convert_to_cirq
    