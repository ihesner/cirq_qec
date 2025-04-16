import cirq as cq
import numpy as np
import pymatching as pm
import stimcirq as sc
import stim

from src.circuit_generation import generate_split_circuit
from src.noisify import *


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


def convert_to_cirq(circuit:stim.Circuit, noise_model=None,
                    simulation_noise=None):

    cirq_circuit = sc.stim_circuit_to_cirq_circuit(circuit)
    clean_cirq_circuit = clean_stim_conversion(cirq_circuit)
    if noise_model is not None:
        if simulation_noise == 'stim':
            clean_cirq_circuit = cirq_noisify_stim(clean_cirq_circuit, 
                                          noise_model=noise_model)
        elif simulation_noise == 'idle':
            clean_cirq_circuit = cirq_noisify_idle(clean_cirq_circuit, 
                                          noise_model=noise_model)
        elif simulation_noise == "with_coherent_error":
            clean_cirq_circuit = cirq_noisify_idle(clean_cirq_circuit, 
                                          noise_model=noise_model)
            clean_cirq_circuit = add_coherent_error(clean_cirq_circuit,)
        else:
            raise AttributeError("Specify a specific cirq simulation method.")
        
    return clean_cirq_circuit
# End convert_to_cirq


def add_coherent_error(cirq_circuit, rotation=.11*np.pi):

    new_moments = []
    for moment in cirq_circuit:
        new_moments.append(moment)

        if len(moment.operations) == 7 and \
        isinstance(moment.operations[0].gate, cq.MeasurementGate):
            new_moments.append(
                cq.Moment(cq.ry(rotation).on(cq.LineQubit(1)))
                )

    return cq.Circuit(new_moments)
# End add_coherent_error
    

def get_expectation(readout='zz', init='z', shots=int(1e4), 
                    noise_model=None, average=None, p_eff=None,
                    rotate=False, factor=1,
                    simulation='stim', arb_init=False,
                    cirq_method=None):
    """_summary_

    Args:
        Just inputs for generate_split_circuit() and noisify(). See those 
        docstrings for context

    Returns:
        Dict: raw, decoded, and post-selected expectation values for given circuit.
    """
    # Defining Circuit
    stim_circuit, mtrack = generate_split_circuit(initialize=init,
                    rotate=rotate, 
                    readout=readout,
                    preselect_mmt=True,
                    init_detectors=True,
                    z_detectors='all',
                    arb_init=arb_init)
    noisy_circuit = noisify(stim_circuit, 
                            noise_model=noise_model, 
                            average=average, p_eff=p_eff, factor=factor,
                            pipelined=True, virtual_z=True)
    
    # Stim/pymatching Tools
    sampler = noisy_circuit.compile_sampler()
    converter = noisy_circuit.compile_m2d_converter()
    dem = noisy_circuit.detector_error_model()
    match = pm.Matching(dem)

    # Gathering data
    if simulation == 'stim':
        mmts = sampler.sample(shots)
        dets = converter.convert(measurements=mmts, 
                                separate_observables=True)[0]
    elif simulation == 'cirq':
        simulator = cq.Simulator()

        cirq_circuit = convert_to_cirq(stim_circuit, noise_model=noise_model,
                                       simulation_noise=cirq_method)
        results = simulator.run(cirq_circuit, repetitions=shots)
        mmts = results.data.to_numpy(dtype=np.bool_)
        dets, log_ops = converter.convert(measurements=mmts, 
                                        separate_observables=True)
        
    else:
        raise KeyError("Enter either 'stim' or 'cirq' for simulation method.")
    
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
# End get_expectation