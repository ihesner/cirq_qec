import stim
import cirq as cq
import numpy as np
from src.qutrit_gates import *
from typing import Iterable


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


def cirq_noisify_stim(clean_circuit:cq.Circuit, noise_model, 
                 use_average=True, pipelined=True):
    """This method takes a noiseless cirq circuit and adds noise based on the \
    noise model provided. This noisification method reimplements the first 
    order stim noise model to make sure the two simulation methods are doing
    the same thing.

    Args:
        clean_circuit (cq.Circuit): Noiseless 
        noise_model (_type_): _description_
        use_average (bool, optional): Uses average value instead of qubit 
            specific values. Defaults to True.
        pipelined (bool, optional): Skip idle noise to qubits not being 
            measured during measurements. Defaults to True.

    Raises:
        NameError: _description_

    Returns:
        cirq.Circuit: Noisified circuit.
    """

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

        # Add Idle Noise on other qubits
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



def cirq_noisify_idle(clean_circuit:cq.Circuit, noise_model, 
                 use_average=True, pipelined=True,
                 timings={"CZ": 140e-9, "M": 400e-9, "else": 48e-9},
                 echo=True):
    """This method takes a noiseless cirq circuit and adds noise based on the 
    noise model provided. This noisification method adds more sophisticated 
    idle noise than just the stim based error model.

    Args:
        clean_circuit (cq.Circuit): Noiseless 
        noise_model (_type_): _description_
        use_average (bool, optional): Uses average value instead of qubit 
            specific values. Defaults to True.
        pipelined (bool, optional): Skip idle noise to qubits not being 
            measured during measurements. Defaults to True.

    Raises:
        NameError: _description_

    Returns:
        cirq.Circuit: Noisified circuit.
    """

    moments = []
    system_qubits = clean_circuit.all_qubits()

    # Likely need to initialize qutrits here.

    for moment in clean_circuit:
        pre_noise, post_noise, extra_noise = [], [], []
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

        # Add Idle Noise on other qubits
        t2 = "t2_echo" if echo else "t2_star"
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
                        # T1 decay
                        amp_damp = cq.AmplitudeDampingChannel(
                            gamma=1-np.exp(-timings['M']/noise_model[qb]['t1'])
                            ).on(qubit)
                        # T2 decay
                        phase_damp = cq.phase_flip(
                                p=1-np.exp(-timings['M']/noise_model[qb][t2])
                                ).on(qubit)
                        post_noise.append(amp_damp)
                        extra_noise.append(phase_damp)

                # Two-qubit Gate
                elif isinstance(op.gate, cq.CZPowGate):
                    # T1 decay
                    amp_damp = cq.AmplitudeDampingChannel(
                        gamma=1-np.exp(-timings['CZ']/noise_model[qb]['t1'])
                        ).on(qubit)
                    # T2 decay
                    phase_damp = cq.phase_flip(
                            p=1-np.exp(-timings['CZ']/noise_model[qb][t2])
                            ).on(qubit)
                    post_noise.append(amp_damp)
                    extra_noise.append(phase_damp)
                # Single qubit gates
                elif type(op.gate) in (cq.YPowGate, cq.XPowGate, 
                               cq.ops.pauli_gates._PauliY, 
                               cq.ops.pauli_gates._PauliX,
                               cq.ops.common_gates.HPowGate):
                    # T1 decay
                    amp_damp = cq.AmplitudeDampingChannel(
                        gamma=1-np.exp(-timings['else']/noise_model[qb]['t1'])
                        ).on(qubit)
                    # T2 decay
                    phase_damp = cq.phase_flip(
                            p=1-np.exp(-timings['else']/noise_model[qb][t2])
                            ).on(qubit)
                    post_noise.append(amp_damp)
                    extra_noise.append(phase_damp)
                else:
                    raise NameError(f"Is noise configured for {op.gate} type gates?")

        moments.append(cq.Moment(pre_noise))
        moments.append(moment)
        moments.append(cq.Moment(post_noise))
        moments.append(cq.Moment(extra_noise))

    return cq.Circuit(moments)
# End cirq_noisify


def qutrit_noisify_stim(clean_circuit:cq.Circuit, noise_model, 
                 use_average=True, pipelined=True,
                 timings={"CZ": 140e-9, "M": 400e-9, "else": 48e-9},
                 echo=True):
    """This method takes a noiseless cirq circuit and adds noise based on the 
    noise model provided. This noisification method adds more sophisticated 
    idle noise than just the stim based error model.

    Args:
        clean_circuit (cq.Circuit): Noiseless 
        noise_model (_type_): _description_
        use_average (bool, optional): Uses average value instead of qubit 
            specific values. Defaults to True.
        pipelined (bool, optional): Skip idle noise to qubits not being 
            measured during measurements. Defaults to True.

    Raises:
        NameError: _description_

    Returns:
        cirq.Circuit: Noisified circuit.
    """

    moments = []
    system_qubits = clean_circuit.all_qubits()

    # Likely need to initialize qutrits here.
    qutrits = []
    for qb in system_qubits:
        qutrits.append(cq.LineQid(qb.x, dimension=3))
    system_qutrits = frozenset(qutrits)

    for moment in clean_circuit:
        pre_noise, post_noise, extra_noise = [], [], []
        operated_qutrits = []

        # Adding Gate Noise
        for op in moment:
            # Keep track of operated qutrits
            for qubit in op.qubits:
                operated_qutrits.append(qubit)
            qb = "average" if use_average else op.qubits[0].x

            # Measurement Noise
            if isinstance(op.gate, cq.MeasurementGate):
                pre_noise.append(QutritBitFlip(
                    p=noise_model[qb]['ro']).on_each(*op.qubits))
                
            # Two Qubit Gate Noise
            elif type(op.gate) in (cq.CZPowGate, cq.CXPowGate,
                                   TwoQutritGate):
                q1, q2 = op.qubits[0].x, op.qubits[1].x
                p = noise_model['average']['two_qubit'] if use_average\
                    else noise_model[q1]['two_qubit'][q2]
                post_noise.append(
                    TwoQutritMixture(
                        cq.depolarize(p, n_qubits=2)).on(*op.qubits)
                )
                
            # Single Qubit Noise
            elif type(op.gate) in (cq.YPowGate, cq.XPowGate, 
                                cq.ops.pauli_gates._PauliY, 
                                cq.ops.pauli_gates._PauliX,
                               cq.ops.common_gates.HPowGate,
                               SingleQutritGate):
                post_noise.append(
                    SingleQutritMixture(
                        cq.depolarize(
                            noise_model[qb]['single'])).on_each(*op.qubits))

            # Virtual Z Noise
            elif type(op.gate) in (cq.ZPowGate, 
                                cq.ops.pauli_gates._PauliZ):
                pass # print("Virtual Z gate")
            else:
                print("Whats Happening Here?", op) # Hopefully nothing

        # Add Idle Noise on other qubits
        t2 = "t2_echo" if echo else "t2_star"
        for qubit in system_qutrits:
            if qubit in operated_qutrits:
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
                        # T1 decay
                        amp_damp = SingleQutritMixture(cq.AmplitudeDampingChannel(
                            gamma=1-np.exp(-timings['M']/noise_model[qb]['t1'])
                            )).on(qubit)
                        # T2 decay
                        phase_damp = SingleQutritMixture(cq.phase_flip(
                                p=1-np.exp(-timings['M']/noise_model[qb][t2])
                                )).on(qubit)
                        post_noise.append(amp_damp)
                        extra_noise.append(phase_damp)

                # Two-qubit Gate
                elif type(op.gate) in [cq.CZPowGate,
                                       TwoQutritGate]:
                    # T1 decay
                    amp_damp = SingleQutritMixture(cq.AmplitudeDampingChannel(
                        gamma=1-np.exp(-timings['CZ']/noise_model[qb]['t1'])
                        )).on(qubit)
                    # T2 decay
                    phase_damp = SingleQutritMixture(cq.phase_flip(
                            p=1-np.exp(-timings['CZ']/noise_model[qb][t2])
                            )).on(qubit)
                    post_noise.append(amp_damp)
                    extra_noise.append(phase_damp)
                # Single qubit gates
                elif type(op.gate) in (cq.YPowGate, cq.XPowGate, 
                               cq.ops.pauli_gates._PauliY, 
                               cq.ops.pauli_gates._PauliX,
                               cq.ops.common_gates.HPowGate,
                               SingleQutritGate):
                    # T1 decay
                    amp_damp = SingleQutritMixture(cq.AmplitudeDampingChannel(
                        gamma=1-np.exp(-timings['else']/noise_model[qb]['t1'])
                        )).on(qubit)
                    # T2 decay
                    phase_damp = SingleQutritMixture(cq.phase_flip(
                            p=1-np.exp(-timings['else']/noise_model[qb][t2])
                            )).on(qubit)
                    post_noise.append(amp_damp)
                    extra_noise.append(phase_damp)
                else:
                    raise NameError(f"Is noise configured for {op.gate} type gates?")

        moments.append(cq.Moment(pre_noise))
        moments.append(moment)
        moments.append(cq.Moment(post_noise))
        moments.append(cq.Moment(extra_noise))

    return cq.Circuit(moments)
# End cirq_noisify




class NoiseModelTest(cq.NoiseModel):
    """A default noise model that adds no noise."""

    def noisy_moments(self, 
                      moments: Iterable[cq.Moment], 
                      system_qubits: Sequence[cq.Qid]):
        return list(moments)
    # End noisy_moments

    def noisy_moment(self, 
                     moment: cq.Moment, 
                     system_qubits: Sequence[cq.Qid]
                     ) -> cq.Moment:
        return moment
    # End noisy_moment

    def noisy_operation(self, operation: cq.Operation) -> cq.Operation:
        return operation
    # End noisy_operation

    def _value_equality_values_(self):
        return None

    def __str__(self) -> str:
        return 'NoiseModelTest'

    def __repr__(self) -> str:
        return 'cirq.NoiseModelTest'

    def _has_unitary_(self) -> bool:
        return True

    def _has_mixture_(self) -> bool:
        return True