import src.circuit_components as cc
import stim
import logging

class measurement_tracker():
    """This measurement tracker makes adding detectors much easier later.
    Just add measurements to this tracker as the circuit is generated, then 
    you can easily just call get_targets() to get the last time a qubit was 
    measured. 

    Example:
    mtrack = measurement_tracker()
    mtrack.add_mmts([1,2,3,4])
    mtrack.add_mmts([1,  3  ])
    mtrack.add_mmts([1,2,3,4])
    mtrack.get_targets(3, [-1, -3]) --> [-2, -8]
    """
    def __init__(self):
        self.mmts = []

    def add_mmts(self, mmts):
        '''
        Adds measurements to measurement history.
        '''
        self.mmts += mmts

    def get_targets(self, index, mmts_ago):
        """Returns target (for use in stim.target_rec() input) for a given 
        index qubit, and how many measurements ago.

        FIXME: This works, but is quite sloppily done. Could be optimized.

        Args:
            index (int or list): qubit index(s) to find previous measurements of
            mmts_ago (int or list): For a given index, how many measurements
                previously of that particular index are you targeting. 

        Raises:
            TypeError: _description_

        Returns:
            list: targets to pass to stim.target_rec() for a detector or observable.
        """
        if isinstance(mmts_ago, int):
            if isinstance(index, int):
                indexes = []
                for m in range(len(self.mmts)):
                    if self.mmts[m] == index:
                        indexes.append(m)
                targets = indexes[mmts_ago] - len(self.mmts)
            elif isinstance(index, list):
                # FIXME: Done quick and sloppy, can be done better
                targets = []
                for i in index:
                    indexes = []
                    for m in range(len(self.mmts)):
                        if self.mmts[m] == i:
                            indexes.append(m)
                    targets.append(indexes[mmts_ago] - len(self.mmts))
            
        elif isinstance(mmts_ago, list):
            if isinstance(index, int):
                indexes = []
                for m in range(len(self.mmts)):
                    if self.mmts[m] == index:
                        indexes.append(m)
                targets = [indexes[m_ago] - len(self.mmts) 
                           for m_ago in mmts_ago]
            elif isinstance(index, list):
                indexes, targets = [], []
                for i in index:
                    for m in range(len(self.mmts)):
                        if self.mmts[m] == i:
                            indexes.append(m)
                    targets += [indexes[m_ago] - len(self.mmts) 
                                for m_ago in mmts_ago]

        else:
            raise TypeError("What are you trying to pass here? (int or list of ints)")
        
        return targets

# End measurement_tracker class #

# Sloppy Test
def test_mtrack():
    mtrack = measurement_tracker()
    mtrack.add_mmts([1,2,3,4, 6])
    mtrack.add_mmts([1,4])
    mtrack.add_mmts([1,2,3,4, 5])

    assert mtrack.get_targets(5, -1) == -1
    assert mtrack.get_targets(6, -1) == -8
    assert mtrack.get_targets(1, [-1, -3]) == [-5, -12]
    assert mtrack.get_targets([1, 2, 3], -1) == [-5, -4, -3]
    assert mtrack.get_targets([1, 4], [-1, -2]) == [-5, -7, -2, -6]
# End test_mtrack


def generate_split_circuit(initialize="z",
                           readout="zz",
                           pauli_updates=None,
                           full_idles=3,
                           split_idles=2,
                           rotate=True,
                           init_detectors=True,
                           z_detectors="all",
                           x_detectors="bulk",
                           include_observable=True,
                           preselect_mmt=True,
                           excite=False,
                           arb_init=False):
    """Generates a stim circuit based on the repetition code split operation 
    for Ilya and Micheal's upcoming split paper. 

    Args:
        initialize (str, optional): Initialization basis. Defaults to "z".
        readout (str, optional): readout basis. Defaults to "zz".
        full_idles (int, optional): Number of stabilizer rounds before
            split operation. Defaults to 3.
        split_idles (int, optional): Number of stabilizer rounds after
            split operation. Defaults to 2.
        rotate (bool, optional): Is there an X (decomposed into Y & Z) included
            during stabilizer measurements for dynamical decoupling reasons.
            Defaults to False.
        init_detectors(bool, optional): Include detectors that compare the 
            initial pre-selecting measurement round to ideal initialization.
            Defaults to False.
        z_detectors(str, optional): Defines what z detectors to include
            None - don't include any Z type detectors
            'bulk' - only include bulk detectors, not initialization or RO
            'all' - Include bulk, as well as RO/Init detectors
            Defaults to "all".
        x_detectors(str, optional): Same as z_detectors arg but with X-detectors
            Defaults to "bulk". (Only bulk detectors implemented currently)
        include_observable(bool, optional): Include the logical observable 
            associated with the readout given. If dealing with a non-deterministic
            observables, useful to set to False. Defaults to True.
        preselect_mmt(bool, optional): Include preselection measurement? 
            Defaults to True.
        excite(bool, optional): Excites to qubit before transforming bases, so 
            if init = 'z', then prepares a |1> state (similarly |-> & |-i> for
            x and y bases). Defaults to True.
        arb_init(bool, optional): Uses arbitrary state generation to initialize
            the qubit state. Defaults to True.

    Raises:
        KeyError: _description_

    Returns:
        stim.Circuit: Returns generated stim circuit for split simulation.
    """
    # Starting with Qubit coords #
    circuit = stim.Circuit()
    mtrack = measurement_tracker()
    
    # Defining qubits coords as [x_index, y_index, t, ?, bulk(0=bulk, 1=non-bulk)]
    # Det coords, x, y, t, x/z, bulk
    for qb in cc.device_coords:
        circuit.append("QUBIT_COORDS", qb, cc.device_coords[qb] + [0, 1])

    # Initialization #
    # preselection measurement
    if preselect_mmt:
        circuit.append("M", cc.targets['all'])
        mtrack.add_mmts(cc.targets["all"])

    # Initialize basis
    if arb_init:
        circuit.append("TICK")
        circuit.append("SQRT_Y", [5, 11, 7, 13])
        if initialize == 'z':
            if excite: 
                circuit.append('X', [9])
        elif initialize == "x":
            circuit.append('SQRT_Y', [9])
            if excite: 
                circuit.append('Z', [9])
        elif initialize == 'y':
            circuit.append('SQRT_X', [9])
            if excite: 
                circuit.append('Z', [9])
        else:
            raise ValueError("Initialization need to be in 'zxy'.")
    else:
        # OLD INIT
        if initialize == "x":
            circuit.append("TICK")
            circuit.append("SQRT_Y", cc.targets['data'])
        elif initialize == 'y':
            circuit.append("TICK")
            circuit.append("SQRT_Y", [5, 11, 7, 13])
            circuit.append("SQRT_X", [9])
        else:
            pass # Already in z state
        # Exciting qubit
        if excite:
            if initialize == 'z':
                circuit.append('X', [3, 9, 15])
            elif initialize == 'x':
                circuit.append('Z', [5, 9, 13])
            elif initialize == 'y':
                circuit.append('Z', [9])

    # Initialization detectors
    if init_detectors:
        for qb in cc.device_coords:
            circuit.append("DETECTOR", 
                    stim.target_rec(mtrack.get_targets(qb, -1)),
                    cc.device_coords[qb] + [0, 1, 0])
        circuit.append("SHIFT_COORDS", [], (0, 0, 1))

    # Idling full surface code (pre-split) #
    # first X cycle
    circuit += cc.stabilizer_round('x', rotate=rotate)
    circuit.append("TICK")
    circuit.append("M", cc.targets['x_stabs'])
    mtrack.add_mmts(cc.targets['x_stabs'])
    circuit.append("SHIFT_COORDS", [], (0, 0, 1))
    # wt-2 init detectors
    if arb_init:
        circuit.append("DETECTOR", stim.target_rec(-2), [2, 6, 0, 0, 0])
        circuit.append("DETECTOR", stim.target_rec(-3), [4, 0, 0, 0, 0])\
    # X type init detectors
    elif x_detectors == 'all' and initialize == 'x':
        for stab in cc.targets['x_stabs']:
            targets = [-1,-2] if preselect_mmt else [-1]
            targets = mtrack.get_targets(stab, targets)
            circuit.append("DETECTOR", 
                            [stim.target_rec(t) for t in targets], 
                            cc.device_coords[stab] + [0, 0, 0]) # TODO: Double check last 2 coords

    # Repeated cycles 
    for i in range(full_idles):
        # Z stabilizer measurement
        circuit += cc.stabilizer_round('z', rotate=rotate)
        circuit.append("TICK")
        circuit.append("M", cc.targets['z_stabs'])
        mtrack.add_mmts(cc.targets['z_stabs'])
        # z-detectors 
        if i == 0 and z_detectors == "all" and initialize == 'z':
            if arb_init:
                circuit.append("DETECTOR", stim.target_rec(-1), [0, 2, 0, 1, 0])
                circuit.append("DETECTOR", stim.target_rec(-2), [6, 4, 0, 1, 0])
            else:
                for stab in cc.targets['z_stabs']:
                    prev_rounds = [-1, -2] if preselect_mmt else [-1]
                    targets = mtrack.get_targets(stab, prev_rounds)
                    if preselect_mmt:
                        targets += mtrack.get_targets(
                            get_adj_qubits(stab, cc.device_coords), -1)
                    circuit.append("DETECTOR", 
                                    [stim.target_rec(t) for t in targets], 
                                    cc.device_coords[stab] + [0, 1, 0])

        # First round, requires initialization readouts on data qubits
        if i == 1 and z_detectors is not None:
            pass# Ignore init detectors because of arb state gen
            for stab in cc.targets['z_stabs']:
                prev_rounds = [-1, -3] if preselect_mmt else [-1]
                targets = mtrack.get_targets(stab, prev_rounds)
                # if i == 0:
                #     targets += mtrack.get_targets(
                #         get_adj_qubits(stab, cc.device_coords), -1)
                circuit.append("DETECTOR", 
                                [stim.target_rec(t) for t in targets], 
                                cc.device_coords[stab] + [0, 1, 1])
        # Bulk
        elif i > 1 and z_detectors is not None:
            for stab in cc.targets['z_stabs']:
                targets = mtrack.get_targets(stab, [-1, -3])
                circuit.append("DETECTOR", 
                                [stim.target_rec(t) for t in targets], 
                                cc.device_coords[stab] + [0, 1, 1])
        # Warn if throwing out first round detectors due to init basis.
        elif i == 0 and z_detectors == "all" and initialize != 'z':
            pass
            # logging.warning("Throwing out first round detectors because "+\
            #               "initializing in x/y results in indeterminitic "+\
            #               "first round detetectors.")
        circuit.append("SHIFT_COORDS", [], (0, 0, 1))

        # X stabilizer measurement
        circuit += cc.stabilizer_round('x', rotate=rotate)
        circuit.append("TICK")
        circuit.append("M", cc.targets['x_stabs'])
        mtrack.add_mmts(cc.targets['x_stabs'])
        # Add Split mmts during last idle round.
        if i == full_idles - 1:
            circuit.append("M", [5, 9, 13])
            mtrack.add_mmts([5, 9, 13])
        # x-detectors
        if x_detectors is not None and i > 0:
            for stab in cc.targets['x_stabs']:
                targets = mtrack.get_targets(stab, [-1, -3])
                circuit.append("DETECTOR", 
                                [stim.target_rec(t) for t in targets], 
                                cc.device_coords[stab] + [0, 0, 1])
        elif i == 0 and x_detectors == 'all':
            for stab in cc.targets['x_stabs']:
                targets = [-1,-3] if preselect_mmt else [-1]
                targets = mtrack.get_targets(stab, targets)
                circuit.append("DETECTOR", 
                                [stim.target_rec(t) for t in targets], 
                                cc.device_coords[stab] + [0, 0, 1])
        circuit.append("SHIFT_COORDS", [], (0, 0, 1))

    # Post-split idles #
    for i in range(split_idles):
        circuit += cc.stabilizer_round('z_split', rotate=rotate)
        if i != split_idles - 1:
            circuit.append("TICK")
            circuit.append("M", cc.targets["z_stabs"])
            mtrack.add_mmts(cc.targets['z_stabs'])
            # add Z-Detectors
            if z_detectors:
                for stab in cc.targets["z_stabs"]:
                    targets = mtrack.get_targets(stab, [-1, -3])
                    # Need to include split ROs in first split idle round
                    if i == 0:
                        targets += mtrack.get_targets(
                                        get_adj_qubits(stab, 
                                        qubit_coords=cc.device_coords, 
                                        only_center=True), -1)
                    circuit.append("DETECTOR", 
                                [stim.target_rec(t) for t in targets],
                                cc.device_coords[stab] + [0, 1, 0])
            circuit.append("SHIFT_COORDS", [], (0, 0, 1))

    # Readout #
    # Prepare left rep code readout
    if readout[0] == 'z':
        pass # Already in desired basis
    elif readout[0] == "x":
        circuit.append("SQRT_Y_DAG", [2, 3, 7])
    elif readout[0] == "y":
        circuit.append("SQRT_Y_DAG", [2, 7])
        circuit.append("SQRT_X_DAG", [3])
    else:
        raise KeyError("Readout for the left code should be in ['z', 'x', 'y']")
    # Prepare right rep code readout
    if readout[1] == 'z':
        pass # Already in desired basis
    elif readout[1] == "x":
        circuit.append("SQRT_Y_DAG", [11, 15, 16])
    elif readout[1] == "y":
        circuit.append("SQRT_Y_DAG", [11, 16])
        circuit.append("SQRT_X_DAG", [15])
    else:
        raise KeyError("Readout for the right code should be in ['z', 'x', 'y']")

    # Prepare Readout OLD
    # if readout == 'zz':
    #     pass # Already in desired basis
    # elif readout == "xx":
    #     circuit.append("SQRT_Y_DAG", [2, 3, 7, 11, 15, 16])
    # elif readout == "yy":
    #     circuit.append("SQRT_Y_DAG", [2, 7, 11, 16])
    #     circuit.append("SQRT_X_DAG", [3, 15])
    # else:
    #     raise KeyError("Readout should be in ['zz', 'xx', 'yy']")

    # Final Readout
    circuit.append("TICK")
    circuit.append("M", cc.targets['all'])
    mtrack.add_mmts(cc.targets['all'])

    # Final Detectors
    if z_detectors is not None:
        # Second to last detectors (last bulk detector)
        for stab in cc.targets["z_stabs"]:
            targets = mtrack.get_targets(stab, [-1, -3])
            circuit.append("DETECTOR", 
                        [stim.target_rec(t) for t in targets],
                        cc.device_coords[stab] + [0, 1, 0])
        circuit.append("SHIFT_COORDS", [], (0, 0, 1))
    if z_detectors == "all":
        # Detectors including data readouts
        for stab in cc.targets['z_stabs']:
            if (readout[0] == 'z' and stab in [1, 8]) \
            or (readout[1] == 'z' and stab in [10,17]):
                targets = mtrack.get_targets(stab, [-1, -2])
                # include data readouts
                targets += mtrack.get_targets(
                                    get_adj_qubits(stab, cc.device_coords, 
                                                    exclude_center=True), -1)
                circuit.append("DETECTOR", 
                            [stim.target_rec(t) for t in targets],
                            cc.device_coords[stab] + [0, 1, 0])

    # Logical Observable #
    # if include_observable:
        # OLD
        # # Center data qubit from each rep code
        # log_op_targets = mtrack.get_targets([3, 15], -1)
        # # Center qubit from split mmt
        # if readout in ["zz", "yy"]:
        #     log_op_targets += mtrack.get_targets([9], -2)
        # if readout in ['xx', "yy"]:
        #     # Other rep code qubits
        #     log_op_targets += mtrack.get_targets([2, 7, 11, 16], -1)
        #     # X-stabilizer comparisons
        #     log_op_targets += mtrack.get_targets(cc.targets["x_stabs"], [-2, -3])
        # circuit.append("OBSERVABLE_INCLUDE", 
        #             [stim.target_rec(t) for t in log_op_targets], 0)
    if pauli_updates is None:
        pauli_updates = readout

    log_op_targets = {}
    # Center data qubit from each rep code
    log_op_targets['combined'] = mtrack.get_targets([3, 15], -1)
    log_op_targets['left'] = mtrack.get_targets([3], -1)
    log_op_targets['right'] = mtrack.get_targets([15], -1)
    # Individual rep code corrections: other rep code qubits
    if pauli_updates[0] in ['x', 'y']:
        log_op_targets['left'] += mtrack.get_targets([2, 7], -1)
        log_op_targets['combined'] += mtrack.get_targets([2, 7], -1)
    if pauli_updates[1] in ['x', 'y']:
        log_op_targets['right'] += mtrack.get_targets([11, 16], -1)
        log_op_targets['combined'] += mtrack.get_targets([11, 16], -1)
    # Center qubit from split mmt (add to left)
    if pauli_updates[0] in ['z', 'y']:
        log_op_targets['left'] += mtrack.get_targets([9], -2)
        log_op_targets['combined'] += mtrack.get_targets([9], -2)
    # X-stabilizer comparisons 
    if pauli_updates[0] in ['x', 'y']:
        log_op_targets['left'] += mtrack.get_targets([4, 12], [-2, -3])
        log_op_targets['combined'] += mtrack.get_targets([4, 12], [-2, -3])
    if pauli_updates[1] in ['x', 'y']:
        log_op_targets['right'] += mtrack.get_targets([6, 14], [-2, -3])
        log_op_targets['combined'] += mtrack.get_targets([6, 14], [-2, -3])
    # Actually appending observable
    if include_observable in [True, "all", "combined"]:
        circuit.append("OBSERVABLE_INCLUDE", 
                    [stim.target_rec(t) for t in log_op_targets['combined']],
                    0)
    if include_observable in ['all', 'left']:
        circuit.append("OBSERVABLE_INCLUDE", 
                    [stim.target_rec(t) for t in log_op_targets['left']],
                    1)
    if include_observable in ['all', 'right']:
        circuit.append("OBSERVABLE_INCLUDE", 
                    [stim.target_rec(t) for t in log_op_targets['right']],
                    2)

    return circuit, mtrack
# End generate_split_circuit


def get_adj_qubits(index:int, qubit_coords, 
                   exclude_center=False,
                   only_center=False):
    """Returns adjacent qubits, with option of excluding center column of 
    data qubits for split specifics.

    Args:
        index (int): _description_
        qubit_coords (dict): taken from circuit_components.py
        exclude_center (bool, optional): returns adjacent qubits excluding 
            center column of data qubits. Defaults to False.
        only_center (bool, optional): Return only adjacent qubits associated 
            with center column of data qubits. Defaults to False.

    Returns:
        _type_: _description_
    """
    index_x, index_y = qubit_coords[index]
    neighbors = []
    for qb in qubit_coords:
        x, y = qubit_coords[qb]
        if abs(index_x - x) <= 1 and abs(index_y - y) <= 1 and qb != index:
            neighbors.append(qb)

    if exclude_center:
        neighbors = [n for n in neighbors if n not in [5, 9, 13]]
    elif only_center:
        neighbors = [n for n in neighbors if n in [5, 9, 13]]

    return neighbors
# End get_adj_qubits\


def generate_split_circuit_d1(initialize="z",
                            readout="zz",
                            data_qubits=[3, 9, 15],
                            anc_qubits=[8, 10],
                            pre_split_idles=1,
                            preselect_mmt=True,
    ):
    
    assert pre_split_idles > 0, "Need at least one entangling idle."
    # Checking Qubits are connected in correct way
    low_data, mid_data, high_data = sorted(data_qubits)
    low_anc, high_anc = sorted(anc_qubits)
    assert check_connection(low_data, low_anc), \
        f'{low_data} and {low_anc} are expected to be connected.'
    assert check_connection(mid_data, low_anc), \
        f'{mid_data} and {low_anc} are expected to be connected.'
    assert check_connection(mid_data, high_anc), \
        f'{mid_data} and {high_anc} are expected to be connected.'
    assert check_connection(high_data, high_anc), \
        f'{high_data} and {high_anc} are expected to be connected.'
    used_qubits = data_qubits + anc_qubits
    
    # Starting with Qubit coords #
    circuit = stim.Circuit()
    mtrack = measurement_tracker()
    
    # Defining qubits coords as [x_index, y_index, t, bulk(0=bulk, 1=non-bulk)]
    for qb in used_qubits:
        circuit.append("QUBIT_COORDS", qb, cc.device_coords[qb])

    # Initialization #
    # preselection measurement
    if preselect_mmt:
        circuit.append("M", used_qubits)
        mtrack.add_mmts(used_qubits)
    # Initialize basis
    if initialize == "x":
        circuit.append("TICK")
        circuit.append("SQRT_Y", data_qubits)
    elif initialize == 'y':
        circuit.append("TICK")
        circuit.append("SQRT_X", [low_data, high_data])
    else:
        pass # Already in z state

    # Initialization detectors
    # if init_detectors:
    #     for qb in cc.device_coords:
    #         circuit.append("DETECTOR", 
    #                 stim.target_rec(mtrack.get_targets(qb, -1)),
    #                 cc.device_coords[qb] + [0, 1, 0])
    #     circuit.append("SHIFT_COORDS", [], (0, 0, 1))
    
    # Repeated cycles 
    for i in range(pre_split_idles):
        circuit.append("TICK")
        # Stab measurement
        circuit.append("SQRT_Y", anc_qubits)
        circuit.append("TICK")
        circuit.append("CZ", [low_data, low_anc,
                              mid_data, high_anc])
        circuit.append("TICK")
        circuit.append("CZ", [mid_data, low_anc,
                              high_data, high_anc])
        circuit.append("TICK")
        circuit.append("SQRT_Y_DAG", anc_qubits)
        circuit.append("TICK")
        circuit.append("M", anc_qubits)
        mtrack.add_mmts(anc_qubits)

        # Detectors???
        circuit.append("SHIFT_COORDS", [], (0, 0, 1))

    # Prepare readout
    if readout == 'zz':
        pass # Already in desired basis
    elif readout == "xx":
        circuit.append("SQRT_Y_DAG", [low_data, high_data])
        circuit.append("TICK")
    elif readout == "yy":
        circuit.append("SQRT_X_DAG", [low_data, high_data])
        circuit.append("TICK")
    else:
        raise KeyError("Readout should be in ['zz', 'xx', 'yy']")

    # Readout
    circuit.append("M", data_qubits)
    mtrack.add_mmts(data_qubits)
    
    # Final Detectors??

    # Logical Observable


    return circuit, mtrack
# End generate_split_circuit_d1
    

def check_connection(qb1, qb2) -> bool:
    """Checks to see if two qubits are connected. 

    Args:
        qb1 (int): qb index
        qb2 (int): qb index

    Returns:
        bool: True if two qubits share a connection on an S17 device.
    """
    x1, y1 = cc.device_coords[qb1]
    x2, y2 = cc.device_coords[qb2]

    return abs(x1 - x2) == 1 and abs(y1 - y2) == 1
# End check_connection
