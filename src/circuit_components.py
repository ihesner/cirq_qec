import stim

targets = {"data": [2, 5, 11, 3, 9, 15, 7, 13, 16],
           "x_stabs": [4, 6, 12, 14],
           "z_stabs": [8, 10, 1, 17],
           "all": [4, 6, 12, 8, 10, 1, 17, 2, 5, 11, 9, 15, 7, 13, 16, 3, 14]}

device_coords = {2: [1, 1], 5: [3, 1], 11: [5, 1], 3: [1, 3], 9: [3, 3],
                 15: [5, 3], 7: [1, 5], 13: [3, 5], 16: [5, 5],
                 6: [4, 0], 4: [2, 2], 14: [4, 4], 12: [2, 6],
                 1: [0, 2], 8: [2, 4], 10: [4, 2], 17: [6, 4]}

def stabilizer_round(stab_type:str, rotate=False):
    """Returns a stabilizer measurement sequence

    Args:
        stab_type (str): Which type of stabilizer measurement, x, z, or 
            post split weight 2 z mmts (z_split).
        rotate (bool, optional): Determines weather the dynamical decoupling
            inspired X (decomposed into Y-Z pulses) included. Defaults to False.

    Raises:
        KeyError: wrong stab type.
    """
    circuit = stim.Circuit()
    if stab_type in ["z", "x"]:
        circuit += stabilizers[stab_type][:7]
        if rotate:
            circuit += stabilizers[stab_type][7:11]
        circuit += stabilizers[stab_type][11:]
    elif stab_type == "z_split":
        circuit += stabilizers["z_split"][:5]
        if rotate:
            circuit += stabilizers[stab_type][5:9]
        circuit += stabilizers[stab_type][9:]
    else:
        raise KeyError("Stabilizer types are: 'z', 'x', 'z_split'")
    
    return circuit
# End stabilizer_round


# Lazy copy paste from provided circuit
stabilizers = {
    "x": 
stim.Circuit("""TICK
SQRT_Y_DAG 11 5 2 9 3 15 16 13 7
SQRT_Y 4 6 12 14
TICK
CZ 12 7 14 9 4 2
TICK
CZ 4 5 12 13 14 15
TICK
Y 11 5 2 9 3 15 16 13 7
TICK
Z 11 5 2 9 3 15 16 13 7
TICK
CZ 6 11 4 9 14 16
TICK
CZ 4 3 6 5 14 13
TICK
SQRT_Y_DAG 6 11 5 4 2 9 3 14 15 16 13 12 7"""),
    "z":
stim.Circuit("""TICK
SQRT_Y 1 8 10 17
TICK
CZ 8 3 10 5 17 15
TICK
CZ 10 11 1 2 8 9
TICK
Y 2 3 9 13 7 5 11 15 16
TICK
Z 2 3 1 9 13 7 5 11 15 16 17
TICK
CZ 10 15 8 13 1 3
TICK
CZ 8 7 10 9 17 16
TICK
SQRT_Y_DAG 1 8 10 17"""),
    "z_split":
stim.Circuit("""TICK
SQRT_Y 1 8 10 17
TICK
CZ 8 3 17 15 10 11 1 2
TICK
Y 2 3 7 11 15 16
TICK
Z 2 3 1 7 8 11 15 10 16 17
TICK
CZ 10 15 1 3 8 7 17 16
TICK
SQRT_Y_DAG 1 8 10 17""")}