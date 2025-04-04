import cirq as cq
import numpy as np


class LeakyCZ(cq.Gate):
    """
    """
    def __init__(self, phi_L=0, phi_C=0):
        super(LeakyCZ, self)
        self.phi_L = phi_L
        self.phi_C = phi_C
    # End __init__
    
    def _qid_shape_(self):
        return (3, 3)

    def _unitary_(self):
        pL = self.phi_L
        pC = self.phi_C
        pi = np.pi
        H = np.array([[0 , 0., 0., 0., 0., 0., 0., 0., 0.],
                      [0., 0 , 0., 0., 0., 0., 0., 0., 0.],
                      [0., 0., 0 , 0., 0., 0., 0., 0., 0.],
                      [0., 0., 0., 0 , 0., 0., 0., 0., 0.],
                      [0., 0., 0., 0., pi, 0., pL, 0., 0.],
                      [0., 0., 0., 0., 0., 0 , 0., 0., 0.],
                      [0., 0., 0., 0., pL, 0., 0 , 0., 0.],
                      [0., 0., 0., 0., 0., 0., 0., pC, 0.],
                      [0., 0., 0., 0., 0., 0., 0., 0., 0]], dtype=complex)
        # TRASH TODO FIXME DELETE
        # Problem is this creates ones in place of all zeros when doing the 
        # np.exp part.
        # H = np.zeros((9,9), dtype=complex)
        # H[self.i][self.i] = pi
        # print(np.round(np.exp(1j * H), decimals=10))
        # return np.round(np.exp(1j * H), decimals=10)
        toreturn = np.eye(9, dtype=complex)
        toreturn[4, 4] = -1
        return toreturn
    # End _unitary
    
    def _circuit_diagram_info_(self, args):
        return '&', "&"
    
### End LeakyCZ class ###

class TwoQutritGate(cq.Gate):
    """
    Generic Qutrit gate that takes a qubit gate on converts it to 
    act on a qutrit only in the 01 subspace.
    """
    def __init__(self, qubit_gate):
        # super(QutritGate, self)
        self.qubit_gate = qubit_gate
    # End __init__
    
    def _qid_shape_(self):
        return (3,3)

    def _unitary_(self):
        g = cq.unitary(self.qubit_gate)
        return np.array(
            [[g[0][0], g[0][1], 0., g[0][2], g[0][3], 0., 0., 0., 0.],
             [g[1][0], g[1][1], 0., g[1][2], g[1][3], 0., 0., 0., 0.],
             [     0.,      0., 0 ,      0.,      0., 0., 0., 0., 0.],
             [g[2][0], g[2][1], 0., g[2][2], g[2][3], 0., 0., 0., 0.],
             [g[3][0], g[3][1], 0., g[3][2], g[3][3], 0., 0., 0., 0.],
             [     0.,      0., 0.,      0.,      0., 0 , 0., 0., 0.],
             [     0.,      0., 0.,      0.,      0., 0., 0 , 0., 0.],
             [     0.,      0., 0.,      0.,      0., 0., 0., 0., 0.],
             [     0.,      0., 0.,      0.,      0., 0., 0., 0., 0]], 
            dtype=complex)
        
    # End _unitary_
    
    def _circuit_diagram_info_(self, args):
        return self.qubit_gate._circuit_diagram_info_(args)
    
### End qutrit_gate class ###



class SingleQutritGate(cq.Gate):
    """
    Generic Qutrit gate that takes a qubit gate on converts it to 
    act on a qutrit only in the 01 subspace.
    """
    def __init__(self, qubit_gate):
        # super(QutritGate, self)
        self.qubit_gate = qubit_gate
    # End __init__
    
    def _qid_shape_(self):
        return (3,)

    def _unitary_(self):
        g = cq.unitary(self.qubit_gate)
        return np.array([[g[0][0], g[0][1], 0.],
                         [g[1][0], g[1][1], 0.],
                         [     0.,      0., 0 ]])
    # End _unitary_
    
    def _circuit_diagram_info_(self, args):
        return self.qubit_gate._circuit_diagram_info_(args)
    
### End qutrit_gate class ###



def qutritify(qubit_circuit:cq.Circuit) -> cq.Circuit:
    qutrit_moments = []
    for moment in qubit_circuit:
        qutrit_ops = []
        for op in moment:
            target_xs = [qb.x for qb in op.qubits]

            targets = [cq.LineQid(x, dimension=3) for x in target_xs]
            # print(op, targets)
            # if isinstance(op.gate, cq.CZPowGate):
            #     qutrit_ops.append(LeakyCZ(phi_L=phi_L, phi_C=phi_C
            #                               ).on_each(targets))
            if isinstance(op.gate, cq.MeasurementGate):
                key = op.gate.key
                qutrit_ops.append(cq.measure(targets, key=key))
            elif len(targets) == 2:
                qutrit_ops.append(TwoQutritGate(op.gate).on(*targets))
            elif len(targets) == 1: 
                qutrit_ops.append(SingleQutritGate(op.gate).on(targets[0]))
            else:
                raise TypeError(f"Not convertion compatable - targets:{targets}")

        qutrit_moments.append(cq.Moment(qutrit_ops))

    return cq.Circuit(qutrit_moments)
# End qutritify
