import cirq as cq
import numpy as np
from typing import Sequence, Tuple

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
             [     0.,      0., 1 ,      0.,      0., 0., 0., 0., 0.],
             [g[2][0], g[2][1], 0., g[2][2], g[2][3], 0., 0., 0., 0.],
             [g[3][0], g[3][1], 0., g[3][2], g[3][3], 0., 0., 0., 0.],
             [     0.,      0., 0.,      0.,      0., 1 , 0., 0., 0.],
             [     0.,      0., 0.,      0.,      0., 0., 1 , 0., 0.],
             [     0.,      0., 0.,      0.,      0., 0., 0., 1., 0.],
             [     0.,      0., 0.,      0.,      0., 0., 0., 0., 1]], 
            dtype=complex)
        
    # End _unitary_
    
    def _circuit_diagram_info_(self, args):
        return self.qubit_gate._circuit_diagram_info_(args)
### End qutrit_gate class ###


class TwoQutritMixture(cq.Gate):
    """
    Generic Qutrit gate that takes a qubit gate on converts it to 
    act on a qutrit only in the 01 subspace.
    """
    def __init__(self, two_qb_mixture):
        # super(QutritGate, self)
        self.two_qb_mixture = two_qb_mixture
    # End __init__
    
    def _qid_shape_(self):
        return (3,3)

    def _mixture_(self):
        new_mixture = []
        for p, g in cq.mixture(self.two_qb_mixture):
            new_mixture.append((p, np.array(
            [[g[0][0], g[0][1], 0., g[0][2], g[0][3], 0., 0., 0., 0.],
             [g[1][0], g[1][1], 0., g[1][2], g[1][3], 0., 0., 0., 0.],
             [     0.,      0., 1 ,      0.,      0., 0., 0., 0., 0.],
             [g[2][0], g[2][1], 0., g[2][2], g[2][3], 0., 0., 0., 0.],
             [g[3][0], g[3][1], 0., g[3][2], g[3][3], 0., 0., 0., 0.],
             [     0.,      0., 0.,      0.,      0., 1 , 0., 0., 0.],
             [     0.,      0., 0.,      0.,      0., 0., 1 , 0., 0.],
             [     0.,      0., 0.,      0.,      0., 0., 0., 1., 0.],
             [     0.,      0., 0.,      0.,      0., 0., 0., 0., 1]], 
            dtype=complex)))
        return new_mixture
    # End _mixture_
    
    def _circuit_diagram_info_(self, args):
        return self.two_qb_mixture._circuit_diagram_info_(args)
    
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


class SingleQutritMixture(cq.Gate):
    """
    Generic Qutrit gate that takes a qubit gate on converts it to 
    act on a qutrit only in the 01 subspace.
    """
    def __init__(self, single_qb_mixture):
        # super(QutritGate, self)

        self.single_qb_mixture = single_qb_mixture
    # End __init__
    
    def _qid_shape_(self):
        return (3,)
    
    def _has_mixture_(self):
        return True

    def _mixture_(self):
        new_mixture = []
        for p, g in cq.mixture(self.single_qb_mixture):
            new_mixture.append((p,             
                    np.array([[g[0][0], g[0][1], 0.],
                              [g[1][0], g[1][1], 0.],
                              [     0.,      0., 1 ]])))
        return new_mixture
    # End _mixture_
    
    def _circuit_diagram_info_(self, args):
        return self.single_qb_mixture._circuit_diagram_info_(args)
    
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



###########################################################
###             Qutrit Noise Definitions                ###
###########################################################


class QutritI(cq.Gate):
    """
    Generic Qutrit gate that takes a qubit gate on converts it to 
    act on a qutrit only in the 01 subspace.
    """
    def __init__(self):
        super(QutritI, self)
    # End __init__
    
    def _qid_shape_(self):
        return (3,)

    def _unitary_(self):
        return np.array([[1, 0, 0],
                         [0, 1, 0],
                         [0, 0, 1]])
    # End _unitary_
    
    def _circuit_diagram_info_(self, args):
        return "I"    
### End qutrit_gate class ###


class QutritX(cq.Gate):
    """
    Generic Qutrit gate that takes a qubit gate on converts it to 
    act on a qutrit only in the 01 subspace.
    """
    def __init__(self):
        # super(QutritX, self)
        pass
    # End __init__
    
    def _qid_shape_(self):
        return (3,)

    def _unitary_(self):
        return np.array([[0, 1, 0],
                         [1, 0, 0],
                         [0, 0, 1]])
    # End _unitary_
    
    def _circuit_diagram_info_(self, args):
        return "X"    
### End qutrit_gate class ###

def qutrit_bit_flip(p):
    qtI = np.array([[1-p, 0, 0],
                    [0, 1-p, 0],
                    [0, 0, 1]])
    qtX = np.array([[0, p, 0],
                    [p, 0, 0],
                    [0, 0, 1]])
    return cq.KrausChannel([qtI, qtX], key="mmt-flip")

    # return cq.MixedUnitaryChannel([(1-p, cq.unitary(QutritI())), 
    #                                (p, cq.unitary(QutritX()))], 
    #                        key="BF")
# End qutrit_bit_flip


class QutritBitFlip(cq.Gate):
    """
    Generic Qutrit gate that takes a qubit gate on converts it to 
    act on a qutrit only in the 01 subspace.
    """
    def __init__(self, p:float):
        # super(QutritX, self)
        self._p = p
        # self.error_probabilities = None
    # End __init__

    def _mixture_(self) -> Sequence[Tuple[float, np.ndarray]]:
        I = np.array([[1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        X = np.array([[0, 1, 0],
                      [1, 0, 0],
                      [0, 0, 1]])
        return [(1 - self._p, I),
              (self._p, X)]
    # End _mixture_ 

    def _has_mixture_(self) -> bool:
        return True
    
    def _value_equality_values_(self):
        return self._p
    
    def _qid_shape_(self):
        return (3,)

    # def _unitary_(self):
    #     return np.array([[0, 1, 0],
    #                      [1, 0, 0],
    #                      [0, 0, 1]])
    # # End _unitary_

    @property
    def P(self) -> float:
        return self._p

    @property
    def n_qubits(self) -> int:
        return 1
    
    def _circuit_diagram_info_(self, args):
        return f"BF(p={self._p})"
### End qutrit_gate class ###