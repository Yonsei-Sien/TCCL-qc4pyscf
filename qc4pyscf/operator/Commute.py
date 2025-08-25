from qiskit.quantum_info import SparsePauliOp


def get_commutator(A, B): #[A, B]
    com = A@B - B@A
    com = SparsePauliOp.simplify(com)
    return com


def get_commutators(A, Bs): #[A, B]
    commutators = [get_commutator(A, op) for op in Bs]
    return commutators