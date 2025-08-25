from qc4pyscf.operator.Second_Quantization import creators_destructors
from qiskit.quantum_info import SparsePauliOp
import numpy as np


def make_rdm1s(cost_func, Estimator, Ansatz, parameter, N_orb, mapping='jordan_wigner'):
    C, D = creators_destructors(2 * N_orb, mapping=mapping)
    H_a_ops = []
    H_b_ops = []
    H_ind_x = []
    H_ind_y = []
    for p in range(N_orb):
        for q in range(p, N_orb):
            H_a_ops.append(SparsePauliOp.simplify((C[p] @ D[q] + C[q] @ D[p])/2))
            H_b_ops.append(SparsePauliOp.simplify((C[N_orb+p] @ D[N_orb+q] + C[N_orb+q] @ D[N_orb+p])/2))
            H_ind_x.append(p)
            H_ind_y.append(q)
    print(f"RDM-1 Initialization Done with {mapping}. #Operator: {len(H_a_ops)*2}")

    H_a = cost_func(parameter, Ansatz, H_a_ops, Estimator)
    H_b = cost_func(parameter, Ansatz, H_b_ops, Estimator)

    mat_a = np.zeros((N_orb, N_orb))
    mat_b = np.zeros((N_orb, N_orb))
    for ind in range(len(H_ind_x)):
        pos_x = H_ind_x[ind]
        pos_y = H_ind_y[ind]
        if pos_x == pos_y:
            mat_a[pos_x,pos_y] = H_a[ind]
            mat_b[pos_x,pos_y] = H_b[ind]
        else:
            mat_a[pos_x,pos_y] = H_a[ind]
            mat_a[pos_y,pos_x] = H_a[ind]
            mat_b[pos_x,pos_y] = H_b[ind]
            mat_b[pos_y,pos_x] = H_b[ind]

    return np.array([mat_a, mat_b])


def make_rdm1(cost_func, Estimator, Ansatz, parameter, N_orb, mapping='jordan_wigner'):
    C, D = creators_destructors(2 * N_orb, mapping=mapping)
    H_ops = []
    H_ind_x = []
    H_ind_y = []
    for p in range(N_orb):
        for q in range(p, N_orb):
            H_ops.append(SparsePauliOp.simplify((C[p] @ D[q] + C[q] @ D[p])/2 + (C[N_orb+p] @ D[N_orb+q] + C[N_orb+q] @ D[N_orb+p])/2))
            H_ind_x.append(p)
            H_ind_y.append(q)
    print(f"RDM-1 Initialization Done with {mapping}. #Operator: {len(H_ops)}")

    H = cost_func(parameter, Ansatz, H_ops, Estimator)

    mat = np.zeros((N_orb, N_orb))
    for ind in range(len(H_ind_x)):
        pos_x = H_ind_x[ind]
        pos_y = H_ind_y[ind]
        if pos_x == pos_y:
            mat[pos_x,pos_y] = H[ind]
        else:
            mat[pos_x,pos_y] = H[ind]
            mat[pos_y,pos_x] = H[ind]

    return np.array(mat)


def mo2ao(RDM, mo_coeff, is_split=True):
    # is_split: Whether RDM exists in alpha-beta splitted
    if is_split:
        D1   = np.einsum('ij,pi,qj->pq', RDM[0], mo_coeff[0], mo_coeff[0])
        D2   = np.einsum('ij,pi,qj->pq', RDM[1], mo_coeff[1], mo_coeff[1])
        return np.array([D1, D2])
    else:
        D0   = np.einsum('ij,pi,qj->pq', RDM, mo_coeff, mo_coeff)
        return D0