from qiskit.circuit import QuantumCircuit, ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qc4pyscf.operator import Second_Quantization


def single(C, D, N_orb, N_alpha, N_beta, alg_type, spin_symm):
    # [Input]           C: SparcePauliOp list of Creation Operator
    # [Input]           D: SparcePauliOp list of Annihilation Operator
    # [Input]       N_orb: Number of total basis functions
    # [Input]     N_alpha: Number of alpha spin electron
    # [Input]      N_beta: Number of beta spin electron
    # [Input]    alg_type: 0 for partial excitation (e.g. UCC) | 1 for all excitation (e.g. ADAPT-VQE)
    # [Input]   spin_symm: spin symmetry applied. If True, same type of alpha and beta spin excitation will have same amplitude

    ex_ops  = []
    ex_ind  = []
    if spin_symm:
        sum_max = N_orb if alg_type else N_alpha
        sum_min = 0 if alg_type else N_alpha
        if N_alpha != N_orb:
            for i in range(0, sum_max):
                for a in range(sum_min, N_orb): # i: occ. | a: virt.
                    T_ai    = C[a] @ D[i]
                    T_bj    = C[N_orb + a] @ D[N_orb + i]
                    T_ai_dg = T_ai.conjugate().transpose()
                    T_bj_dg = T_bj.conjugate().transpose()
                    T_pre   = T_ai + T_bj - T_ai_dg - T_bj_dg
                    T       = SparsePauliOp.simplify(T_pre)
                    ex_ops.append(T)
                    ex_ind.append(f"a{a}//a{i} (+) b{a}//b{i}")
        else:
            print('Excitation Denied')
    else:
        sum_max_a = N_orb if alg_type else N_alpha
        sum_min_a = 0 if alg_type else N_alpha
        sum_max_b = N_orb if alg_type else N_beta
        sum_min_b = 0 if alg_type else N_beta
        if N_alpha != N_orb:
            for i in range(0, sum_max_a):
                for a in range(sum_min_a, N_orb): # i: occ. | a: virt.
                    T_ai    = C[a] @ D[i]
                    T_ai_dg = T_ai.conjugate().transpose()
                    T_pre   = T_ai - T_ai_dg
                    T       = SparsePauliOp.simplify(T_pre)
                    ex_ops.append(T)
                    ex_ind.append(f"a{a}//a{i}")
        else:
            print('Alpha Excitation Denied')
        
        if N_beta != N_orb:
            for j in range(0, sum_max_b):
                sum_min = 0 if alg_type else N_beta
                for b in range(sum_min_b, N_orb): # i: occ. | a: virt.
                    T_bj    = C[N_orb + b] @ D[N_orb + j]
                    T_bj_dg = T_bj.conjugate().transpose()
                    T_pre   = T_bj - T_bj_dg
                    T       = SparsePauliOp.simplify(T_pre)
                    ex_ops.append(T)
                    ex_ind.append(f"b{b}//b{j}")
        else:
            print('Beta Excitation Denied')
    return [ex_ops, ex_ind]


def double(C, D, N_orb, N_alpha, N_beta, alg_type, spin_symm):
    # [Input]           C: SparcePauliOp list of Creation Operator
    # [Input]           D: SparcePauliOp list of Annihilation Operator
    # [Input]       N_orb: Number of total basis functions
    # [Input]     N_alpha: Number of alpha spin electron
    # [Input]      N_beta: Number of beta spin electron
    # [Input]    alg_type: 0 for partial excitation (e.g. UCC) | 1 for all excitation (e.g. ADAPT-VQE)
    # [Input]   spin_symm: spin symmetry applied. If True, same type of alpha and beta spin excitation will have same amplitude

    ex_ops  = []
    ex_ind  = []
    if alg_type:
        for i in range(0, N_orb):
            for j in range(0, N_orb):
                for a in range(0, N_orb): 
                    for b in range(0, N_orb): # ij: occ. | ab: virt.
                        T_aa        = C[a] @ C[b] @ D[j] @ D[i]
                        T_ab        = C[a] @ C[N_orb+b] @ D[N_orb+j] @ D[i]
                        T_ba        = C[N_orb+a] @ C[b] @ D[j] @ D[N_orb+i]
                        T_bb        = C[N_orb+a] @ C[N_orb+b] @ D[N_orb+j] @ D[N_orb+i]

                        if spin_symm:
                            T_11_pre    = T_aa + T_bb - T_aa.conjugate().transpose() - T_bb.conjugate().transpose()
                            T_12_pre    = T_ab + T_ba - T_ab.conjugate().transpose() - T_ba.conjugate().transpose()
                            T_11        = SparsePauliOp.simplify(T_11_pre)
                            T_12        = SparsePauliOp.simplify(T_12_pre)

                            ex_ops.append(T_11)
                            ex_ops.append(T_12)
                            ex_ind.append(f"a{a}/a{b}//a{j}/a{i} (+) b{a}/b{b}//b{j}/b{i}")
                            ex_ind.append(f"a{a}/b{b}//b{j}/a{i} (+) b{a}/a{b}//a{j}/b{i}")
                        else:
                            T_11_a_pre    = T_aa - T_aa.conjugate().transpose()
                            T_11_b_pre    = T_bb - T_bb.conjugate().transpose()
                            T_12_a_pre    = T_ab - T_ab.conjugate().transpose()
                            T_12_b_pre    = T_ba - T_ba.conjugate().transpose()
                            T_11_a        = SparsePauliOp.simplify(T_11_a_pre)
                            T_11_b        = SparsePauliOp.simplify(T_11_b_pre)
                            T_12_a        = SparsePauliOp.simplify(T_12_a_pre)
                            T_12_b        = SparsePauliOp.simplify(T_12_b_pre)

                            ex_ops.append(T_11_a)
                            ex_ops.append(T_11_b)
                            ex_ops.append(T_12_a)
                            ex_ops.append(T_12_b)
                            ex_ind.append(f"a{a}/a{b}//a{j}/a{i}")
                            ex_ind.append(f"b{a}/b{b}//b{j}/b{i}")
                            ex_ind.append(f"a{a}/b{b}//b{j}/a{i}")
                            ex_ind.append(f"b{a}/a{b}//a{j}/b{i}")
    else:
        if spin_symm:
            if N_alpha >= 1 and N_orb - N_alpha >= 1: 
                for i in range(0, N_alpha):
                    for j in range(0, N_alpha):
                        for a in range(N_alpha, N_orb): 
                            for b in range(N_alpha, N_orb): # ij: occ. | ab: virt.
                                T_aa        = C[a] @ C[b] @ D[j] @ D[i]
                                T_ab        = C[a] @ C[N_orb + b] @ D[N_orb + j] @ D[i]
                                T_ba        = C[N_orb + a] @ C[b] @ D[j] @ D[N_orb + i]
                                T_bb        = C[N_orb + a] @ C[N_orb + b] @ D[N_orb + j] @ D[N_orb + i]
                                T_11_pre    = T_aa + T_bb - T_aa.conjugate().transpose() - T_bb.conjugate().transpose()
                                T_12_pre    = T_ab + T_ba - T_ab.conjugate().transpose() - T_ba.conjugate().transpose()
                                T           = SparsePauliOp.simplify(T_11_pre + T_12_pre)

                                ex_ops.append(T)
                                ex_ind.append(f"a{a}/a{b}//a{j}/a{i} (+) b{a}/b{b}//b{j}/b{i} (+) a{a}/b{b}//b{j}/a{i} (+) b{a}/a{b}//a{j}/b{i}")
            else:
                print('Double Excitation Denied')

        else:
            # Alpha - Alpha
            if N_orb - N_alpha >= 2 and N_alpha >= 2: 
                for i in range(0, N_alpha-1):
                    for j in range(i+1, N_alpha):
                        for a in range(N_alpha, N_orb-1): 
                            for b in range(a+1, N_orb): # ij: occ. | ab: virt.
                                T_aa        = C[a] @ C[b] @ D[j] @ D[i]
                                T_11_pre    = T_aa - T_aa.conjugate().transpose()
                                T_11        = SparsePauliOp.simplify(T_11_pre)

                                ex_ops.append(T_11)
                                ex_ind.append(f"a{a}/a{b}//a{j}/a{i}")
            else:
                print('Alpha - Alpha Excitation Denied')
            # Alpha - Beta
            if N_orb != N_alpha and N_orb != N_beta and N_alpha >= 1 and N_beta >= 1:
                for i in range(0, N_alpha):
                    for j in range(0, N_beta):
                        for a in range(N_alpha, N_orb): 
                            for b in range(N_beta, N_orb): # ij: occ. | ab: virt.
                                T_ab        = C[a] @ C[N_orb + b] @ D[N_orb + j] @ D[i]
                                T_11_pre    = T_ab - T_ab.conjugate().transpose()
                                T_11        = SparsePauliOp.simplify(T_11_pre)

                                ex_ops.append(T_11)
                                ex_ind.append(f"a{a}/b{b}//b{j}/a{i}")
            else:
                print('Alpha - Beta Excitation Denied')
            # Beta - Beta
            if N_orb - N_beta >= 2 and N_beta >= 2:
                for i in range(0, N_beta-1):
                    for j in range(i+1, N_beta):
                        for a in range(N_beta, N_orb-1): 
                            for b in range(a+1, N_orb): # ij: occ. | ab: virt.
                                T_bb        = C[N_orb + a] @ C[N_orb + b] @ D[N_orb + j] @ D[N_orb + i]
                                T_11_pre    = T_bb - T_bb.conjugate().transpose()
                                T_11        = SparsePauliOp.simplify(T_11_pre)

                                ex_ops.append(T_11)
                                ex_ind.append(f"b{a}/b{b}//b{j}/b{i}")
            else:
                print('Beta - Beta Excitation Denied')
    return [ex_ops, ex_ind]


def UCC(ex_code, N_orb, N_alpha, N_beta, mapping, spin_symm):
    # Initialization
    C, D            = Second_Quantization.creators_destructors(2 * N_orb, mapping)
    ex_code         = ex_code.lower().strip()
    op_pool         = []
    op_ind          = []

    # Generate Excitation Operator Pool
    if 's' in ex_code:
        ex_ops      = single(C, D, N_orb, N_alpha, N_beta, 0, spin_symm)
        op_pool    += ex_ops[0]
        op_ind     += ex_ops[1]
    if 'd' in ex_code:
        ex_ops      = double(C, D, N_orb, N_alpha, N_beta, 0, spin_symm)
        op_pool    += ex_ops[0]
        op_ind     += ex_ops[1]

    return [op_pool, op_ind]


def ADAPT_VQE(ex_code, N_orb, N_alpha, N_beta, mapping, spin_symm):
    # Initialization
    C, D            = Second_Quantization.creators_destructors(2 * N_orb, mapping)
    ex_code         = ex_code.lower().strip()
    op_pool         = []
    op_ind          = []

    # Generate Excitation Operator Pool
    if 's' in ex_code:
        ex_ops      = single(C, D, N_orb, N_alpha, N_beta, 1, spin_symm)
        op_pool    += ex_ops[0]
        op_ind     += ex_ops[1]
    if 'd' in ex_code:
        ex_ops      = double(C, D, N_orb, N_alpha, N_beta, 1, spin_symm)
        op_pool    += ex_ops[0]
        op_ind     += ex_ops[1]

    return [op_pool, op_ind]