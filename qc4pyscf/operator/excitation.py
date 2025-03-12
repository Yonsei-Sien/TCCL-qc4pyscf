from qiskit.quantum_info import SparsePauliOp
from qc4pyscf.operator import second_q


def single(C, D, nspace_orb, nalpha, nbeta, dagger=True):
    # [Input]           C: SparcePauliOp list of Creation Operator
    # [Input]           D: SparcePauliOp list of Annihilation Operator
    # [Input]  nspace_orb: Number of total basis functions
    # [Input]      nalpha: Number of alpha spin electron
    # [Input]       nbeta: Number of beta spin electron
    # [Input]      dagger: Whether make T dagger list or not
    # [Output]       list: List of Single Excitation Operators, and its dagger

    T = []
    T_dagger = []

    if nalpha != nspace_orb:
        for i in range(0, nalpha):
            for a in range(nalpha, nspace_orb): # i: occ. | a: virt.
                Tai = C[a] @ D[i]
                T.append(Tai)
                if dagger:
                    T_dagger.append(Tai.conjugate().transpose())
    else:
        print('All alpha spin orbitals are filled.')

    if nbeta != nspace_orb:
        for j in range(0, nbeta):
            for b in range(nbeta, nspace_orb): # j: occ. | b: virt.
                Tbj = C[nspace_orb+b] @ D[nspace_orb+j]
                T.append(Tbj)
                if dagger:
                    T_dagger.append(Tbj.conjugate().transpose())
    else:
        print('All beta spin orbitals are filled.')
        
    return [T, T_dagger]


def double(C, D, nspace_orb, nalpha, nbeta, dagger=True):
    # [Input]           C: SparcePauliOp list of Creation Operator
    # [Input]           D: SparcePauliOp list of Annihilation Operator
    # [Input]  nspace_orb: Number of total basis functions
    # [Input]      nalpha: Number of alpha spin electron
    # [Input]       nbeta: Number of beta spin electron
    # [Input]      dagger: Whether make T dagger list or not
    # [Output]       list: List of Double Excitation Operators, and its dagger

    T = []
    T_dagger = []

    # alpha - alpha
    if nspace_orb - nalpha >= 2:
        for i in range(0, nalpha-1):
            for j in range(i+1, nalpha):
                for a in range(nalpha, nspace_orb-1): # ij: occ. | ab: virt.
                    for b in range(a+1, nspace_orb):
                        Taa = C[a] @ C[b] @ D[j] @ D[i]
                        T.append(Taa)
                        if dagger:
                            T_dagger.append(Taa.conjugate().transpose())
    else:
        print('Alpha - alpha excitation denied')

    # alpha - beta
    if nalpha != nspace_orb and nbeta != nspace_orb:
        for i in range(0, nalpha):
            for j in range(0, nbeta):
                for a in range(nalpha, nspace_orb): # ij: occ. | ab: virt.
                    for b in range(nbeta, nspace_orb):
                        Tab = C[a] @ C[nspace_orb + b] @ D[nspace_orb + j] @ D[i]
                        T.append(Tab)
                        if dagger:
                            T_dagger.append(Tab.conjugate().transpose())
    else:
        print('Alpha - beta excitation denied')

    # beta - beta
    if nspace_orb - nbeta >= 2:
        for i in range(0, nbeta-1):
            for j in range(i+1, nbeta):
                for a in range(nbeta, nspace_orb-1): # ij: occ. | ab: virt.
                    for b in range(a+1, nspace_orb):
                        Tbb = C[nspace_orb + a] @ C[nspace_orb + b] @ D[nspace_orb + j] @ D[nspace_orb + i]
                        T.append(Tbb)
                        if dagger:
                            T_dagger.append(Tbb.conjugate().transpose())
    else:
        print('Beta - beta excitation denied')
        
    return [T, T_dagger]


def triple(C, D, nspace_orb, nalpha, nbeta, dagger=True):
    # [Input]           C: SparcePauliOp list of Creation Operator
    # [Input]           D: SparcePauliOp list of Annihilation Operator
    # [Input]  nspace_orb: Number of total basis functions
    # [Input]      nalpha: Number of alpha spin electron
    # [Input]       nbeta: Number of beta spin electron
    # [Input]      dagger: Whether make T dagger list or not
    # [Output]       list: List of STriple Excitation Operators, and its dagger

    T = []
    T_dagger = []

    # alpha - alpha - alpha
    if nspace_orb - nalpha >= 3:
        for i in range(0, nalpha-2):
            for j in range(i+1, nalpha-1):
                for k in range(j+1, nalpha):
                    for a in range(nalpha, nspace_orb-2): # ijk: occ. | abc: virt.
                        for b in range(a+1, nspace_orb-1):
                            for c in range(b+1, nspace_orb):
                                Taaa = C[a] @ C[b] @ C[c] @ D[k] @ D[j] @ D[i]
                                T.append(Taaa)
                                if dagger:
                                    T_dagger.append(Taaa.conjugate().transpose())
    else:
        print('Alpha - alpha - alpha excitation denied')

    # alpha - alpha - beta
    if nspace_orb - nalpha >= 2 and nbeta != nspace_orb:
        for i in range(0, nalpha-1):
            for j in range(i+1, nalpha):
                for k in range(0, nbeta):
                    for a in range(nalpha, nspace_orb-1): # ijk: occ. | abc: virt.
                        for b in range(a+1, nspace_orb):
                            for c in range(nbeta, nspace_orb):
                                Taab = C[a] @ C[b] @ C[nspace_orb+c] @ D[nspace_orb+k] @ D[j] @ D[i]
                                T.append(Taab)
                                if dagger:
                                    T_dagger.append(Taab.conjugate().transpose())
    else:
        print('Alpha - alpha - beta excitation denied')

    # alpha - beta - beta
    if nspace_orb - nbeta >= 2 and nalpha != nspace_orb:
        for i in range(0, nalpha):
            for j in range(0, nbeta-1):
                for k in range(j+1, nbeta):
                    for a in range(nalpha, nspace_orb): # ijk: occ. | abc: virt.
                        for b in range(nbeta, nspace_orb-1):
                            for c in range(b+1, nspace_orb):
                                Tabb = C[a] @ C[nspace_orb+b] @ C[nspace_orb+c] @ D[nspace_orb+k] @ D[nspace_orb+j] @ D[i]
                                T.append(Tabb)
                                if dagger:
                                    T_dagger.append(Tabb.conjugate().transpose())
    else:
        print('Alpha - beta - beta excitation denied')

    # beta - beta - beta
    if nspace_orb - nbeta >= 3:
        for i in range(0, nbeta-2):
            for j in range(i+1, nbeta-1):
                for k in range(j+1, nbeta):
                    for a in range(nbeta, nspace_orb-2): # ijk: occ. | abc: virt.
                        for b in range(a+1, nspace_orb-1):
                            for c in range(b+1, nspace_orb):
                                Tbbb = C[nspace_orb+a] @ C[nspace_orb+b] @ C[nspace_orb+c] @ D[nspace_orb+k] @ D[nspace_orb+j] @ D[nspace_orb+i]
                                T.append(Tbbb)
                                if dagger:
                                    T_dagger.append(Tbbb.conjugate().transpose())
    else:
        print('Beta - beta - beta excitation denied')

    return [T, T_dagger]


def excitation(ex_code, nspace_orb, nalpha, nbeta, mapping="jordan_wigner", dagger=True):
    # [Input]     ex_code: Excitation code. String with 's', 'd', or 't'
    # [Input]  nspace_orb: Number of total basis functions
    # [Input]      nalpha: Number of alpha spin electron
    # [Input]       nbeta: Number of beta spin electron
    # [Input]     mapping: Mapping Method. Currently supports 'jordan_wigner' only
    # [Input]      dagger: Whether make T dagger list or not
    # [Output]       list: List of excitation operators. Both T and T dagger.

    ex_lst = ex_code.lower().strip()
    T = []
    T_dagger = []
    C, D = second_q.creators_destructors(2*nspace_orb, mapping)
    if 's' in ex_lst:
        singles = single(C, D, nspace_orb, nalpha, nbeta, dagger)
        T += singles[0]
        T_dagger += singles[1]

    if 'd' in ex_lst:
        doubles = double(C, D, nspace_orb, nalpha, nbeta, dagger)
        T += doubles[0]
        T_dagger += doubles[1]
        
    if 't' in ex_lst:
        triples = triple(C, D, nspace_orb, nalpha, nbeta, dagger)
        T += triples[0]
        T_dagger += triples[1]
    return [T, T_dagger]