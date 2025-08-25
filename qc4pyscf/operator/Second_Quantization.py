from qiskit.quantum_info import SparsePauliOp
import numpy as np


def identity(n):
    return SparsePauliOp.from_list([("I" * n, 1)])


def creators_destructors(n, mapping="jordan_wigner"): # Where only "Mapping" works
    c_list = []
    if mapping == "jordan_wigner":
        for p in range(n):
            if p == 0:
                l, r = "I" * (n - 1), ""
            elif p == n - 1:
                l, r = "", "Z" * (n - 1)
            else:
                l, r = "I" * (n - p - 1), "Z" * p
            cp = SparsePauliOp.from_list([(l + "X" + r, 0.5), (l + "Y" + r, -0.5j)])
            c_list.append(cp)
    else:
        raise ValueError("Unsupported mapping.")
    d_list = [cp.conjugate().transpose() for cp in c_list]
    return c_list, d_list


def cholesky(V, eps):
    no = V.shape[0]
    chmax, ng = 20 * no, 0
    W = V.reshape(no**2, no**2)
    L = np.zeros((no**2, chmax))
    Dmax = np.diagonal(W).copy()
    nu_max = np.argmax(Dmax)
    vmax = Dmax[nu_max]
    while vmax > eps:
        L[:, ng] = W[:, nu_max]
        if ng > 0:
            L[:, ng] -= np.dot(L[:, 0:ng], (L.T)[0:ng, nu_max])
        L[:, ng] /= np.sqrt(vmax)
        Dmax[: no**2] -= L[: no**2, ng] ** 2
        ng += 1
        nu_max = np.argmax(Dmax)
        vmax = Dmax[nu_max]
    L = L[:, :ng].reshape((no, no, ng))
    print(
        "accuracy of Cholesky decomposition ",
        np.abs(np.einsum("prg,qsg->prqs", L, L) - V).max(),
    )
    return L, ng


def gen_H(qc_mf, mapper='jordan_wigner', cd_acc=1e-6):
    restric_checker = True
    spin    = qc_mf.spin
    ecore   = qc_mf.e0_core
    h1e     = qc_mf.mo_h1e
    h2e     = qc_mf.mo_h2e

    if spin != 0:
        restric_checker = False
        print("System is not restricted. Accuracy of hamiltonian can be lower. Please check your 'cd_acc' input.")

    if restric_checker:
        N_orb, _ = h1e.shape
                
        # 2. Make creation & annihilation operators with Pauli operators
        C, D = creators_destructors(2 * N_orb, mapping=mapper)
        Exc = []
        for p in range(N_orb):
            Excp = [C[p] @ D[p] + C[N_orb + p] @ D[N_orb + p]]
            for r in range(p + 1, N_orb):
                Excp.append(
                    C[p] @ D[r]
                    + C[N_orb + p] @ D[N_orb + r]
                    + C[r] @ D[p]
                    + C[N_orb + r] @ D[N_orb + p]
                )
            Exc.append(Excp)

        # 3. Do low-rank decomposition of the 2e Hamiltonian for operator optimization
        Lop, ng = cholesky(h2e, cd_acc)
        t1e = h1e - 0.5 * np.einsum("pxxr->pr", h2e) # decomposed 2e term's S + original 1e term

        # 4. Sum all operators with coefficients we got above
        H = ecore * identity(2 * N_orb)

        # one-body term
        for p in range(N_orb):
            for r in range(p, N_orb):
                H += t1e[p, r] * Exc[p][r - p] # r-p follows the order of the operator list

        # two-body term
        for g in range(ng):
            Lg = 0 * identity(2 * N_orb)
            for p in range(N_orb):
                for r in range(p, N_orb):
                    Lg += Lop[p, r, g] * Exc[p][r - p]
            H += 0.5 * Lg @ Lg

        return H.chop().simplify()
        
    else:
        h1e_a = h1e[0]
        h1e_b = h1e[1]
        h2e_a = h2e[0]
        h2e_b = h2e[1]
        N_orb, _ = h1e_a.shape

        # 2. Make creation & annihilation operators with Pauli operators
        C, D = creators_destructors(2 * N_orb, mapping=mapper)
        Exc_a = []
        Exc_b = []
        for p in range(N_orb):
            Excp_a = [C[p] @ D[p]]
            Excp_b = [C[N_orb+p] @ D[N_orb+p]]
            for r in range(p + 1, N_orb):
                Excp_a.append(
                    C[p] @ D[r]
                    + C[r] @ D[p]
                )
                Excp_b.append(
                    C[N_orb + p] @ D[N_orb + r]
                    + C[N_orb + r] @ D[N_orb + p]
                )
            Exc_a.append(Excp_a)
            Exc_b.append(Excp_b)

        # 3. Do low-rank decomposition of the 2e Hamiltonian for operator optimization
        Lop_a, ng_a = cholesky(h2e_a, cd_acc)
        t1e_a = h1e_a - 0.5 * np.einsum("pxxr->pr", h2e_a)
        Lop_b, ng_b = cholesky(h2e_b, cd_acc)
        t1e_b = h1e_b - 0.5 * np.einsum("pxxr->pr", h2e_b)

        # 4. Sum all operators with coefficients we got above
        H = ecore * identity(2 * N_orb)
        # one-body term
        for p in range(N_orb):
            for r in range(p, N_orb):
                H += t1e_a[p, r] * Exc_a[p][r - p]
                H += t1e_b[p, r] * Exc_b[p][r - p]
        # two-body term ~ Someday it will broke...
        for g in range(ng_a):
            Lg = 0 * identity(2 * N_orb)
            for p in range(N_orb):
                for r in range(p, N_orb):
                    Lg += Lop_a[p, r, g] * Exc_a[p][r - p]
                    Lg += Lop_b[p, r, g] * Exc_b[p][r - p]
            H += 0.5 * Lg @ Lg

        return H.chop().simplify()