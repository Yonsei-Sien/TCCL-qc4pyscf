from qc4pyscf.operator import second_q
from pyscf import ao2mo
import numpy as np


def HF(mol, mf, mapping='jordan_wigner', cd_acc=1e-6):
    # [Input]     mol: gto.mol object
    # [Input]      mf: scf.RHF or scf.UHF object. ROHF not supported
    # [Input]  cd_acc: Accuracy of Cholesky Decomposition during generating Hamiltonian
    # [Input] mapping: Mapping Method. Currently only jordan_wigner
    # [Output]      H: Second quantized Hamiltonian based on MO

    mf_coef = mf.mo_coeff
    h1e_ao = mf.get_hcore()
    ecore = mol.energy_nuc()

    if mol.spin != 0:
        mf_coef_a = mf_coef[0]
        mf_coef_b = mf_coef[1]

        h1e_a = np.einsum('ia,jb,ij->ab',mf_coef_a,mf_coef_a,h1e_ao)
        h1e_b = np.einsum('ia,jb,ij->ab',mf_coef_b,mf_coef_b,h1e_ao)

        mo_eri_a = ao2mo.kernel(mol, mf_coef_a)
        mo_eri_b = ao2mo.kernel(mol, mf_coef_b)

        h2e_a = ao2mo.restore(1, mo_eri_a, mf_coef_a.shape[1])
        h2e_b = ao2mo.restore(1, mo_eri_b, mf_coef_b.shape[1])

        H = second_q.gen_hamiltonian(1, ecore, [h1e_a, h1e_b], [h2e_a, h2e_b], mapping, cd_acc)
        return H
    else:
        h1e = np.einsum('ia,jb,ij->ab',mf_coef,mf_coef,h1e_ao)

        mo_eri = ao2mo.kernel(mol, mf_coef)
        h2e = ao2mo.restore(1, mo_eri, mf_coef.shape[1])

        H = second_q.gen_hamiltonian(0, ecore, h1e, h2e, mapping, cd_acc)
        return H