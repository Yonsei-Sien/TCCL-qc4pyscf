from pyscf import ao2mo
import numpy as np


class transpile:
    def __init__(self, mol, mf, is_chkfile=False):
        self.mol        = mol
        self.spin       = mol.spin
        self.e0_core    = mol.energy_nuc()
        self.ao_e1_kin  = mol.intor('int1e_kin')
        self.ao_e1_nuc  = mol.intor('int1e_nuc')
        self.ao_h1e     = self.ao_e1_kin + self.ao_e1_nuc

        if is_chkfile:
            self.ao_mo_coeff    = mf['mo_coeff']
            self.ao_mo_occ      = mf['mo_occ']
        else:
            self.ao_mo_coeff    = mf.mo_coeff
            self.ao_mo_occ      = mf.mo_occ
        
        if self.spin != 0:
            h1e_a       = np.einsum('ia,jb,ij->ab', self.ao_mo_coeff[0], self.ao_mo_coeff[0], self.ao_h1e)
            h1e_b       = np.einsum('ia,jb,ij->ab', self.ao_mo_coeff[1], self.ao_mo_coeff[1], self.ao_h1e)
            
            mo_eri_a    = ao2mo.kernel(self.mol, self.ao_mo_coeff[0])
            mo_eri_b    = ao2mo.kernel(self.mol, self.ao_mo_coeff[1])

            h2e_a       = ao2mo.restore(1, mo_eri_a, self.ao_mo_coeff[0].shape[1])
            h2e_b       = ao2mo.restore(1, mo_eri_b, self.ao_mo_coeff[1].shape[1])

            self.mo_h1e = [h1e_a, h1e_b]
            self.mo_h2e = [h2e_a, h2e_b]
        else:
            h1e         = np.einsum('ia,jb,ij->ab', self.ao_mo_coeff, self.ao_mo_coeff, self.ao_h1e)
            mo_eri      = ao2mo.kernel(self.mol, self.ao_mo_coeff)
            h2e         = ao2mo.restore(1, mo_eri, self.ao_mo_coeff.shape[1])

            self.mo_h1e = h1e
            self.mo_h2e = h2e

    

    