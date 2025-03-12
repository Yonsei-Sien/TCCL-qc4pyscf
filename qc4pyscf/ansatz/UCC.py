from qiskit.circuit import QuantumCircuit, ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
from qc4pyscf.operator import excitation, hamiltonian
from qiskit.quantum_info import SparsePauliOp
from qc4pyscf.ansatz import Initial

class UCC:
    def __init__(self, mol, mf, ex_code='sd', mapping='jordan_wigner', cd_acc=1e-6):
        self.mol = mol
        self.mf = mf
        self.ex_code = ex_code
        self.mapping = mapping
        self.cd_acc = cd_acc
    

    def build(self):
        self.ansatz = self.gen_ansatz()
        self.H = self.gen_hamiltonian()
        print("Ansatz & Hamiltonian Build Done")


    def gen_ansatz(self, ex_code=None, mol=None, mapping=None):
        if mol == None:
            mol = self.mol
        if ex_code == None:
            ex_code = self.ex_code
        if mapping == None:
            mapping = self.mapping
        
        nspace_orb = len(mol.intor('int1e_ovlp'))
        nalpha = mol.nelec[0]
        nbeta = mol.nelec[1]
        init_circuit = Initial.initial(nspace_orb, nalpha, nbeta)
        exciation_ops_lst = excitation.excitation(ex_code, nspace_orb, nalpha, nbeta, mapping, True)
        exciation_ops = self.reduce_excitation(exciation_ops_lst[0], exciation_ops_lst[1])
        pv = ParameterVector("CC_Amp", len(exciation_ops))

        circuit = QuantumCircuit(2*nspace_orb)
        for i in range(len(exciation_ops)):
            evo = PauliEvolutionGate(exciation_ops[i], time=pv[i])
            circuit.append(evo, range(2*nspace_orb))

        ansatz = init_circuit.compose(circuit)
        self.ansatz = ansatz
        return ansatz
    

    def gen_hamiltonian(self, mol=None, mf=None, mapping=None, cd_acc=None):
        if mol == None:
            mol = self.mol
        if mf == None:
            mf = self.mf
        if mapping == None:
            mapping = self.mapping
        if cd_acc == None:
            cd_acc = self.cd_acc

        if str(type(mf)) == "<class 'pyscf.scf.uhf.UHF'>" or str(type(mf)) == "<class 'pyscf.scf.hf.RHF'>" or str(type(mf)) == "<class 'pyscf.scf.rohf.ROHF'>":
            H = hamiltonian.HF(mol, mf, mapping, cd_acc)
            self.H = H
            return H
        
        else:
            print('Warning: Unsupported Initial State.')

    
    def reduce_excitation(self, T, T_dagger):
        result = []

        for i in range(len(T)):
            reduced_op = SparsePauliOp.simplify(T[i]-T_dagger[i])
            reduced_op.coeffs = reduced_op.coeffs * 1.0j
            result.append(reduced_op)

        return result