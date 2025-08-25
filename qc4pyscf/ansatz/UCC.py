from qiskit.circuit import QuantumCircuit, ParameterVector, Parameter
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qc4pyscf.operator import Second_Quantization, Excitation, RDM
from qc4pyscf.tools import Transpiler, Evaluate
from qc4pyscf.ansatz import Initial
from scipy.optimize import minimize
import numpy as np
import time


mapper_code = {'jordan_wigner': 'Jordan-Wigner'}


class UCC:
    def __init__(self, 
                 mol                , mf                        , 
                 ex_code='sd'       , mapping='jordan_wigner'   , 
                 cd_acc=1e-6        , max_iteration=200000      , spin_symm=True,
                 amplitudes=[]):
        # Input Values
        self.mol            = mol           # gto.Mole Object
        self.mf             = mf            # PySCF's scf or dft Object
        self.ex_code        = ex_code       # Excitation Code: (sdt)
        self.mapping        = mapping       # Mapping: (jordan_wignar)
        self.cd_acc         = cd_acc        # Cholesty Decomposition Accuracy
        self.max_iteration  = max_iteration # Minimizer Maximum Iteration
        self.amplitudes     = amplitudes    # Amplitudes of Excitation Evolvers

        # Saving Values
        self.H              = None          # 2nd Quantized Hamiltonian
        self.ansatz         = None          # Ansatz Quantum Circuit
        self.energy         = 0             # Energy at Last Cycle
        self.op_ind         = None          # Excitation Operator Indices
        self.op_pool        = None          # Qubit Operator Pools

        # Auto-Loaded Values
        self.N_amplitudes   = 0                                                     # Number of Amplitudes
        self.N_orb          = len(self.mol.intor('int1e_ovlp'))                     # Number of Spacial Orbital
        self.N_alpha        = self.mol.nelec[0]                                     # Number of Alpha-Spin Electron
        self.N_beta         = self.mol.nelec[1]                                     # Number of Beta-Spin Electron
        self.spin_symm      = spin_symm if self.N_alpha == self.N_beta else False   # Excitation Spin Symmetry (It changes into False if the system is open shell)
        self.tp_obj         = Transpiler.transpile(mol, mf)                         # Transpiled Object for 2nd Quantization
        

    def build(self):
        self.H          = Second_Quantization.gen_H(self.tp_obj, self.mapping, self.cd_acc)
        self.gen_ops()
        self.ansatz     = self.gen_ansatz(self.op_pool)
        print("Build Done")

    
    def gen_ops(self):
        excitation_ops      = Excitation.UCC(self.ex_code, self.N_orb, self.N_alpha, self.N_beta, self.mapping, self.spin_symm)
        self.op_pool        = excitation_ops[0]
        self.op_ind         = excitation_ops[1]
        self.N_amplitudes   = len(self.op_pool)
    

    def gen_ansatz(self, op_pool):
            ansatz          = Initial.initial(self.N_orb, self.N_alpha, self.N_beta)
            cc_amps         = ParameterVector("CC Amp.", len(op_pool))
            for ind, op in enumerate(op_pool):
                op.coeffs  *= 1.0j
                op_evolver  = PauliEvolutionGate(op, time=cc_amps[ind])
                ansatz.append(op_evolver, range(2 * self.N_orb))

            return ansatz
    

    def run(self, cost_func, Estimator, minimize_algorithm='COBYLA', on_ansatz_ftn=None):
        # [Input]      cost_func     : Cost function. Input will be (parameters, ansatz, H, estimator)
        # [Input]      Estimator     : Qiskit estimator object
        # [Input] minimize_algorithm : Minimize algorithm for scipy minimizer
        # [Input]    on_ansatz_ftn   : Function for ansatz operated between creation of the ansatz and operation. e.g.) passmanager

        print("=========================================================="
              f"\n     Minimize Algorithm     : {minimize_algorithm}"
              f"\n  Max Minimizing Iteration  : {self.max_iteration}"
              f"\n      Excitation Code       : {self.ex_code.upper()}"
              f"\n     Mapping Algorithm      : {mapper_code[self.mapping]}"
              f"\n        C.D. Accuracy       : {self.cd_acc:.1e}"
              f"\n     Initial Conditions     : {"True" if len(self.amplitudes) == self.N_amplitudes else "False"}"
              "\n==========================================================")

        time_start = time.time()
        cc_amps         = []
        if len(self.amplitudes) == self.N_amplitudes:
            cc_amps     = self.amplitudes
        else: 
            cc_amps     = [0.0] * self.N_amplitudes
            if len(self.amplitudes) != 0:
                print("Initial Amplitude Application Failed. Reset Amplitudes")

        if on_ansatz_ftn != None:
            self.ansatz = on_ansatz_ftn(self.ansatz)

        res             = minimize( cost_func, 
                                    cc_amps, 
                                    args    = (self.ansatz, self.H, Estimator), 
                                    method  = minimize_algorithm, 
                                    options = {'maxiter': self.max_iteration})
        self.energy     = getattr(res,'fun')
        self.amplitudes = getattr(res,'x').tolist()

        time_end = time.time()
        print(f"Calculated Energy: {self.energy} Hartree | Minimize Success: {getattr(res,'success')} | Computation Time: {time_end - time_start:.3f}s")
        return self.energy
    
    def energy_tot(self, cost_func, Estimator, on_ansatz_ftn=None):
        energy = Evaluate.expectation(cost_func, Estimator, self.ansatz, self.H, self.amplitudes, on_ansatz_ftn)
        return energy


    def quasi_distribution(self, Sampler, shot=100000):
        quasi_dists = Evaluate.quasi_distribution(Sampler, self.ansatz, self.amplitudes, shot)
        return quasi_dists
    

    def make_rdm1(self, cost_func, Estimator):
        reduced_density_matrices = RDM.make_rdm1(cost_func, Estimator, self.ansatz, self.amplitudes, self.N_orb, self.mapping)
        return reduced_density_matrices


    def make_rdm1s(self, cost_func, Estimator):
        reduced_density_matrices = RDM.make_rdm1s(cost_func, Estimator, self.ansatz, self.amplitudes, self.N_orb, self.mapping)
        return reduced_density_matrices