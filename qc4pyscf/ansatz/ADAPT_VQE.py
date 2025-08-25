from qiskit.circuit import QuantumCircuit, ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qc4pyscf.operator import Second_Quantization, Excitation, Commute, RDM
from qc4pyscf.tools import Transpiler, Evaluate
from qc4pyscf.ansatz import Initial
from scipy.optimize import minimize
import numpy as np
import time
import json


mapper_code = {'jordan_wigner': 'Jordan-Wigner'}


class ADAPT_VQE:
    def __init__(self, 
                 mol                , mf                        , 
                 ex_code='sd'       , mapping='jordan_wigner'   , 
                 cd_acc=1e-6        , max_cycle=50              , max_iteration=200000,
                 energy_conv=1e-6   , grad_conv=1e-2            , conv_type=0,
                 initial=None       , chkfile=None              , spin_symm=True):
        # Input Values
        self.mol            = mol           # gto.Mole Object
        self.mf             = mf            # PySCF's scf or dft Object
        self.ex_code        = ex_code       # Excitation Code: (sdt)
        self.mapping        = mapping       # Mapping: (jordan_wignar)
        self.cd_acc         = cd_acc        # Cholesky Decomposition Accuracy
        self.energy_conv    = energy_conv   # Energy Converge Tolerance
        self.grad_conv      = grad_conv     # 2-Norm of Gradient Vector Converge Tolerance
        self.conv_type      = conv_type     # 0: Grad | 1: Energy | 2: Grad and Energy
        self.max_cycle      = max_cycle     # ADAPT-VQE Maximum Cycle
        self.max_iteration  = max_iteration # Minimizer Maximum Iteration
        self.initial        = initial       # Initial Ansatz JSON Directory
        self.chkfile        = chkfile       # Chkfile Directory
        self.spin_symm      = spin_symm     # Spin Symmetry

        # Saving Values
        self.is_converged   = False         # Convergence of Algorithm
        self.ansatz         = None          # Ansatz Quantum Circuit
        self.H              = None          # 2nd Quantized Hamiltonian
        self.op_pool        = None          # Qubit Operator Pools
        self.op_ind         = None          # Fermionic Indices of Operator Pools
        self.op_comm        = None          # H-Commutated Operator Pools (= [H, A])
        self.op_mod_coef    = None          # Coefficient-Modified Excitation Operator Pools
        self.op_used        = None          # Index Data of Used Evolvers
        self.amplitudes     = None          # Amplitudes of Evolvers
        self.energy         = 0             # Energy at Last Cycle
        self.grad_norm      = 0             # 2-Norm of Gradient Vector at Last Cycle

        # Auto-Loaded Values
        self.N_orb          = len(self.mol.intor('int1e_ovlp'))         # Number of Spacial Orbital
        self.N_alpha        = self.mol.nelec[0]                         # Number of Alpha-Spin Electron
        self.N_beta         = self.mol.nelec[1]                         # Number of Beta-Spin Electron
        self.tp_obj         = Transpiler.transpile(self.mol, self.mf)   # Transpiled Object for 2nd Quantization


    def build(self):
        self.H              = Second_Quantization.gen_H(self.tp_obj, self.mapping, self.cd_acc)
        self.gen_ops()
        self.gen_ansatz()
        print(f"Build Done. #Operator in Operator Pool: {len(self.op_pool)}")


    def gen_ansatz(self):
        ansatz = Initial.initial(self.N_orb, self.N_alpha, self.N_beta)    # Initial Circuit (Occ. / Virt.)
        if self.initial != None:
            initial_data = None
            with open(f"{self.initial}.json", 'r') as f:
                initial_data = json.load(f)

            self.amplitudes = initial_data['Amplitudes']
            self.op_used    = initial_data['Used_Operators']
            amp_pv = ParameterVector("ADAPT-VQE Amp.", len(self.op_used))

            for op in range(len(self.op_used)):
                try:
                    used_op_index           = self.op_ind.index(self.op_used[op])
                    target_operator         = self.op_pool[used_op_index]
                    target_operator.coeffs *= 1.0j
                    target_evolver          = PauliEvolutionGate(target_operator, time=amp_pv[op])
                    ansatz.append(target_evolver, range(2*self.N_orb))
                except ValueError:
                    print("Appropriate Excitation Doesn't Exist in Current Operator Pool. Skip the Evolver")
        self.ansatz = ansatz
                    

    def gen_ops(self):
        excitation_ops      = Excitation.ADAPT_VQE(self.ex_code, self.N_orb, self.N_alpha, self.N_beta, self.mapping, self.spin_symm)
        self.op_pool        = excitation_ops[0]
        self.op_ind         = excitation_ops[1]
        self.op_comm        = Commute.get_commutators(self.H, self.op_pool)
        self.op_mod_coef    = self.coefficient_modifier(self.op_pool)


    def coefficient_modifier(self, Ops):
        modified_ops = []
        for op in Ops:
            op.coeffs *= 1.0j
            modified_ops.append(op)
        return modified_ops


    def convergence_checker(self, dE, grad): # 1: converged | 0: not converged
        if self.conv_type == 0: # grad
            return 1 if grad < self.grad_conv else 0
        elif self.conv_type == 1: # dE
            return 1 if dE < self.energy_conv else 0
        else:
            return 1 if (grad < self.grad_conv and dE < self.energy_conv) else 0
        

    def run(self, cost_func, Estimator, minimize_algorithm='BFGS', on_ansatz_ftn=None):
        # [Input]      cost_func     : Cost function. Input will be (parameters, ansatz, H, estimator)
        # [Input]      Estimator     : Qiskit estimator object
        # [Input] minimize_algorithm : Minimize algorithm for scipy minimizer
        # [Input]    on_ansatz_ftn   : Function for ansatz operated between creation of the ansatz and operation. e.g.) passmanager
        
        print("=========================================================="
              f"\n     Minimize Algorithm     : {minimize_algorithm}"
              f"\n  Max Minimizing Iteration  : {self.max_iteration}"
              f"\n     Max ADAPT-VQE Cycle    : {self.max_cycle}"
              f"\n      Excitation Code       : {self.ex_code.upper()}"
              f"\n     Mapping Algorithm      : {mapper_code[self.mapping]}"
              f"\n        C.D. Accuracy       : {self.cd_acc:.2e}"
              f"\n     Converge Tolerance     : E({self.energy_conv:.2e}) | Grad({self.grad_conv:.2e})"
              f"\n      Convergence Type      : {self.conv_type} (0: Grad | 1: Energy | 2: Grad and Energy)"
              f"\n         Checkfile          : {self.chkfile}"
              "\n==========================================================")

        # Initialization
        cycle       = 0
        E_old       = 0
        E_new       = 1
        grad_norm   = 2*self.grad_conv

        amp_pv      = ParameterVector("ADAPT-VQE Amp.", 0)
        amp_lst     = []
        op_used     = []

        zero_observable_checker = []
        zero_observable = SparsePauliOp('I' * (2 * self.N_orb)) * 0
        for ind, op in enumerate(self.op_comm):
            if op == zero_observable:
                op.coeffs = 1
                zero_observable_checker.append(ind)
        time_start_alg  = time.time()

        while not(self.convergence_checker(abs(E_old-E_new), grad_norm)) and self.max_cycle > cycle:
            # 0th: Initial Setting
            time_start = time.time()
            if on_ansatz_ftn != None:
                self.ansatz = on_ansatz_ftn(self.ansatz)
            E_old = E_new

            # 1st: Gradient Calculate
            grads               = cost_func(amp_lst, self.ansatz, self.op_comm, Estimator) # Estimator 1
            for zero_op in zero_observable_checker:
                grads[zero_op] =0
            grad_norm           = np.linalg.norm(grads)
            grad_largest_ind    = grads.argmax()
            op_used.append(self.op_ind[grad_largest_ind])

            # 2nd: State Evolve
            amp_pv.resize(cycle+1)
            amp_lst.append(0)
            target_operator         = self.op_mod_coef[grad_largest_ind]
            target_evolver          = PauliEvolutionGate(target_operator, time=amp_pv[cycle])
            self.ansatz.append(target_evolver, range(2*self.N_orb))

            # 3rd: Energy Minimizing & Proper Parameter Searching 
            if on_ansatz_ftn != None:
                self.ansatz = on_ansatz_ftn(self.ansatz)

            master_cost_func = cost_func
                
            res     = minimize( master_cost_func, 
                                amp_lst, 
                                args=(self.ansatz, self.H, Estimator), 
                                method=minimize_algorithm, 
                                options={'maxiter': self.max_iteration}) # Estimator 2
            E_new   = getattr(res,'fun')
            amp_lst = getattr(res,'x').tolist()
            
            time_end = time.time()
            print(f"#Evolver: {cycle:>5} | Grad_Norm: {grad_norm:.5f} | Used Op.: {self.op_ind[grad_largest_ind]:^32} | Minimize: {"Success" if getattr(res,'success') else "Failed "} | Energy: {E_new:.9f} AU | dE: {E_new-E_old:.9f} AU | Computation Time for Cycle: {time_end - time_start:.3f}s")

            if self.chkfile != None:
                self.amplitudes = amp_lst
                self.op_used    = op_used
                self.save(f"{self.chkfile}")
            cycle += 1

        time_end_alg  = time.time()
        if self.convergence_checker(abs(E_new-E_old), grad_norm):
            self.converged=True
            print(f"System Converged. E = {E_new} Hartree | Total Computation Time: {time_end_alg - time_start_alg:.3f}s")
        else:
            print(f"System Not Converged. E = {E_new} Hartree | Total Computation Time: {time_end_alg - time_start_alg:.3f}s")
        self.energy     = E_new
        self.amplitudes = amp_lst
        self.grad_norm  = grad_norm
        self.op_used    = op_used

        return self.energy
    

    def save(self, dir):
        data = {'Amplitudes':self.amplitudes, 'Used_Operators':self.op_used}
        with open(f"{dir}.json", 'w') as f:
            json.dump(data, f, indent=4)
        print('Save Done')
    

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

