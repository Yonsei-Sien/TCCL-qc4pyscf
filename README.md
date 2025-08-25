# 🚀 qc4pyscf
By G. Kim

**qc4pyscf** is a Python toolkit that connects **PySCF** (quantum chemistry) with **Qiskit** (quantum computing).  
It provides ansatz constructors (UCC, ADAPT-VQE), operator utilities, and tools for hybrid quantum–classical simulations.

---

## ✨ Features
- UCC and ADAPT-VQE ansatz generation  
- Automatic Hamiltonian construction from PySCF  
- Excitation operator pools (singles, doubles)  
- Commutator evaluation and operator selection  
- Reduced density matrices (1-RDM, spin-resolved 1-RDM)  
- Qiskit-compatible expectation value and quasi-distribution evaluation  

---

## 📦 Installation

```bash
git clone https://github.com/Yonsei-Sien/TCCL-qc4pyscf.git
cd qc4pyscf
pip install -e .
```

### < Requirements >
- qiskit
- pyscf
- numpy
- scipy

---

## 🔥 Quick Example
```
from qiskit.primitives import StatevectorEstimator as Estimator
from qc4pyscf.ansatz.ADAPT_VQE import ADAPT_VQE
from qc4pyscf.ansatz.UCC import UCC
from pyscf import gto,scf,fci

mol     = gto.M(atom = "H 0.000000 0.000000 0; H 0.000000 0.000000 0.720000", 
                basis = "sto3g",
                charge = 0,
                spin = 0,
                verbose = 0)

# 2. Make HF object and run
mf      = scf.RHF(mol)
E_hf    = mf.kernel()

# Code for reference
mf_fci  = fci.FCI(mf)
E_fci   = mf_fci.kernel()[0]

print(f" HF: {E_hf:.3f} Hartree\nFCI: {E_fci:.3f} Hartree\n\nCorrelation Energy: {E_fci - E_hf:.3f} Hartree")

def cost_func(params, ansatz, H, estimator):
    res = estimator.run([(ansatz, H, params)]).result()
    energy = res[0].data.evs
    return energy

#########<<< UCCSD

qc_mf = UCC(mol           = mol,                # gto.Mole Object
            mf            = mf,                 # PySCF's scf or dft Object
            ex_code       = 'sd',               # Excitation Code: (sd)
            mapping       = 'jordan_wigner',    # Mapping: (jordan_wignar)
            cd_acc        = 1e-6,               # Cholesty Decomposition Accuracy
            max_iteration = 200000,             # Minimizer Maximum Iteration
            amplitudes    = [],                 # Amplitudes of Excitation Evolvers
            spin_symm     = True)               # Excitation Spin Symmetry (It changes into False if the system is open shell)
            
qc_mf.build()

qc_mf.run(  cost_func,                      # Cost function to use
            Estimator(),                    # Quantum estimator
            minimize_algorithm='COBYLA',    # Mimization Algorithm for scipy minimize
            on_ansatz_ftn=None)             # Function for ansatz between building and operating e.g. pm.run()

#########<<< ADAPT-VQE

qc_mf = ADAPT_VQE(  mol           = mol,                # gto.Mole Object
                    mf            = mf,                 # PySCF's scf or dft Object
                    ex_code       = 'sd',               # Excitation Code: (sd)
                    mapping       = 'jordan_wigner',    # Mapping: (jordan_wignar)
                    cd_acc        = 1e-6,               # Cholesky Decomposition Accuracy
                    energy_conv   = 1e-6,               # Energy Converge Tolerance
                    grad_conv     = 1e-2,               # 2-Norm of Gradient Vector Converge Tolerance
                    conv_type     = 0,                  # 0: Grad | 1: Energy | 2: Grad and Energy
                    max_cycle     = 50,                 # ADAPT-VQE Maximum Cycle
                    max_iteration = 200000,             # Minimizer Maximum Iteration
                    initial       = None,               # Initial Ansatz JSON Directory
                    chkfile       = None,               # Chkfile Directory
                    spin_symm     = True)               # Spin Symmetry
            
# 2. Build QC Ingredients
qc_mf.build()

# 3. Run UCC Algorithm
qc_mf.run(  cost_func,                      # Cost function to use
            Estimator(),                    # Quantum estimator
            minimize_algorithm='BFGS',      # Mimization Algorithm for scipy minimize
            on_ansatz_ftn=None)             # Function for ansatz between building and operating e.g. pm.run()

```
---

## 🧩 Module Reference  

### qc4pyscf.ansatz.Initial.initial(N_orb, N_alpha, N_beta)
Construct a Hartree–Fock reference state.

Arguments
- N_orb (int)     – Number of spatial orbitals.
- N_alpha (int)   – Number of α-spin electrons.
- N_beta (int)    – Number of β-spin electrons.

Returns
- QuantumCircuit – Circuit with X gates applied to occupied orbitals.

---

### qc4pyscf.ansatz.UCC.UCC (class)

UCC(mol, mf, ex_code='sd', mapping='jordan_wigner',
    cd_acc=1e-6, max_iteration=200000,
    spin_symm=True, amplitudes=[])
    
< Arguments >
- mol (pyscf.gto.Mole)      – Molecule object.
- mf (pyscf.scf.HF or DFT)  – PySCF HF or KS object.
- ex_code (str)             – Excitation type ('s', 'd', 'sd').
- mapping (str)             – Fermion-to-qubit mapping ('jordan_wigner').
- cd_acc (float)            – Cholesky decomposition accuracy.
- max_iteration (int)       – Optimizer max iterations.
- spin_symm (bool)          – Enforce spin symmetry (auto-disabled for open-shell).
- amplitudes (list[float])  – Initial guess for amplitudes.

< Attributes >

- H (SparsePauliOp)         – Hamiltonian.
- ansatz (QuantumCircuit)   – UCC ansatz circuit.
- energy (float)            – Final energy.
- op_pool (list)            – Operator pool.
- op_ind (list)             – Operator indices.
- amplitudes (list[float])  – Excitation Amplitudes.

< Methods >

- build()                                                                    – Build Hamiltonian and ansatz.
- run(cost_func, Estimator, minimize_algorithm='COBYLA', on_ansatz_ftn=None) – Optimize amplitudes.
- energy_tot(cost_func, Estimator, on_ansatz_ftn=None)                       – Compute total energy.
- quasi_distribution(Sampler, shot=100000)                                   – Sample distribution.
- make_rdm1(...), make_rdm1s(...)                                            – Build reduced density matrices.

---

### ADAPT_VQE.ADAPT_VQE (class)

ADAPT_VQE(mol, mf, ex_code='sd', mapping='jordan_wigner',
          cd_acc=1e-6, max_cycle=50, max_iteration=200000,
          energy_conv=1e-6, grad_conv=1e-2, conv_type=0,
          initial=None, chkfile=None, spin_symm=True)
          
< Arguments >

- mol, mf             – PySCF molecule and (HF or KS) object. objects.
- ex_code (str)       – Excitation type ('s', 'd', 'sd').
- mapping (str)       – Fermion-to-qubit mapping.
- cd_acc (float)      – Cholesky Decomposition accuracy.
- max_cycle (int)     – Maximum ADAPT-VQE iterations.
- max_iteration (int) – Max optimizer iterations.
- energy_conv (float) – Energy convergence threshold.
- grad_conv (float)   – Gradient convergence threshold.
- conv_type (int)     – Convergence check type (0: grad, 1: energy, 2: both).
- initial (str)       – Path to saved JSON initial ansatz.
- chkfile (str)       – Path to save intermediate results.
- spin_symm (bool)    – Spin symmetry enforcement.

< Attributes >

- ansatz (QuantumCircuit)  – Current ansatz.
- op_pool (list)           – Operator pool.
- op_ind (list)            – Operator indices.
- energy (float)           – Current energy.
- grad_norm (float)        – Norm of gradient.
- amplitudes (list[float]) – Optimized amplitudes.
- op_used (list)           – Indices of selected operators.
- is_converged (bool)      – Whether algorithm converged.

< Methods >

- build()                                                                  – Build Hamiltonian, operators, and ansatz.
- run(cost_func, Estimator, minimize_algorithm='BFGS', on_ansatz_ftn=None) – Run ADAPT-VQE loop.
- save(dir)                                                                – Save amplitudes/operators to JSON.
- energy_tot(...), quasi_distribution(...)                                 – Evaluate results.
- make_rdm1(...), make_rdm1s(...)                                          – Build reduced density matrices

---

### qc4pyscf.operator

#### Second_Quantization.gen_H(qc_mf, mapper='jordan_wigner', cd_acc=1e-6)
Build Hamiltonian from PySCF HF or KS object.

#### Second_Quantization.creators_destructors(n, mapping)
Generate creation/annihilation operators.

#### Excitation.UCC(...), Excitation.ADAPT_VQE(...)
Build excitation operator pools.

#### Commute.get_commutator(A, B), get_commutators(A, Bs)
Compute commutators (operators in classical computer).

#### RDM.make_rdm1(...), make_rdm1s(...), mo2ao(...)
Build reduced density matrices.

---

### qc4pyscf.tools

#### Evaluate.expectation(cost_func, Estimator, Ansatz, H, parameter=None, on_ansatz_ftn=None)
Evaluate expectation value.

#### Evaluate.quasi_distribution(Sampler, Ansatz, parameter, shot=100000, dir=None)
Sample quasi-distribution (saves JSON if dir is given).

#### Transpiler.transpile(mol, mf, is_chkfile=False) (class)
Convert PySCF AO/MO integrals into Hamiltonian form.

Attributes: mo_h1e, mo_h2e, e0_core, spin.
