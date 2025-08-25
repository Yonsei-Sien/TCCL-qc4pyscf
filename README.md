# 🚀 qc4pyscf

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
git clone https://github.com/Yonsei-Sien/qc4pyscf.git
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
123
```
---

## 🧩 Module Reference  

### qc4pyscf.ansatz.Initial.initial(N_orb, N_alpha, N_beta)
Construct a Hartree–Fock reference state.

Arguments
- N_orb (int) – Number of spatial orbitals.
- N_alpha (int) – Number of α-spin electrons.
- N_beta (int) – Number of β-spin electrons.

Returns
- QuantumCircuit – Circuit with X gates applied to occupied orbitals.

---

### qc4pyscf.ansatz.UCC.UCC (class)

UCC(mol, mf, ex_code='sd', mapping='jordan_wigner',
    cd_acc=1e-6, max_iteration=200000,
    spin_symm=True, amplitudes=[])
    
< Arguments >
- mol (pyscf.gto.Mole) – Molecule object.
- mf (pyscf.scf.HF or DFT) – PySCF mean-field solution.
- ex_code (str) – Excitation type ('s', 'd', 'sd').
- mapping (str) – Fermion-to-qubit mapping ('jordan_wigner').
- cd_acc (float) – Cholesky decomposition accuracy.
- max_iteration (int) – Optimizer max iterations.
- spin_symm (bool) – Enforce spin symmetry (auto-disabled for open-shell).
- amplitudes (list[float]) – Initial guess for amplitudes.

< Attributes >

- H (SparsePauliOp) – Hamiltonian.
- ansatz (QuantumCircuit) – UCC ansatz circuit.
- energy (float) – Final energy.
- op_pool (list) – Operator pool.
- op_ind (list) – Operator indices.
- amplitudes (list[float]) – Excitation Amplitudes.

< Methods >

- build() – Build Hamiltonian and ansatz.
- run(cost_func, Estimator, minimize_algorithm='COBYLA', on_ansatz_ftn=None) – Optimize amplitudes.
- energy_tot(cost_func, Estimator, on_ansatz_ftn=None) – Compute total energy.
- quasi_distribution(Sampler, shot=100000) – Sample distribution.
- make_rdm1(...), make_rdm1s(...) – Build reduced density matrices.

---

### ADAPT_VQE.ADAPT_VQE (class)

ADAPT_VQE(mol, mf, ex_code='sd', mapping='jordan_wigner',
          cd_acc=1e-6, max_cycle=50, max_iteration=200000,
          energy_conv=1e-6, grad_conv=1e-2, conv_type=0,
          initial=None, chkfile=None, spin_symm=True)
          
< Arguments >

- mol, mf – PySCF molecule and mean-field objects.
- ex_code (str) – Excitation type ('s', 'd', 'sd').
- mapping (str) – Fermion-to-qubit mapping.
- cd_acc (float) – Cholesky accuracy.
- max_cycle (int) – Maximum ADAPT iterations.
- max_iteration (int) – Max optimizer iterations.
- energy_conv (float) – Energy convergence threshold.
- grad_conv (float) – Gradient convergence threshold.
- conv_type (int) – Convergence check type (0: grad, 1: energy, 2: both).
- initial (str) – Path to saved JSON initial ansatz.
- chkfile (str) – Path to save intermediate results.
- spin_symm (bool) – Spin symmetry enforcement.

< Attributes >

- ansatz (QuantumCircuit) – Current ansatz.
- op_pool (list) – Operator pool.
- op_ind (list) – Operator indices.
- energy (float) – Current energy.
- grad_norm (float) – Norm of gradient.
- amplitudes (list[float]) – Optimized amplitudes.
- op_used (list) – Indices of selected operators.
- is_converged (bool) – Whether algorithm converged.

< Methods >

- build() – Build Hamiltonian, operators, and ansatz.
- run(cost_func, Estimator, minimize_algorithm='BFGS', on_ansatz_ftn=None) – Run ADAPT-VQE loop.
- save(dir) – Save amplitudes/operators to JSON.
- energy_tot(...), quasi_distribution(...) – Evaluate results.
- make_rdm1(...), make_rdm1s(...) – Build reduced density matrices

---

### qc4pyscf.operator

#### Second_Quantization.gen_H(qc_mf, mapper='jordan_wigner', cd_acc=1e-6)
Build Hamiltonian from PySCF mean-field.

#### Second_Quantization.creators_destructors(n, mapping)
Generate creation/annihilation operators.

#### Excitation.UCC(...), Excitation.ADAPT_VQE(...)
Build excitation operator pools.

#### Commute.get_commutator(A, B), get_commutators(A, Bs)
Compute commutators.

#### RDM.make_rdm1(...), make_rdm1s(...), mo2ao(...)
Build reduced density matrices.

### qc4pyscf.tools

#### Evaluate.expectation(cost_func, Estimator, Ansatz, H, parameter=None, on_ansatz_ftn=None)
Evaluate expectation value.

#### Evaluate.quasi_distribution(Sampler, Ansatz, parameter, shot=100000, dir=None)
Sample quasi-distribution (saves JSON if dir is given).

#### Transpiler.transpile(mol, mf, is_chkfile=False) (class)
Convert PySCF AO/MO integrals into Hamiltonian form.

Attributes: mo_h1e, mo_h2e, e0_core, spin.
