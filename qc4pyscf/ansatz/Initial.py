from qiskit.circuit import QuantumCircuit

def initial(N_orb, N_alpha, N_beta):
    circuit = QuantumCircuit(2*N_orb)
    for i in range(N_alpha):
        circuit.x(i)
    for j in range(N_beta):
        circuit.x(j + N_orb)
    return circuit