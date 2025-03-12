from qiskit.circuit import QuantumCircuit

def initial(nspace_orb, nalpha, nbeta):
    circuit = QuantumCircuit(2*nspace_orb)
    for i in range(nalpha):
        circuit.x(i)
    for j in range(nbeta):
        circuit.x(j + nspace_orb)
    return circuit