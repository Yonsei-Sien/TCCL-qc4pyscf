import json

def expectation(cost_func, Estimator, Ansatz, H, parameter=None, on_ansatz_ftn=None):
    # Initial Setting
    if on_ansatz_ftn != None:
        Ansatz = on_ansatz_ftn(Ansatz)
            
    # Calculate Energy
    energy = cost_func(Ansatz, H, Estimator) if parameter == None else cost_func(parameter, Ansatz, H, Estimator)
    print(f"Calculated Energy: {energy} Hartree")

    return energy


def quasi_distribution(Sampler, Ansatz, parameter, shot=100000, dir=None):
    Ansatz.measure_all()
    quasi_dists = Sampler.run(Ansatz, parameter, shots=shot).result().quasi_dists[0]

    states = {}
    for key in quasi_dists:
        bin_key         = str(bin(int(key))[2:].zfill(Ansatz.num_qubits))
        states[bin_key] = quasi_dists[key]

    if dir != None:
        with open(f"{dir}.json", 'w') as f:
            json.dump(states, f, indent=4)
        print('Save Done')
    return quasi_dists