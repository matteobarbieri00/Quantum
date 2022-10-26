from qutip import *
import random
import numpy as np
import matplotlib.pyplot as plt

bra_0 = Qobj([[1, 0]])
bra_1 = Qobj([[0, 1]])

P_00 = Qobj([[1, 0], [0, 0]])
P_11 = Qobj([[0, 0], [0, 1]])
P_01 = Qobj([[0, 1], [0, 0]])
P_10 = Qobj([[0, 0], [0, 1]])
cnot_12 = tensor(P_00, qeye(2), qeye(2)) + tensor(P_11, sigmax(), qeye(2))
cnot_13 = tensor(P_00, qeye(2), qeye(2)) + tensor(P_11, qeye(2), sigmax())
cnot_14 = tensor(P_00, qeye(2), qeye(2), qeye(2), qeye(2)) + \
    tensor(P_11, qeye(2), qeye(2), sigmax(), qeye(2))
cnot_24 = tensor(qeye(2), P_00, qeye(2), qeye(2), qeye(2)) + \
    tensor(qeye(2), P_11, qeye(2), sigmax(), qeye(2))
cnot_25 = tensor(qeye(2), P_00, qeye(2), qeye(2), qeye(2)) + \
    tensor(qeye(2), P_11, qeye(2), qeye(2), sigmax())
cnot_35 = tensor(qeye(2), qeye(2), P_00, qeye(2), qeye(2)) + \
    tensor(qeye(2), qeye(2), P_11, qeye(2), sigmax())
X_1 = tensor(sigmax(), qeye(2), qeye(2))
X_2 = tensor(qeye(2), sigmax(), qeye(2))
X_3 = tensor(qeye(2), qeye(2), sigmax())
X_1_tot = tensor(sigmax(), qeye(2), qeye(2), qeye(2), qeye(2))
X_2_tot = tensor(qeye(2), sigmax(), qeye(2), qeye(2), qeye(2))
X_3_tot = tensor(qeye(2), qeye(2), sigmax(), qeye(2), qeye(2))
X_4_tot = tensor(qeye(2), qeye(2), qeye(2), sigmax(), qeye(2))
X_5_tot = tensor(qeye(2), qeye(2), qeye(2), qeye(2), sigmax())
P2_00 = tensor(basis(2, 0), basis(2, 0)) * tensor(bra_0, bra_0)
P2_10 = tensor(basis(2, 1), basis(2, 0)) * tensor(bra_1, bra_0)
P2_01 = tensor(basis(2, 0), basis(2, 1)) * tensor(bra_0, bra_1)
P2_11 = tensor(basis(2, 1), basis(2, 1)) * tensor(bra_1, bra_1)


def find_probs(select):
    if select == 0:
        return tensor(qeye(2), qeye(2), qeye(2), P2_00)
    elif select == 1:
        return tensor(qeye(2), qeye(2), qeye(2), P2_10)
    elif select == 2:
        return tensor(qeye(2), qeye(2), qeye(2), P2_01)
    elif select == 3:
        return tensor(qeye(2), qeye(2), qeye(2), P2_11)


fidelity_not_encoded = []
fidelity_not_encoded_errors = []
fidelity_encoded = []
fidelity_encoded_errors = []
probabilities = []
shots = 1000
for i in range(101):
    p = i/100
    fidelities_ne = []
    probabilities.append(p)
    for t in range(shots):
        #theta = random.uniform(0, np.pi)
        #phi = np.pi/2
        psi_ = basis(2, 0)
        phi_ = basis(2, 0)
        # psi_ = np.cos(theta / 2) * basis(2, 0) + \
        #    np.exp(complex(0, phi)) * basis(2, 1)
        # phi_ = np.cos(theta / 2) * basis(2, 0) + \
        #    np.exp(complex(0, phi)) * basis(2, 1)
        if random.random() <= p:
            phi_ = sigmax() * phi_
        fidelities_ne.append(metrics.fidelity(phi_, psi_))
    fidelity_not_encoded.append(np.mean(fidelities_ne))
    fidelity_not_encoded_errors.append(np.std(fidelities_ne)/np.sqrt(shots))

    fidelities_e = []
    for t in range(shots):
        psi = basis(2, 0)
        #theta = random.uniform(0, np.pi)
        #phi = np.pi/2
        # psi = np.cos(theta / 2) * basis(2, 0) + \
        #    np.exp(complex(0, phi)) * basis(2, 1)
        # encoding the system and introducing the ancilla system
        s_1 = basis(2, 0)
        s_2 = basis(2, 0)
        sys = tensor(psi, s_1, s_2)
        sys_encoded = cnot_12 * cnot_13 * sys
        a_0 = basis(2, 0)
        a_1 = basis(2, 0)
        tot_sys = tensor(sys_encoded, a_0, a_1)
        tot_sys_ = tensor(sys_encoded, a_0, a_1)

        if random.random() < p:
            tot_sys = X_1_tot * tot_sys
        if random.random() < p:
            tot_sys = X_2_tot * tot_sys
        if random.random() < p:
            tot_sys = X_3_tot * tot_sys

        # detection
        tot_sys = cnot_14 * tot_sys
        tot_sys = cnot_24 * tot_sys
        tot_sys = cnot_25 * tot_sys
        tot_sys = cnot_35 * tot_sys

        # correction
        p_detected = []
        for i in range(4):
            Proj = find_probs(i)
            p_detected.append(tot_sys.dag() * Proj * tot_sys)

        scalar_1 = Qobj([[1]])

        if p_detected[0] == scalar_1:
            pass
        elif p_detected[1] == scalar_1:
            tot_sys = X_1_tot * tot_sys
            tot_sys = X_4_tot * tot_sys
        elif p_detected[2] == scalar_1:
            tot_sys = X_3_tot * tot_sys
            tot_sys = X_5_tot * tot_sys
        elif p_detected[3] == scalar_1:
            tot_sys = X_2_tot * tot_sys
            tot_sys = X_4_tot * tot_sys
            tot_sys = X_5_tot * tot_sys

        fidelities_e.append(metrics.fidelity(tot_sys_, tot_sys))
    fidelity_encoded.append(np.mean(fidelities_e))
    # fidelity_encoded_errors.append(np.sqrt(shots))
    fidelity_encoded_errors.append(np.std(fidelities_e)/np.sqrt(shots))

with open('NotEncoded.txt', 'w') as file:
    for i in range(101):
        file.write(str(probabilities[i]))
        file.write('\t')
        file.write(str(fidelity_not_encoded[i]))
        file.write('\t')
        file.write(str(fidelity_not_encoded_errors[i]))
        file.write('\n')

with open('Encoded.txt', 'w') as file:
    for i in range(101):
        file.write(str(probabilities[i]))
        file.write('\t')
        file.write(str(fidelity_encoded[i]))
        file.write('\t')
        file.write(str(fidelity_encoded_errors[i]))
        file.write('\n')

#plt.plot(probabilities, fidelity_not_encoded,fmt='.-')
#plt.plot(probabilities, fidelity_encoded, fmt='.-')


def f(x):
    return 1-x


def f_(x):
    return (1-x)*(1-x)*(1-x) + 3*x*(1-x)*(1-x)


x = np.linspace(0, 1, 100)

plt.errorbar(probabilities, fidelity_not_encoded,
             fidelity_not_encoded_errors, fmt='.-')
plt.errorbar(probabilities, fidelity_encoded,
             fidelity_encoded_errors, fmt='.-')
plt.legend(['Not encoded', 'Encoded'])
plt.plot(x, f(x), color='red')
plt.plot(x, f_(x), color='green')


plt.grid()

plt.xlabel('Probability')
plt.ylabel('Fidelity')
plt.show()
