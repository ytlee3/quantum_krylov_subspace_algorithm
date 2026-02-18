from measure_model import hubbard_edges, tensor_product, measure_XX_ov, measure_ZZ, precondition_diag
from scipy.sparse import csr_matrix, kron 
import numpy as np 
from qiskit.circuit.library import XXPlusYYGate, CPhaseGate, SwapGate
from qiskit.quantum_info.operators import Operator

I = csr_matrix([[1,0],[0,1]])
swap = csr_matrix(Operator(SwapGate()).data)
II = kron(I, I)

def statevector_circuit(n, state, t, U, dt, krylov_d):
    '''
    TESTING: Half filling, <IZ>+<ZI> = 0 is applied, need to modified this for general condition if you want 
    Args: 
        n (int): size of the system 
        state (QuantumCircuit): quantum circuit with H gate (ancilia) and controlled gate for state preparation
        t (float): hopping term for XX+YY
        U (float): self energy term 
        dt (float): time interval for trotterization 
        krylov_d (int): Krylov dimension
    Return:
        result (list): ground state energy for each krylov dimension
        H (array), S(array): Krylov H and Krylov S 
    '''
    up_edges, down_edges, site_edges = hubbard_edges(n)
    # 2nd order trotterization trotterizatoin # 
    site_half_layer = kron(tensor_product(csr_matrix(CPhaseGate(-dt*U/2).to_matrix()), 4),I)
    hop_half_layer =  kron(kron(I, tensor_product(csr_matrix(XXPlusYYGate(-2*t*dt/2).to_matrix()),3)), II)
    swap_layer = kron(tensor_product(swap, 4), I)
    hop_full_layer =  kron(kron(I, tensor_product(csr_matrix(XXPlusYYGate(-2*t*dt).to_matrix()),3)), II)
    H_kry, S_kry = np.zeros(shape=(krylov_d,krylov_d),dtype= np.complex128), np.zeros(shape=(krylov_d,krylov_d),dtype= np.complex128)
    for i in range(krylov_d):
        expval_XX, expval_overlap = measure_XX_ov(up_edges, down_edges, state, n)
        expval_ZZ = measure_ZZ(site_edges, state, n)
        energy = -t* expval_XX + (expval_overlap*4+expval_ZZ)*U/4
        state = hop_full_layer @ (swap_layer @ (hop_half_layer @ (site_half_layer @ state)))
        state = site_half_layer @ (hop_half_layer @ (swap_layer @ state))  
        ### building entire H and S based on toeplitz property
        for j in range(krylov_d):
            for k in range(krylov_d):
                if (j-k) == -i: 
                    H_kry[j][k], S_kry[j][k] = energy, expval_overlap
                elif (j-k) == i:
                    H_kry[j][k], S_kry[j][k] = energy.conj(), expval_overlap.conj()

    result = precondition_diag(krylov_d, H_kry, S_kry)
    return result[0], H_kry, S_kry