import numpy as np 
from scipy.sparse import csr_matrix, kron 
from measure_model import tensor_product
from qiskit.circuit.library import XXPlusYYGate, CPhaseGate, SwapGate
from qiskit.quantum_info.operators import Operator
import scipy as sp 
from exact import H_exact

def numerical_quantum_krylov_diag(n, initial_config, t, U, n_occ, dt, krylov_d, toeplitz=False): 
    '''
    Numerically tranform the hamiltonian into Krylov basis and diagonalize it 
    Args: 
        n (int): size of the system 
        initial_config (bitstring): bitstring state 
        t (float): hopping term for XX+YY
        U (float): self energy term 
        n_occ (int): occupancy 
        dt (float): time interval for trotterization 
        krylov_d (int): Krylov dimension
        toeplitz (bool): whether conserve the toeplitz structure or not 
    Return:
        numerical_krylov_result (list): ground state energy for each krylov dimension
        H (array), S(array): Krylov H and Krylov S 
    '''   

    psi0 = np.zeros(shape=(2**n,1), dtype=complex)
    psi0[int(initial_config,2)] = 1 
    csr_statevector = csr_matrix(psi0)

    ####### krylov basis with 2nd order trotterization ######
    I = csr_matrix([[1,0],[0,1]])
    site_half = csr_matrix(CPhaseGate(-dt*U/2).to_matrix())
    hop_half = csr_matrix(XXPlusYYGate(-2*t*dt/2).to_matrix())
    hop_full = csr_matrix(XXPlusYYGate(-2*t*dt).to_matrix())
    swap = csr_matrix(Operator(SwapGate()).data)
    swap_layer = tensor_product(swap, n//2)
    site_half_layer = tensor_product(site_half, n//2)
    hop_half_layer = (kron(kron(I,tensor_product(hop_half, n//2-1)),I))
    swap_full_layer =  swap_layer @ kron(kron(I, tensor_product(hop_full, n//2-1)), I) @ swap_layer
    krylov_basis = []
    for _ in range(krylov_d):
        krylov_basis.append(csr_statevector.toarray())
        csr_statevector =  site_half_layer @ (hop_half_layer @ (swap_full_layer @ (hop_half_layer @ (site_half_layer @ csr_statevector))))

    ######## Trasformation to Krylov basis #######
    Ham = H_exact(n, t, U, n_occ)
    v = np.array(krylov_basis).T  # corresponding shape = (dim, krylov), which is (512, 35) # span{u^{itH}|psi>} 
    v = v.reshape(-1,krylov_d)
    H_kry = np.conj(v).T @ Ham @ v 
    S_kry = np.conj(v).T @ v 
    results =[]
    ####### construct the krylov H and S based on toeplitz structure #######
    if toeplitz: 
        H_kry_toeplitz = np.zeros(shape=(krylov_d,krylov_d), dtype=complex)
        S_kry_toeplitz = np.zeros(shape=(krylov_d,krylov_d), dtype=complex)
        for i in range(krylov_d):
            for j in range(krylov_d):
                H_kry_toeplitz[i][j] = H_kry[0,np.abs(i-j)]
                S_kry_toeplitz[i][j] = S_kry[0,np.abs(i-j)]   

        
        for d in range(1,krylov_d+1):
            S_kry_d, H_kry_d = S_kry_toeplitz[:d,:d], H_kry_toeplitz[:d,:d] # d dimension space < D 
            vals, vecs = sp.linalg.eigh(H_kry_d, S_kry_d)
            results.append(min(vals))
        return results, H_kry_toeplitz, S_kry_toeplitz
    ####### construct the krylov H and S with trotterization for each element #######
    ####### NOTE: under trotterization, Krylov H and S won't be toeplitz matrix #######
    else: 
        threshold_rescale = 0.1
        for d in range(1,krylov_d+1):
            S_kry_d, H_kry_d = S_kry[:d,:d], H_kry[:d,:d] # d dimension space < D 
            s_vals, s_vecs = sp.linalg.eigh(S_kry_d) # dim = (d,d)
            s_vecs = s_vecs.T 
            reg_transform = [vec for val, vec in zip(s_vals, s_vecs) if val > d * threshold_rescale*1E-24 ]
            if len(reg_transform) != 0:
                reg_transform = np.array(reg_transform).T # dim = (sub_d, sub_d
                H_reg = reg_transform.T.conj() @ H_kry_d @ reg_transform # project/regularize to subspace 
                S_reg = reg_transform.T.conj() @ S_kry_d @ reg_transform # project/regularize to subspace 
                # Solve generalized eigenvalue problem
                vals, vecs = sp.linalg.eigh(H_reg, S_reg)
                results.append(min(vals))
            else: 
                print("D = "+str(d)+' WHOLE SUBSPACE ILL-CONDITIONED')
                continue
        return results, H_kry, S_kry


    

