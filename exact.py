### exact ###
from more_itertools import distinct_permutations
from measure_model import hubbard_edges
from qiskit.quantum_info import SparsePauliOp
import numpy as np 
import scipy as sp 

def elem(dim, i):
    '''
    making a pure statevector with only one bitstring with P=1 
    Args: 
        dim (int): dimension of the quantum space (2^n)
        i (int): index of the bitstring (int(bitstring, 2))
    Returns:
        vec (array)
    '''
    vec = np.zeros(dim)
    vec[i] = 1.0
    return vec 

def H_exact(n, t, U, n_occ, diagonalize=False):
    '''
    Evaluating the exact ground state by diagonalizing the Hamiltonian 
    Args: 
        n (int): size of the system 
        t (float): hopping term for XX+YY
        U (float): self energy term 
        n_occ (int): occupancy
        diagonalization (bool): whether diagonalize the hamiltonian or not  
    Returns:
        if diagonalize:
            gs_energy (array): ground state energy, gs_psi (array): ground state wavefunction 
        else:
            Ham (sparse_matrix): Hamiltonian
    '''
    up_edges, down_edges, site_edges = hubbard_edges(n)
    # hopping term 
    H = SparsePauliOp.from_sparse_list([('XX', e, -t/2) for e in up_edges], num_qubits=n)
    H += SparsePauliOp.from_sparse_list([('YY', e, -t/2) for e in up_edges], num_qubits=n)
    H += SparsePauliOp.from_sparse_list([('XX', e, -t/2) for e in down_edges], num_qubits=n)
    H += SparsePauliOp.from_sparse_list([('YY', e, -t/2) for e in down_edges], num_qubits=n)
    H += SparsePauliOp.from_sparse_list([('II', e, U/4) for e in site_edges], num_qubits=n)
    H += SparsePauliOp.from_sparse_list([('IZ', e, -U/4) for e in site_edges], num_qubits=n)
    H += SparsePauliOp.from_sparse_list([('ZI', e, -U/4) for e in site_edges], num_qubits=n)
    H += SparsePauliOp.from_sparse_list([('ZZ', e, U/4) for e in site_edges], num_qubits=n)
    Ham = H.to_matrix(sparse=True)
    if diagonalize:
        ### reduce the space based on the occupation ### 
        symm_sec_vecs = np.array([elem(2**n, int(''.join(v), 2)) for v in distinct_permutations('0'*(n - n_occ) + '1'*n_occ)]).T
        tilde_H = symm_sec_vecs.T @ Ham @ symm_sec_vecs
        ## basis are still orthogonal
        gs_energy, gs_psi = sp.sparse.linalg.eigsh(tilde_H, which='SA', k=1)
        return gs_energy, gs_psi
    else: 
        return Ham