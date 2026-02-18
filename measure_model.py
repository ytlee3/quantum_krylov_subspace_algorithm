
from scipy.sparse import csr_matrix, kron
import numpy as np 
from qiskit.quantum_info import SparsePauliOp
import scipy as sp 

### CX XZ CX = YY ###
### CX YZ CX = -XY ###
### CX XI CX = XX ###
### CX YI CX = YX ###
### CX XX CX = XI ###
### CX YX CX = YI ###

y_rotate = csr_matrix([[1/np.sqrt(2), -1j*1/np.sqrt(2)],[-1j*1/np.sqrt(2), 1/np.sqrt(2)]])
H = csr_matrix([[1,1],[1,-1]]/np.sqrt(2))
I = csr_matrix([[1,0],[0,1]])

def tensor_product(matrix, times):
    '''
    Do tensor product multiple times
    Args: 
        matrix (array): target matrix
        times (int): how many times you want to multiply it 
    Returns:
        result (matrix)
    '''
    if times == 0:
        return csr_matrix([1])
    else:
        result = matrix
        for _ in range(times-1):
            result = kron(result, matrix)
        return result

def hubbard_edges(n): 
    '''
    Getting the edge of the 1D Hubbard model with spin ladder configuration 
    0-3-4-7 spin up
    | | | | 
    1-2-5-6 spin down 
    up_edges (0,3)(3,4)(4,7)
    down_edges (1,2)(2,5)(5,6)
    site_edges (0,1)(2,3)(4,5)(6,7)
    Args: 
        n (int): size of system
    Return:
        up_edges(list), down_edges(list), site_edges (list)
    '''
    up_edges, down_edges, site_edges = [],[],[]
    for i in range(n//2):
        if i+1 < n:
            site_edges.append((i*2,i*2+1))
    up, down  = [], []
    for i in range(n//4):
        up.append(4*i), up.append(4*i+3)
        down.append(4*i+1), down.append(4*i+2)

    up_edges= [(up[i], up[i+1]) for i in range(len(up)-1)]
    down_edges=[(down[i], down[i+1]) for i in range(len(down)-1)]
    return up_edges, down_edges, site_edges


def transform_XX(paulis, n):
    new_paulis = ''
    for i in range(n):
        if i % 2 == 0: 
            new_paulis += paulis[i]
        else: 
            if paulis[i] == 'Z':
                new_paulis += 'I'
            else:
                new_paulis += 'Z'
    new_paulis += 'Z'
    return new_paulis

def transform_ZZ(paulis, n):
    new_paulis = ''
    for i in range(n):
        if i % 2 == 0: 
            new_paulis += paulis[i]
        else: 
            new_paulis += 'Z'
    new_paulis += 'Z'
    return new_paulis


def measure_XX_ov(up_edges, down_edges, state, n):
    global H
    '''
    Evaluating the <XX> term, <YY> equals to <YY> based on symmtry. 
    The overlap term is <IIIIIIIII>, but the final measurement can be combined with the final control gate as 'IZIZIZIZZ' 
    The real (imaginary) part is by measuring in X (Y) basis 
    Args: 
        up_edges (list): index of spin up coupling 
        down_edges (list): index of spin down coupling 
        state (matrix): wavefunction
        n (int): size of the system
    Return:
        Re(<XX>)+Im(<XX>), Re(<S>)+Im(<S>)
    '''
    H_XX_rotate = SparsePauliOp.from_sparse_list([('ZZ', e, 1) for e in up_edges], num_qubits=n)
    H_XX_rotate += SparsePauliOp.from_sparse_list([('ZZ', e, 1) for e in down_edges], num_qubits=n)
    Paulis_XX = [transform_XX(str(H_XX_rotate[i].paulis[0]), n) for i in range(len(H_XX_rotate))]
    Paulis_XX = [SparsePauliOp(i, 1).to_matrix(sparse=True) for i in Paulis_XX]
    ### real part: control in X others in X ####
    real_state = tensor_product(H, n+1) @ state ## H rotate
    ### real part: control in Y others in X ####
    imag_state = kron(tensor_product(H, n), y_rotate) @ state ## rx(np.pi/2) rotate 
    real_XX, imag_XX = 0, 0
    for i in Paulis_XX:
        real_XX += (real_state.getH() @ i @ real_state).toarray()[0]
        imag_XX += (imag_state.getH() @ i @ imag_state).toarray()[0]
    ### overlap ### 
    overlap_paulis = SparsePauliOp('IZIZIZIZZ',1).to_matrix(sparse=True)
    real_ov_state = tensor_product(H,n+1) @ state
    real_overlap = (real_ov_state.getH() @ overlap_paulis @ real_ov_state).toarray()[0]
    ### 
    imag_ov_state = kron(tensor_product(H,n), y_rotate) @ state #
    imag_overlap = (imag_ov_state.getH() @ overlap_paulis @ imag_ov_state).toarray()[0]
    
    return np.real(real_XX)[0]+1j*np.real(imag_XX)[0], np.real(real_overlap)[0]+1j*np.real(imag_overlap)[0]


def measure_ZZ(site_edges, state, n):
    '''
    Evaluating the <ZZ> term
    for half filling condition, <IZ>+<ZI>= 0
    The real (imaginary) part is by measuring in X (Y) basis 
    Args: 
        site_edges (list): index of site coupling (ZZ)
        state (matrix): wavefunction
        n (int): size of the system
    Return:
        Re(<ZZ>)+Im(<ZZ>)
    '''
    global y_rotate
    H_ZZ = SparsePauliOp.from_sparse_list([('ZZ', e, 1) for e in site_edges], num_qubits=n)
    Paulis_ZZ = [transform_ZZ(str(H_ZZ[i].paulis[0]), n) for i in range(len(H_ZZ))]
    Paulis_ZZ  = [SparsePauliOp(i,1).to_matrix(sparse=True) for i in Paulis_ZZ]
    real_ZZ, imag_ZZ = 0, 0
    HI = kron(I,H)
    YI = kron(I,y_rotate)
    # #assymmetry 
    for i, j in enumerate(Paulis_ZZ):
        # for X is IXIXIXZYY, IXIXZYIXY,.....
    #    YI = kron(I, y_rotate)
        XI_1, XI_2= tensor_product(HI, i), tensor_product(HI, len(Paulis_ZZ)-1-i)
        real_basis_rotate = kron(kron(kron(XI_2, YI), XI_1), y_rotate)
        real_rotated_state = real_basis_rotate @ state
        real_zz = (real_rotated_state.getH() @ j @ real_rotated_state).toarray()[0]
        real_ZZ += real_zz
        imag_basis_rotate = kron(kron(kron(XI_2, YI), XI_1), H)
        imag_rotated_state = imag_basis_rotate @ state
        imag_zz = (imag_rotated_state.getH() @ j @ imag_rotated_state).toarray()[0]
        ### is minus since ### CX YZ CX = -XY 
        imag_ZZ -= imag_zz
    return np.real(real_ZZ[0])+1j*np.real(imag_ZZ[0])


def precondition_diag(krylov_d, H_kry, S_kry):
    '''
    Precondition for stablizing the diagonalization
    Args: 
        krylov_d (int): Krylov dimension 
        H_kry (matrix/array): Krylov Hamiltonian 
        S_kry (matrix/array): Overlap matrix
    Return:
        results (array): ground state energy for each Krylov dimension
        [val, vecs, H_reg, S_reg]: eigenvalues, eigenvectors, regulated-H, regulated-S
    '''
    results =[]
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
    return results, [vals, vecs, H_reg, S_reg]