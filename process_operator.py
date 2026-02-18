from qiskit.circuit.library import SwapGate
from scipy.sparse import csr_matrix, kron
from qiskit.quantum_info.operators import Operator
import numpy as np 

cnot_0 = csr_matrix(np.array([[0,0,1,0],[0,1,0,0],[1,0,0,0],[0,0,0,1]]))
cnot = csr_matrix(np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]]))
cnot_inverse = csr_matrix(np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]))
swap = csr_matrix(Operator(SwapGate()).data)
##### real/imaginary part of XX ####
# rx(np.pi/2) gate
y_rotate = csr_matrix([[1/np.sqrt(2), -1j*1/np.sqrt(2)],[-1j*1/np.sqrt(2), 1/np.sqrt(2)]])
#### real and imanginary part share the same pauli operator after basis chagne
X = csr_matrix(np.array([[0,1],[1,0]]))
Y = csr_matrix(np.array([[0,-1j],[1j,0]]))
Z = csr_matrix(np.array([[1,0],[0,-1]]))
I = csr_matrix([[1,0],[0,1]])
H = csr_matrix([[1,1],[1,-1]]/np.sqrt(2))
ZZ = kron(Z, Z)
XX = kron(X, X)
II = kron(I, I)

Id_8 = []
for i in range(8):
    Id_8.append(csr_matrix(np.eye(2**i)))
cx_0_8 = csr_matrix(np.eye(2**9))
for i in range(7):
    I1 = Id_8[i]
    I2 = Id_8[7-i]
    swap_net = kron(kron(I2,swap),I1)
    cx_0_8 = swap_net @ cx_0_8
cx_0_8 = kron(cnot, Id_8[7]) @ cx_0_8
for i in range(7):
    I1 = Id_8[i+1]
    I2 = Id_8[7-(i+1)]
    swap_net = kron(kron(I1,swap),I2)
    cx_0_8 = swap_net @ cx_0_8


def get_cx(i, inverse): 
    '''
    get the CNOT gate between ith and i+1 th qubit in the 8-qubit case 
    Args: 
        i (int): index of the qubit
        inverse (bool): either the CNOT is ith to i+1th or i+1th to ith
    Return:
        cx (csr_matrix): CNOT between ith and i+1th
    '''
    global Id_8, cnot, cnot_inverse
    I1 = Id_8[i]
    I2 = Id_8[7-i]
    if inverse: 
        cx = kron(kron(I2, cnot_inverse), I1)
    else: 
        cx = kron(kron(I2, cnot), I1 )
    return cx

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
    
def get_cx_general(target_gate, n): 
    '''
    getting general controlled gate across n qubits by applying swap net 
    Args: 
        target_gate (csr_matrix): controlled gate 
        n (int) : number of qubits acrossed
    Returns:
        controlled_gate (matrix)
    '''
    spacing = np.abs(n-2)
    controlled_gate = csr_matrix(np.eye(2**n))
    for i in range(int(spacing)): 
        I1, I2 = csr_matrix(np.eye(2**i)), csr_matrix(np.eye(2**(n-i-2)))
        swap_layer = kron(kron(I2, swap), I1) 
        controlled_gate = swap_layer @ controlled_gate
    controlled_gate = kron(target_gate, tensor_product(I,spacing)) @ controlled_gate 
    for i in range(int(spacing)):
        I1, I2 = csr_matrix(np.eye(2**(i+1))), csr_matrix(np.eye(2**(n-i-3)))
        swap_layer = kron(kron(I1, swap), I2)
        controlled_gate = swap_layer @ controlled_gate
    return controlled_gate