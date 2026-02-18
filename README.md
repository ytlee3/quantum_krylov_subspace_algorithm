# Quantum Krylov Subspace Algorithm

Krylov method is a subspace method designed for solving the $Ax=b$ problem. To solve this, $x \in \mathrm{Span} [b, Ab, A^2 b, A^3 b, \ldots]$, and several numerical methods are based on this concept, including Lanczos and conjugate gradient.

## Quantum version 
We aim to solve the ground-state energy of a given Hamiltonian $\langle \psi | H | \psi \rangle$, where we construct the wavefunction as 

$$
|\psi\rangle = \mathrm{span} [ H^j |\psi_0\rangle ], j = 0,1,2,\ldots, D-1
$$ 

To make it parqametrizable on quantum computers, we use real-time evolution to constrcut the krylov subspace instead of the entire Hamiltonian as follow 

$$
|\psi\rangle = \mathrm{span} [ U^j |\psi_0\rangle ], U=e^{-iHt}, j = 0,1,2,\ldots, D-1
$$ 

By measuring the $\tilde{H_{jk}} = \langle \psi_j | H | \psi_k \rangle$ and the overlap integral $\tilde{S_{jk}} = \langle \psi_j | \psi_k \rangle$, we transform the $\langle \psi | H | \psi \rangle$ problem into a classically solvable generalized eigenvalue problem $\tilde{H} c = E \tilde{S} c$.
