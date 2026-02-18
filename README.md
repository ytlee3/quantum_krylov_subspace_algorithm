# Quantum Krylov Subspace Algorithm

Krylov method is a subspace method designed for solving the $Ax=b$ problem. To solve this, $x \in \mathrm{Span} [b, Ab, A^2 b, A^3 b, \ldots]$, and several numerical methods are based on this concept, including Lanczos and conjugate gradient.

## Quantum version 
We aim to solve the ground-state energy of a given Hamiltonian $\langle \psi | H | \psi \rangle$, where we construct the wavefunction as $$|\psi\rangle = \mathrm{span},[ H^j |\psi_0\rangle ],; j = 0,1,2,\ldots, D-1$$ (Krylov subspace).
