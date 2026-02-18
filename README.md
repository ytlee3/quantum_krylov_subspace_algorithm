
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

By measuring the $\tilde{H_{jk}} = \langle \psi_j | H | \psi_k \rangle$ and the overlap integral $\tilde{S_{jk}} = \langle \psi_j | \psi_k \rangle$, we transform the $\langle \psi | H | \psi \rangle$ problem into a classically solvable generalized eigenvalue problem $\tilde{H} c = E \tilde{S} c$, whose dimension equals the number of Krylov basis states.

## Quantum circuit and derivation 

The quantum circuit for measuring the  $\tilde{H_{jk}} = \langle \psi_j | H | \psi_k \rangle$ and $\tilde{S_{jk}} = \langle \psi_j | \psi_k \rangle$ are from the paper: https://doi.org/10.1038/s41467-025-59716-z.  Under exact time evolution, 

$$
\tilde{H} = \langle \psi_j | H | \psi_k \rangle = \langle \psi | U^{ijtH} H U^{-iktH} |\psi \rangle= \langle \psi | H U^{-i(k-j)tH} |\psi \rangle 
$$

and the overlapt integral $\tilde{S}$ can also descirbed as $\langle \psi_j| \psi_k \rangle = \langle \psi|U^{-i(k-j)Ht}|\psi \rangle$, which is the measurement of $I^{\otimes{N}}$ for the $\tilde{H}$
The overall quantum circuit 

<p align="center">
<img width="499" height="153" alt="circuit" src="https://github.com/user-attachments/assets/5e30a7bd-3c65-44b8-be1b-7f26cbdd58a3" />

</p>

where measuring in the $X$ ($Y$) basis yields the real (imaginary) part of the results $\langle O U_{k-j} \rangle$. 

Let's derive the wavefunction step by step 

1. After applying Hadamard gate, we have $|\psi \rangle = \frac{1}{\sqrt{2}} (|0\rangle |0\rangle^n + |1\rangle |0\rangle^n)$
2. Preparing the target state on $|1\rangle$ (ancilia), we have $|\psi \rangle = \frac{1}{\sqrt{2}} (|0\rangle |0\rangle^n + |1\rangle |\psi_0\rangle)$ 
3. Apply time-evolution operator: $|\psi \rangle =  \frac{1}{\sqrt{2}} (U_{k-j} |0\rangle |0\rangle^n +U_{k-j}   |1\rangle |\psi_0\rangle) =  \frac{1}{\sqrt{2}} (e^{i\phi}|0\rangle |0\rangle^n +U_{k-j}   |1\rangle |\psi_0\rangle) $.  NOTE: $U_{k-j}|0\rangle^n = e^{i \phi} |0\rangle^n$

This phase factor is easy to calculate classically

4. Preparing the target state on $|0\rangle$ (ancilia), we have $|\psi \rangle = \frac{1}{\sqrt{2}} ( e^{i\phi}|0\rangle |\psi_0\rangle +U_{k-j} |1\rangle |\psi_0\rangle) $
 

## Keep in mind ...
