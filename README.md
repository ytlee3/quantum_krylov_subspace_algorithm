
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

5. For measurement in X basis (ancilia), we have to rewrite the ancilia wavefunction as $|+\rangle$ and $|-\rangle$

$|\psi \rangle = \frac{1}{\sqrt{2}} ( e^{i\phi}|0\rangle |\psi_0\rangle +U_{k-j} |1\rangle |\psi_0\rangle) = \frac{1}{2} (|+\rangle (e^{i\phi} |\psi_0\rangle + U_{k-j} |\psi_0\rangle) + |-\rangle (e^{i\phi} |\psi_0\rangle - U_{k-j} |\psi_0\rangle) $ 

$P_+ = \frac{1}{4} (\langle \psi_0 | e^{-i \phi} + \langle \psi_0 | U^\dagger_{k-j}) H (e^{i\phi} |\psi_0\rangle + U_{k-j} |\psi_0\rangle) = \frac{1}{4} (\langle H \rangle + e^{-i\phi} \langle H U_{k-j} \rangle + e^{i\phi} \langle U_{k-j}^\dagger H \rangle + \langle U_{k-j}^\dagger H U_{k-j} \rangle) $
 
$P_- = \frac{1}{4} (\langle \psi_0 | e^{-i \phi} - \langle \psi_0 | U^\dagger_{k-j}) H (e^{i\phi} |\psi_0\rangle - U_{k-j} |\psi_0\rangle) = \frac{1}{4} (\langle H \rangle - e^{-i\phi} \langle H U_{k-j} \rangle - e^{i\phi} \langle U_{k-j}^\dagger H \rangle + \langle U_{k-j}^\dagger H U_{k-j} \rangle)$

$P_+ - P_- =  \frac{1}{2} ( e^{-i\phi} \langle H U_{k-j} \rangle + e^{i\phi} \langle U_{k-j}^\dagger H \rangle ) =  \frac{1}{2} ( e^{-i\phi} \langle H U_{k-j} \rangle +  (e^{-i\phi} \langle H U_{k-j}\rangle )^\dagger) = Re[e^{-i\phi}\langle H U_{k-j} \rangle] $

6. For measurement in X basis (ancilia), we have to rewrite the ancilia wavefunction as $|+i\rangle$ and $|-i\rangle$
$|\psi \rangle = \frac{1}{\sqrt{2}} ( e^{i\phi}|0\rangle |\psi_0\rangle +U_{k-j} |1\rangle |\psi_0\rangle) = \frac{1}{2} (|+i\rangle (e^{i\phi} |\psi_0\rangle - i U_{k-j} |\psi_0\rangle) + |-i\rangle (e^{i\phi} |\psi_0\rangle +i U_{k-j} |\psi_0\rangle) $

$P_{+i} = \frac{1}{4} (\langle \psi_0 | e^{-i \phi} + i \langle \psi_0 | U^\dagger_{k-j}) H (e^{i\phi} |\psi_0\rangle - i U_{k-j} |\psi_0\rangle) = \frac{1}{4} (\langle H \rangle -i e^{-i\phi} \langle H U_{k-j} \rangle + i e^{i\phi} \langle U_{k-j}^\dagger H \rangle + \langle U_{k-j}^\dagger H U_{k-j} \rangle) $

$P_{-i} = \frac{1}{4} (\langle \psi_0 | e^{-i \phi}  - i \langle \psi_0 | U^\dagger_{k-j}) H (e^{i\phi} |\psi_0\rangle + i U_{k-j} |\psi_0\rangle) = \frac{1}{4} (\langle H \rangle + i e^{-i\phi} \langle H U_{k-j} \rangle - i e^{i\phi} \langle U_{k-j}^\dagger H \rangle + \langle U_{k-j}^\dagger H U_{k-j} \rangle) $

$P_{+i}-P_{-i} = \frac{1}{2} ( -i e^{-i \phi} \langle H U_{k-j} \rangle + i e^{i \phi} \langle U^\dagger_{k-j} H \rangle) =\frac{1}{2} ( -i e^{-i \phi} \langle H U_{k-j} \rangle + i (e^{-i \phi} \langle H U_{k-j} \rangle)^\dagger) = Im[e^{-i\phi}\langle H U_{k-j} \rangle] $

Notably, the implementation of the controlled gates can be combined with the measurement, i.e., we do not actually need to implement those gates in certain cases.

## Fermi-Hubbard Model  

$H = -t \sum_{\langle i, j \rangle, \sigma} (a^\dagger_{i, \sigma} a_{j, \sigma} + h.c.) + U \sum_i n_{i, \uparrow} n_{i, \downarrow} = \frac{-t}{2} \sum_{\langle i, j \rangle, \sigma} X_iX_j + Y_iY_j + \sum_i  \frac{U}{4} (1-Z_{i, \uparrow})(1-Z_{i, \downarrow})$

## Keep in mind ...

1. The effectiveness of quantum Krylov diagonalization (QKD) depends on the overlap between the initial state and the true ground state.
2. The controlled state-preparation gate on $|0\rangle$ can be modified to use other states; it does not need to be $|\psi_0\rangle$.
3. When implementing the final measured Pauli operators, it can be combined with the controlled gate. for example  $CX @ XZ @ CX = YY, CX @ YZ @C X = -XY, CX @ XI @ CX = XX, CX @ YI @ CX = YX, CX @ XX @ CX = XI, CX @ YX @ CX = YI$
