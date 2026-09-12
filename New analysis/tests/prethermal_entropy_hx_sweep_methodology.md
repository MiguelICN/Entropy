# Methodology for Long-Time Entanglement-Entropy Sweeps in the Mixed-Field Ising Model

## 1. Purpose

This workflow is designed to study prethermalization through the time evolution of the half-chain von Neumann entanglement entropy of a **single fixed random product state (RPS)** under a family of mixed-field Ising Hamiltonians.

The main control parameter is the transverse field \(h_x\), sampled logarithmically through its exponent,

\[
h_x = 10^\alpha,
\qquad
\alpha=-2,-1.9,\ldots,-0.1,0,
\]

while the other Hamiltonian parameters are held fixed.

For the current implementation,

\[
L=8,\qquad J=1,\qquad h_z=1.
\]

The principal observable is

\[
S_A(t)
=
-\mathrm{Tr}\bigl[\rho_A(t)\log\rho_A(t)\bigr],
\]

with the subsystem \(A\) chosen as the left half of the chain,

\[
L_A=L/2.
\]

The calculation is intended to resolve a possible **two-stage relaxation**:

1. an initial growth of entanglement toward a quasi-stationary value;
2. a parametrically slower drift toward the final finite-size equilibrium regime.

The separation between these two stages is the main numerical signature of prethermalization that this workflow is intended to expose.

---

# 2. Overall computational strategy

For every value of \(h_x\), the calculation follows the sequence

\[
H(h_x)
\rightarrow
\{H_+(h_x),H_-(h_x)\}
\rightarrow
\text{separate diagonalizations}
\rightarrow
|\psi(t)\rangle
\rightarrow
\text{Schmidt singular values}
\rightarrow
S_A(t).
\]

The same initial RPS is used for every Hamiltonian in the sweep.

This is essential. If a new random initial state were generated for every \(h_x\), differences between entropy curves would contain both:

- changes caused by the Hamiltonian;
- sample-to-sample fluctuations of the initial state.

By fixing the initial state once, the parameter sweep isolates the dynamical effect of changing \(h_x\).

---

# 3. Fixed initial random product state

The initial state is generated only once:

```mathematica
SeedRandom[seed];

psi0 = Developer`ToPackedArray[
    N[RandomChainProductState[L]]
];
```

A fixed random seed is used so that the complete calculation is reproducible.

The state has the tensor-product form

\[
|\psi_0\rangle
=
\bigotimes_{j=1}^{L}
|\phi_j\rangle,
\]

and therefore its bipartite entanglement entropy is initially

\[
S_A(0)=0
\]

up to numerical roundoff.

This provides an immediate diagnostic of the complete evolution and entropy pipeline.

---

# 4. Reflection symmetry

The Hamiltonian is assumed to possess spatial reflection symmetry about the center of the chain.

The reflection operator acts as

\[
\mathcal R
|s_1s_2\ldots s_L\rangle
=
|s_Ls_{L-1}\ldots s_1\rangle,
\]

with

\[
\mathcal R^2=I.
\]

Therefore the Hilbert space decomposes into two reflection sectors,

\[
\mathcal H
=
\mathcal H_+
\oplus
\mathcal H_-,
\]

where

\[
\mathcal R|\psi_\pm\rangle
=
\pm|\psi_\pm\rangle.
\]

The code constructs orthonormal bases

\[
B_+,\qquad B_-,
\]

whose columns span the even and odd reflection sectors.

The projected Hamiltonians are

\[
H_+
=
B_+^\dagger H B_+,
\]

and

\[
H_-
=
B_-^\dagger H B_-.
\]

---

# 5. Why the symmetry sectors are diagonalized separately

The prethermal regime under study is expected to depend sensitively on spectral structure and degeneracies.

If the full Hamiltonian were diagonalized without resolving reflection symmetry, an eigensolver could return arbitrary linear combinations of eigenvectors belonging to different reflection sectors whenever their eigenvalues are exactly degenerate.

This does not change the physical spectrum, but it obscures the symmetry content of the eigenvectors and complicates the interpretation of:

- exact degeneracies;
- near-degeneracies;
- level crossings;
- symmetry-protected crossings;
- spectral statistics;
- the relation between spectral structure and relaxation times.

For this reason the workflow diagonalizes

\[
H_+
\]

and

\[
H_-
\]

independently.

The resulting spectra are retained separately as

```mathematica
evalsP
evalsM
```

for every value of \(h_x\).

This allows the later entropy analysis to be compared directly with the parity-resolved spectrum.

---

# 6. Sorted eigensystems

Each projected Hamiltonian is diagonalized as

```mathematica
{evals, evecs} = Eigensystem[Normal[M]];
```

and the eigenpairs are ordered by increasing energy.

Mathematica stores the eigenvectors returned by `Eigensystem` as rows.

Therefore, if

\[
V=
\begin{pmatrix}
v_1^T\\
v_2^T\\
\vdots
\end{pmatrix},
\]

the coefficients of a state in this eigenbasis are calculated through

\[
c_n
=
\langle v_n|\psi\rangle,
\]

implemented as

```mathematica
coeff = Conjugate[evecs] . psi;
```

and reconstruction uses

```mathematica
Transpose[evecs] . coeff
```

rather than `evecs . coeff`.

---

# 7. Projecting the fixed RPS into reflection sectors

Because the reflection bases depend only on \(L\), they are constructed once before the \(h_x\) sweep.

The fixed initial state is projected as

\[
|\psi_{0,+}\rangle
=
B_+^\dagger|\psi_0\rangle,
\]

\[
|\psi_{0,-}\rangle
=
B_-^\dagger|\psi_0\rangle.
\]

The corresponding code is

```mathematica
psi0P = ConjugateTranspose[Bp] . psi0;
psi0M = ConjugateTranspose[Bm] . psi0;
```

A normalization check is performed:

\[
\|\psi_{0,+}\|^2
+
\|\psi_{0,-}\|^2
=
1.
\]

These projected initial states are also independent of \(h_x\) and therefore need to be calculated only once.

---

# 8. Hamiltonian sweep

The parameter values are generated from

```mathematica
hxExponents = N[Range[-2, 0, 1/10]];
hxList = N[10^hxExponents];
```

giving

\[
h_x
=
10^{-2},
10^{-1.9},
10^{-1.8},
\ldots,
10^{-0.1},
10^0.
\]

There are 21 parameter values in total.

This sampling is uniform in

\[
\log_{10} h_x,
\]

not in \(h_x\) itself.

This is useful for prethermalization because relaxation times can change by orders of magnitude when the perturbation becomes small.

---

# 9. Operations repeated for every \(h_x\)

For each value of the transverse field, the script performs the following operations.

## 9.1 Construct the Hamiltonian

\[
H(h_x)
=
\texttt{IsingHamiltonian}[h_x,h_z,J,L].
\]

Only \(h_x\) changes between runs.

---

## 9.2 Project into reflection sectors

The Hamiltonian is projected as

\[
H_\pm(h_x)
=
B_\pm^\dagger
H(h_x)
B_\pm.
\]

---

## 9.3 Diagonalize both sectors

The script calculates

\[
H_+v_n^{(+)}
=
E_n^{(+)}v_n^{(+)},
\]

and

\[
H_-v_n^{(-)}
=
E_n^{(-)}v_n^{(-)}.
\]

The diagonalization time is recorded separately because it may become relevant when larger system sizes are studied.

---

## 9.4 Expand the same RPS in the new eigenbases

For each Hamiltonian, the expansion coefficients change because its eigenvectors change:

\[
c_n^{(+)}
=
\langle v_n^{(+)}|\psi_{0,+}\rangle,
\]

\[
c_n^{(-)}
=
\langle v_n^{(-)}|\psi_{0,-}\rangle.
\]

These are implemented through

```mathematica
coeffP = Conjugate[evecsP] . psi0P;
coeffM = Conjugate[evecsM] . psi0M;
```

---

# 10. Reconstructing the complete eigenbasis

The parity-resolved eigenvectors are lifted back into the complete computational Hilbert space:

\[
V_+
=
B_+ V_+^{(\mathrm{sector})},
\]

\[
V_-
=
B_- V_-^{(\mathrm{sector})}.
\]

They are then concatenated into

\[
V_{\mathrm{full}}
=
\left(
V_+\;V_-
\right).
\]

Similarly,

\[
E_{\mathrm{full}}
=
E_+\cup E_-,
\]

and

\[
c_{\mathrm{full}}
=
c_+\cup c_-.
\]

This does **not** undo the symmetry-resolved diagonalization.

The eigenvectors were obtained independently within their correct symmetry sectors. Joining them afterward is only a computational convenience for evolving the full physical state.

---

# 11. Reconstruction diagnostic

Before any long-time evolution is attempted, the code checks

\[
\left\|
V_{\mathrm{full}}
c_{\mathrm{full}}
-
\psi_0
\right\|.
\]

This is stored as

```mathematica
reconstructionError
```

and should be near machine precision.

A typical acceptable scale is

\[
10^{-12}
\text{--}
10^{-15},
\]

depending on matrix size and numerical conditioning.

A large reconstruction error indicates a problem with:

- eigenvector orientation;
- sector projection;
- ordering;
- normalization;
- numerical precision.

The long-time calculation should not be trusted if this diagnostic fails.

---

# 12. Real-eigenbasis optimization

For the current mixed-field Ising Hamiltonian and reflection basis, the Hamiltonian matrices are real.

Therefore their eigenvectors can also be chosen real.

The code verifies

\[
\max_{ij}
|\operatorname{Im} V_{ij}|
\]

and aborts if it exceeds

\[
10^{-10}.
\]

If the test passes, the full eigenvector matrix is stored as a real packed matrix,

\[
V=\operatorname{Re}V_{\mathrm{full}}.
\]

This enables the complex matrix-vector multiplication

\[
\psi
=
V(a_R+i a_I)
\]

to be replaced by two real matrix-vector products,

\[
\psi_R=Va_R,
\]

\[
\psi_I=Va_I.
\]

This is generally faster than repeatedly multiplying a real matrix by a complex vector.

---

# 13. Packed numerical arrays

All numerically intensive arrays are converted to machine-precision packed arrays using

```mathematica
Developer`ToPackedArray
```

where possible.

The relevant objects include

- eigenvalues;
- eigenvectors;
- expansion coefficients;
- trigonometric phase arrays;
- entropy trajectories.

Packed arrays reduce interpreter overhead and allow numerical linear algebra routines to operate on contiguous machine data.

---

# 14. Time grid

The current optimized sweep uses a uniform linear grid

\[
t_n=t_{\min}+n\Delta t.
\]

The present parameters are

\[
t_{\min}=0,
\qquad
\Delta t=1,
\qquad
t_{\max}=10^5.
\]

The number of time steps is

\[
N_t
=
\frac{t_{\max}-t_{\min}}{\Delta t}.
\]

The script checks explicitly that the selected `tmin`, `tmax`, and `dt` define an integer number of steps.

The current implementation is intentionally optimized for a fixed \(\Delta t\).

A later production version can generalize the interface to arbitrary explicit time grids, including logarithmic grids, without changing the physical methodology.

---

# 15. Spectral time evolution

After diagonalization,

\[
|\psi(t)\rangle
=
\sum_n
c_n
e^{-iE_nt}
|E_n\rangle.
\]

A direct implementation would recalculate

\[
e^{-iE_nt}
\]

for every eigenvalue and every time point.

That is unnecessarily expensive for a uniform grid.

For

\[
t_{n+1}=t_n+\Delta t,
\]

the eigenbasis amplitudes obey

\[
a_j(t+\Delta t)
=
e^{-iE_j\Delta t}
a_j(t).
\]

Define

\[
q_j
=
e^{-iE_j\Delta t}
=
\cos(E_j\Delta t)
-i\sin(E_j\Delta t).
\]

The arrays

\[
\cos(E_j\Delta t)
\]

and

\[
\sin(E_j\Delta t)
\]

are computed once per Hamiltonian.

---

# 16. C-compiled phase recurrence

Write

\[
a_j=a_{R,j}+ia_{I,j}.
\]

Then

\[
a_{R,j}'
=
a_{R,j}\cos(E_j\Delta t)
+
a_{I,j}\sin(E_j\Delta t),
\]

and

\[
a_{I,j}'
=
a_{I,j}\cos(E_j\Delta t)
-
a_{R,j}\sin(E_j\Delta t).
\]

This recurrence is implemented in the C-compiled function

```mathematica
phaseStepC
```

with

```mathematica
CompilationTarget -> "C"
```

and

```mathematica
RuntimeOptions -> "Speed".
```

This removes repeated transcendental-function evaluation from the main time loop.

---

# 17. Periodic exact reseeding

A long recurrence over \(10^5\), \(10^6\), or \(10^7\) time steps accumulates floating-point roundoff.

Therefore the recurrence is periodically reset using the exact spectral expression.

At a reset time \(t_r\),

\[
a_j(t_r)
=
c_j e^{-iE_jt_r}.
\]

The corresponding real and imaginary components are

\[
a_{R,j}
=
c_{R,j}\cos(E_jt_r)
+
c_{I,j}\sin(E_jt_r),
\]

\[
a_{I,j}
=
c_{I,j}\cos(E_jt_r)
-
c_{R,j}\sin(E_jt_r).
\]

This operation is implemented by the C-compiled function

```mathematica
reseedC.
```

The current value is

```mathematica
resetEvery = 100000;
```

so the recurrence is refreshed after every \(10^5\) steps.

This combines the speed of recurrence with protection against uncontrolled long-time phase drift.

---

# 18. Reconstruction of the physical state

At every time step the eigenbasis amplitudes are converted back to computational-basis amplitudes.

Because \(V\) is real,

\[
\psi_R(t)
=
V a_R(t),
\]

\[
\psi_I(t)
=
V a_I(t).
\]

The full state is

\[
|\psi(t)\rangle
=
\psi_R(t)
+
i\psi_I(t).
\]

The matrix-vector products are left to Wolfram Language's optimized native numerical linear-algebra backend rather than being reimplemented inside `Compile`.

---

# 19. Entanglement entropy by Schmidt decomposition

For a pure state, the fastest exact route to the bipartite von Neumann entropy is to avoid constructing the reduced density matrix explicitly.

For the bipartition

\[
A|B,
\]

reshape the state vector into

\[
M(t)
\in
\mathbb C^{2^{L_A}\times2^{L_B}}.
\]

For the half-chain cut,

\[
L_A=L_B=L/2.
\]

The singular-value decomposition is

\[
M
=
U\Sigma V^\dagger,
\]

with singular values

\[
s_\alpha.
\]

The eigenvalues of the reduced density matrix are

\[
p_\alpha=s_\alpha^2.
\]

Therefore,

\[
S_A(t)
=
-\sum_\alpha
p_\alpha\log p_\alpha.
\]

The code uses

```mathematica
SingularValueList[mat]
```

to obtain only the Schmidt singular values.

No reduced density matrix is explicitly formed.

---

# 20. Why the SVD itself is not C-compiled

`SingularValueList` is already executed by optimized native numerical linear-algebra routines.

Wrapping it inside `Compile` does not turn the SVD into a new C implementation.

Therefore the workflow leaves

```mathematica
SingularValueList
```

outside the compiled functions.

Only the inexpensive scalar entropy reduction

\[
-\sum_\alpha p_\alpha\log p_\alpha
\]

is compiled as

```mathematica
entropyFromSVC.
```

---

# 21. Entropy threshold for numerical zero

The compiled entropy reduction ignores Schmidt probabilities satisfying

\[
p_\alpha
\leq10^{-14}.
\]

This avoids evaluating

\[
p\log p
\]

for tiny values that are numerically indistinguishable from zero.

This threshold should be regarded as a numerical regularization parameter, not as a physical cutoff.

---

# 22. Memory-management strategy

Long trajectories can contain millions of entropy samples.

For example, a machine-real trajectory with

\[
10^7
\]

samples requires approximately

\[
8\times10^7\text{ bytes}
\approx80\text{ MB}
\]

before file-format overhead.

Keeping 21 such trajectories simultaneously in memory would be unnecessary.

The workflow therefore saves every \(h_x\) trajectory immediately after its calculation and then clears the large \(h_x\)-dependent arrays.

This allows the parameter sweep to proceed with approximately the memory footprint of a single Hamiltonian run.

---

# 23. Output structure

All results are stored in

```text
prethermal_entropy_hx_sweep/
```

relative to the notebook/script directory.

The directory contains

```text
manifest.wxf
run_01.wxf
run_02.wxf
...
run_21.wxf
run_summary.csv
run_summary.wxf
```

---

# 24. Manifest

The file

```text
manifest.wxf
```

contains the quantities shared by all runs:

- \(L\);
- \(J\);
- \(h_z\);
- random seed;
- the fixed initial state `psi0`;
- the list of \(h_x\) exponents;
- the list of \(h_x\) values;
- `tmin`;
- `tmax`;
- `dt`;
- number of time steps;
- reseeding interval;
- reflection-sector dimensions.

Saving the initial state itself is important for exact reproducibility.

---

# 25. Information saved for every \(h_x\)

Every file

```text
run_XX.wxf
```

contains an association with

- run index;
- \(L\);
- \(J\);
- \(h_z\);
- \(h_x\);
- exponent \(\alpha=\log_{10}h_x\);
- `tmin`;
- `tmax`;
- `dt`;
- number of time steps;
- full entropy trajectory;
- even-parity eigenvalues;
- odd-parity eigenvalues;
- initial-state reconstruction error;
- maximum imaginary component of the reconstructed eigenbasis;
- diagonalization time;
- entropy-evolution time;
- total runtime.

The parity-resolved eigenvalues are retained because the prethermalization mechanism may be directly related to degeneracies and near-degeneracies.

---

# 26. Timing diagnostics

For each value of \(h_x\), three timings are stored:

## Diagonalization time

Time required to diagonalize

\[
H_+
\]

and

\[
H_-.
\]

## Evolution time

Time required to compute the complete entropy trajectory after the eigensystem is already known.

## Total time

Wall-clock time for the complete \(h_x\) run.

These diagnostics will identify the dominant bottleneck as \(L\) and \(t_{\max}\) are increased.

For small \(L\) and extremely long trajectories, the entropy evolution is expected to dominate.

For larger \(L\), the symmetry-resolved diagonalizations and SVD calculations will become increasingly important.

---

# 27. Expected qualitative entropy evolution

The precise result must be established numerically.

However, the calculation is designed to test the following prethermalization scenario.

For sufficiently small \(h_x\), if the \(h_x=0\) Hamiltonian possesses the strong degeneracy structure responsible for the approximate dynamical constraint, then introducing a small \(h_x\) should weakly split or mix that structure.

One then expects two characteristic dynamical scales.

---

## 27.1 Early-time entanglement growth

Starting from

\[
S_A(0)=0,
\]

interactions generate entanglement and the entropy initially increases.

The early-time dynamics can remain relatively similar over a range of small \(h_x\), especially if the dominant local energy scales are controlled by \(J\) and \(h_z\).

---

## 27.2 Prethermal plateau or quasi-plateau

If a prethermal regime exists, the entropy should enter an intermediate-time window in which

\[
S_A(t)
\]

changes much more slowly than during its initial growth.

Schematically,

\[
S_A(t)
\approx
S_{\mathrm{pre}}
\]

for

\[
t_{\mathrm{local}}
\ll
t
\ll
\tau_{\mathrm{pre}}.
\]

The curve need not be perfectly flat.

For a finite system it may show:

- oscillations;
- finite-size recurrences;
- slow drift;
- beating between nearby frequencies.

The relevant signature is a parametrically reduced rate of relaxation over an extended time window.

---

## 27.3 Late-time entropy drift

At sufficiently long times, the perturbation \(h_x\) can destroy the approximate constraint responsible for the intermediate regime.

The entropy should then leave the prethermal plateau and approach its final finite-size stationary regime.

Thus the schematic behavior is

\[
0
\rightarrow
S_{\mathrm{pre}}
\rightarrow
S_{\mathrm{late}}.
\]

This is the expected two-step relaxation pattern.

---

# 28. Dependence on \(h_x\)

The primary quantity to compare across runs is the duration of the intermediate regime.

If \(h_x\) is the perturbation responsible for breaking the prethermal constraint, the generic expectation is

\[
h_x\downarrow
\quad\Longrightarrow\quad
\tau_{\mathrm{pre}}\uparrow.
\]

Therefore the curves should become progressively more separated in time as \(h_x\) decreases.

The smallest values,

\[
h_x\sim10^{-2},
\]

are the most likely to require very long simulations before the final departure from the prethermal regime can be observed.

The largest value,

\[
h_x=1,
\]

is expected to show much weaker separation between early and late relaxation scales if the prethermal behavior is indeed controlled by small \(h_x\).

This statement is a hypothesis to test, not a result assumed by the code.

---

# 29. Possible scaling of the prethermal lifetime

A later stage of the project should extract a characteristic prethermal lifetime

\[
\tau_{\mathrm{pre}}(h_x).
\]

Potential forms to test include a power law,

\[
\tau_{\mathrm{pre}}
\propto
h_x^{-\gamma},
\]

or a more rapid dependence such as

\[
\tau_{\mathrm{pre}}
\sim
\exp\left(\frac{c}{h_x^\nu}\right).
\]

The present sweep does not assume either scaling.

The logarithmic sampling of \(h_x\) is intentionally suitable for distinguishing such possibilities once a reliable definition of \(\tau_{\mathrm{pre}}\) has been established.

---

# 30. Entropy-based definition of the departure time

A useful later definition is to compare the perturbed entropy trajectory to an appropriate reference evolution.

For example, define a departure time as the earliest time satisfying

\[
|S_{h_x}(t)-S_{\mathrm{ref}}(t)|
>
\delta,
\]

for some fixed entropy tolerance \(\delta\).

Alternatively, one may define the end of the prethermal plateau from:

- the derivative \(dS/dt\);
- a smoothed entropy curve;
- deviation from an intermediate-time average;
- deviation from the \(h_x=0\) or weak-perturbation reference trajectory.

The final definition should be fixed before fitting scaling laws for

\[
\tau_{\mathrm{pre}}(h_x).
\]

---

# 31. What would constitute convincing evidence of prethermalization

A single entropy curve with a shoulder is not sufficient.

A stronger numerical case would show all of the following:

1. **Two distinguishable relaxation stages** in \(S_A(t)\).

2. **Systematic movement of the second relaxation scale** as \(h_x\) is varied.

3. A prethermal lifetime

   \[
   \tau_{\mathrm{pre}}(h_x)
   \]

   that increases as the perturbation becomes weaker.

4. A reasonably stable intermediate entropy value over the relevant time window.

5. Consistency of the behavior with the parity-resolved spectral evolution.

6. Persistence, at least qualitatively, when the system size \(L\) is varied.

7. Persistence across more than one RPS once the single-state methodology has been validated.

The present script addresses items 1--5 for one fixed initial state.

---

# 32. Role of degeneracies

A central motivation for retaining the parity-resolved spectrum is that the prethermal timescale may originate in the splitting of exact or near-exact degeneracies.

At \(h_x=0\), suppose a degenerate manifold contains states

\[
|n,a\rangle
\]

with the same energy

\[
E_n.
\]

A perturbation

\[
h_x V
\]

can split these energies,

\[
E_{n,a}(h_x)
=
E_n+\delta E_{n,a}(h_x).
\]

The associated long timescale is then controlled by small frequency differences,

\[
\Delta E
=
|E_{n,a}-E_{n,b}|.
\]

A characteristic dephasing time is

\[
t_{\mathrm{dephase}}
\sim
\frac{1}{\Delta E}.
\]

Therefore increasingly small splittings as \(h_x\rightarrow0\) can naturally produce increasingly long dynamical timescales.

The saved `evalsP` and `evalsM` allow this hypothesis to be examined directly rather than inferring it solely from entropy dynamics.

---

# 33. Why parity-resolved degeneracies must be interpreted carefully

Two distinct types of equality can occur:

### Degeneracy within one parity sector

\[
E_n^{(+)}
=
E_m^{(+)}
\]

or

\[
E_n^{(-)}
=
E_m^{(-)}.
\]

These represent genuine degeneracies inside the same symmetry sector.

### Crossing between parity sectors

\[
E_n^{(+)}
=
E_m^{(-)}.
\]

Since the two states belong to different symmetry sectors, such levels can cross without hybridizing.

These two cases have different physical implications.

The separate diagonalizations prevent them from being mixed together numerically.

---

# 34. Expected finite-size effects

The current default system size is

\[
L=8.
\]

This is useful for developing and validating the long-time methodology, but finite-size effects will be substantial.

Possible effects include:

- strong temporal oscillations;
- discrete-frequency beating;
- exact or near recurrences;
- deviations from smooth relaxation;
- plateau values that differ appreciably from thermodynamic expectations;
- late-time fluctuations that remain large.

Therefore one should not interpret every plateau-like feature as prethermalization.

The central evidence should be the **systematic parameter dependence of the timescale**, not merely the visual presence of a flat interval.

---

# 35. Page entropy as a reference scale

For a random pure state on a bipartite Hilbert space, the Page entropy gives a useful reference for the entanglement expected from a typical highly entangled state.

For a half-chain cut,

\[
d_A=d_B=2^{L/2}.
\]

The Page value is not itself a prediction for the prethermal plateau.

It is instead useful as a benchmark for judging how close the late-time state is to typical full-Hilbert-space entanglement.

A genuine prethermal plateau may remain below this scale before a later increase.

---

# 36. Quantities that should eventually be plotted

The most useful first comparison is

\[
S_A(t)
\]

for all \(h_x\) values on the same axes.

Because the total time range may span many decades, both of the following representations are useful.

## Linear-time plot

Best for resolving:

- initial entanglement growth;
- short-time oscillations;
- local relaxation.

## Logarithmic-time plot

Best for resolving:

- long prethermal plateaus;
- separation of relaxation times;
- parameter-dependent late-time departure.

Even though the present simulation uses a linear \(\Delta t=1\) grid, the stored trajectory can be plotted on a logarithmic horizontal axis afterward.

---

# 37. Recommended post-processing

The raw entropy should always be retained.

For visualization, one may additionally calculate a moving average or Gaussian-smoothed entropy,

\[
\bar S_A(t),
\]

to make slow trends visible beneath finite-size oscillations.

However:

- smoothing should never replace the raw data;
- the smoothing window must be much shorter than the suspected prethermal timescale;
- extracted transition times should be checked for stability against the smoothing window.

The saved WXF trajectories preserve the unsmoothed data for this purpose.

---

# 38. Future extension: arbitrary time grids

The present high-performance evolution uses a fixed uniform time step.

The physical spectral formula is nevertheless valid for arbitrary times,

\[
|\psi(t)\rangle
=
\sum_n
c_n e^{-iE_nt}|E_n\rangle.
\]

A future general routine can accept an explicit `timeGrid`.

For a uniform grid it can automatically use the fast recurrence.

For a logarithmic or irregular grid it can calculate

\[
c_ne^{-iE_nt}
\]

directly at each requested time.

Thus the present optimization does not constrain the conceptual methodology to linear-time sampling.

---

# 39. Future extension: coarse-grained parallelization

The current implementation intentionally optimizes the sequential core before introducing parallel execution.

The natural next parallelization strategy is to divide a long interval into large independent chunks.

For example,

\[
[0,10^7]
\]

can be decomposed as

\[
[0,T_c],
[T_c,2T_c],
\ldots.
\]

At the beginning of each chunk the amplitudes are calculated exactly using

\[
c_ne^{-iEt_0},
\]

and the fixed-step recurrence is then used internally.

This strategy has several advantages:

- very small parallel scheduling overhead;
- no need for kernels to communicate during the chunk;
- automatic suppression of long-term recurrence drift;
- straightforward restart of failed or interrupted chunks;
- compatibility with future cluster execution.

This is preferable to assigning individual time points to separate parallel kernels.

---

# 40. Future extension: several random product states

The current workflow deliberately uses one RPS.

Once the dynamical regime is identified, one can repeat the complete parameter sweep for several independent product states,

\[
|\psi_0^{(r)}\rangle.
\]

One can then study

\[
\langle S_A(t)\rangle_{\mathrm{RPS}}
\]

and its sample-to-sample variance.

This will determine whether the prethermal regime is generic within the RPS ensemble or strongly dependent on the particular initial condition.

The single-RPS calculation should be validated first because it is considerably easier to diagnose and optimize.

---

# 41. Immediate expected outcome of the current sweep

With the current parameter choice,

\[
L=8,\qquad
J=1,\qquad
h_z=1,
\]

and

\[
10^{-2}\leq h_x\leq1,
\]

the primary expected result is a family of 21 entropy trajectories,

\[
S_A(t;h_x).
\]

The main question is not simply whether their final entropy differs.

The important question is whether decreasing \(h_x\) produces an increasingly long interval during which the entropy remains close to an intermediate dynamical value before its late-time evolution becomes visible.

If this occurs systematically, the data should show a hierarchy

\[
\tau_{\mathrm{pre}}(10^{-2})
>
\tau_{\mathrm{pre}}(10^{-1.9})
>
\cdots
>
\tau_{\mathrm{pre}}(1),
\]

at least over the range where a well-defined prethermal regime exists.

The exact ordering can be affected by finite-size resonances and therefore should be verified rather than imposed.

---

# 42. Main diagnostics to inspect after the sweep

For each \(h_x\), verify:

### Initial entropy

\[
S_A(0)\approx0.
\]

### Reconstruction error

\[
\|V_{\mathrm{full}}c_{\mathrm{full}}-\psi_0\|
\ll1.
\]

### Norm conservation

Although not currently saved at every time step, spot checks should satisfy

\[
\langle\psi(t)|\psi(t)\rangle
\approx1.
\]

### Reflection-resolved spectrum

Inspect

\[
E_n^{(+)}(h_x)
\]

and

\[
E_n^{(-)}(h_x)
\]

for degeneracy splitting.

### Intermediate-time entropy

Search for a slowly evolving window.

### Departure time

Determine whether the late-time departure shifts systematically with \(h_x\).

### Finite-size oscillations

Distinguish periodic or quasiperiodic recurrences from a genuine slow relaxation envelope.

---

# 43. Summary of the methodology

The complete numerical methodology can be summarized as

\[
\boxed{
\begin{aligned}
&\text{Generate one fixed RPS}\\
&\downarrow\\
&\text{Construct reflection bases}\\
&\downarrow\\
&\text{For every }h_x:\\
&\qquad H(h_x)\\
&\qquad\downarrow\\
&\qquad H_+\oplus H_-\\
&\qquad\downarrow\\
&\qquad\text{separate exact diagonalization}\\
&\qquad\downarrow\\
&\qquad\text{expand the same RPS}\\
&\qquad\downarrow\\
&\qquad\text{packed spectral representation}\\
&\qquad\downarrow\\
&\qquad\text{C-compiled phase recurrence}\\
&\qquad\downarrow\\
&\qquad|\psi(t)\rangle\\
&\qquad\downarrow\\
&\qquad\text{half-chain SVD}\\
&\qquad\downarrow\\
&\qquad S_A(t)\\
&\qquad\downarrow\\
&\qquad\text{save trajectory + parity spectra}
\end{aligned}
}
\]

The methodological priorities are:

1. preserve the reflection symmetry sectors exactly;
2. use the same initial state for the complete parameter sweep;
3. retain parity-resolved spectra because degeneracy structure is central;
4. calculate entanglement directly from Schmidt singular values;
5. remove unnecessary transcendental evaluations from the long-time loop;
6. periodically reseed phases to control numerical drift;
7. save each long trajectory immediately to control memory usage;
8. establish the optimized sequential kernel before parallelizing;
9. use the \(h_x\)-dependence of the entropy relaxation timescale as the primary diagnostic of prethermalization.

The present script should therefore be regarded as the **single-state, symmetry-resolved, long-time baseline calculation** from which the later parallel, multi-state, larger-\(L\), and quantitative lifetime-scaling analyses will be built.
