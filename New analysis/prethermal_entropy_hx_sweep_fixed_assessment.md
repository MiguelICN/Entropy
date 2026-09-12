# Updated assessment: corrected parallel, sector-batched entropy workflow

Date: 12 September 2026

Files compared:

- `prethermal_entropy_hx_sweep_parallel_sector_batched_TO_BE_FIXED.m`
- `prethermal_entropy_hx_sweep_unique_outputs_CORRECT_WORKFLOW.m`

Corrected implementation: `prethermal_entropy_hx_sweep_parallel_sector_batched_FIXED.m`.

## 1. Outcome

The corrected file preserves the fixed random product state, reflection-sector diagonalization, parallel sweep over hx, sector-sized reconstruction, batched phase evolution, and half-chain Schmidt entropy. It retains the linear/logarithmic grid selector and unique output directories. It restores the working reference's default end time, `tmax = 10^5`, and its full-window moving average of the `{time, entropy}` coordinates.

The central corrections concern Wolfram evaluation and worker initialization, rather than a change to the physical model. Result assignments now survive trailing semicolons inside timing blocks; numerical output is checked before indexing; the private package path is inserted explicitly into worker initialization; and the norm of the actual evolved batches is checked.

**Targeted Wolfram numerical tests passed, including 36 boundary cases and a 100,002-sample trajectory. Native-C execution and a genuine multi-kernel production run could not be verified in the remote evaluator.** The supplied private QMB package remains an external dependency, used unchanged by the delivered file.

## 2. What caused, or could cause, the Null errors?

### 2.1 What the attached source actually contains

The attached `TO_BE_FIXED` file already uses the syntactically valid form

```mathematica
{evolutionSeconds, entropy} = AbsoluteTiming[
    Switch[...]
];
```

Its `linearEntropyEvolutionSector` ends with the expression `entropy`, without an internal terminating semicolon. Its sweep timing likewise returns `ParallelMap[...]` directly.

I executed the original numerical definitions and original `runOneHx` serially on an independent four-spin Ising test model. They returned `"success" -> True` and a saved entropy vector of length 21. **Consequently, the attached text alone did not reproduce the reported Null trajectory failure.** A differing live definition, an edited timing cell, or stale worker definitions remains possible; none can be identified conclusively without the exact failing kernel state.

### 2.2 A verified Wolfram mechanism for producing Null

The following pattern is unsafe:

```mathematica
{seconds, values} = AbsoluteTiming[
    values = computation[];
];
```

The semicolon makes the timed expression return `Null`. The outer assignment then overwrites `values` with that `Null`, even though `computation[]` may have succeeded. A later `values[[indices]]` produces a `Part` error.

A Wolfram test with `Range[3]` reproduced this exactly. The corrected form is:

```mathematica
seconds = First[AbsoluteTiming[
    values = computation[];
]];
```

This keeps the result assigned inside the timing block. The delivered file consistently applies this arrangement to diagonalization, evolution, per-run export, and the overall sweep. The final semicolon is now safe in these blocks. [Wolfram AbsoluteTiming documentation](https://reference.wolfram.com/language/ref/AbsoluteTiming.html)

This is a preventive correction to an established failure mechanism, not a claim that a missing semicolon was found in the supplied text.

## 3. Specific corrections

| Area | Previous behavior | Corrected behavior |
|---|---|---|
| Timed results | Destructuring could overwrite results with Null if the timed body was edited to end in a semicolon | Assign the result inside `AbsoluteTiming`; take only the elapsed time outside |
| Entropy result validation | Indexed `entropy` without first checking its structure | Require a result association and a finite machine-real vector of exactly `nSamples` entries |
| Invalid grid dispatch | No explicit default in the per-run evolution selector | Returns a named failure for an invalid grid type |
| Worker package setup | `ParallelEvaluate[Get[qmbInitPath]]` before explicitly supplying that variable to workers | `With` inserts the actual path before evaluation on workers; package initialization is checked |
| Fixed-state distribution | `psi0` was a dependency of `runOneHx`, but absent from the explicit fixed-input list | Explicitly distribute `psi0` alongside its sector projections |
| Worker numerical readiness | No startup test of both batched evolution branches | Run a tiny linear and arbitrary-time smoke test on each worker before the sweep |
| No available workers | Could continue without the requested worker environment | Stop with a clear message if zero workers launch |
| Norm validation | Norm checks used independently reseeded states only | Also check the norm of every reconstructed production sample, reduced to one maximum per run |
| Reset boundary | Advanced once into a reset boundary, then discarded that value | Skip this redundant update; reseed normally at the next batch |
| Validation indices | Included the sample before the first reset but not necessarily the reset itself | Include both sides of that boundary, clipped to the trajectory length |
| Input validation | Positive batch/kernel/reset values could still be nonintegers | Require positive integers for these controls |
| Short logarithmic interval | Rounding could produce zero logarithmic intervals | Use at least one interval |
| Basis diagnostics | Orthogonality and projected norm errors were printed only | Enforce acceptance before the sweep |
| Export result | Relied only on file existence | Require a string export result and file existence before declaring a run successful |
| Compatibility metadata | The batched format omitted `nSteps` | Save it alongside `nSamples` |
| Default trajectory | Batched file used `tmax = 10^6` | Restore the working file's `tmax = 10^5` |
| Moving average | Causal startup average at the original time coordinates | Restore `MovingAverage[rawData, 200]` over complete windows, as in the working reference |

Explicit package loading matters because `ParallelEvaluate` does not automatically distribute all definitions needed by its argument. The `With` substitution avoids relying on the worker having a value for `qmbInitPath`. Definition distribution follows package loading. [Wolfram ParallelEvaluate documentation](https://reference.wolfram.com/language/ref/ParallelEvaluate.html)

The missing explicit `psi0` entry was not, by itself, proof of a bug: recursive definition distribution can supply dependencies. Listing it explicitly makes this essential input unambiguous.

## 4. Preserved computational workflow

The corrected file keeps the following order:

1. Load your QMB package and set the saved notebook's directory.
2. Define and compile reusable numerical routines.
3. Set parameters and construct the selected time grid.
4. Create a unique parameter/run output directory.
5. Build the reflection bases and generate one fixed RPS.
6. Define one locally scoped complete hx run.
7. Initialize workers, distribute definitions, and verify numerical readiness.
8. Run independent hx values with `ParallelMap`.
9. Save summaries and export plots.

Every Hamiltonian-dependent eigensystem, coefficient vector, and reconstruction matrix is created inside its own `runOneHx` call after constructing that Hamiltonian. No random state is regenerated inside the sweep.

The working reference was used to recover successful workflow conventions, not to discard the requested parallel and batched implementation.

### Mathematical content

For sector eigenvector matrices Wplus and Wminus, the state is still

\[
\psi(t)=B_+W_+\bigl(c_+e^{-iE_+t}\bigr)
       +B_-W_-\bigl(c_-e^{-iE_-t}\bigr).
\]

Both terms are added coherently in the computational basis before computing the entropy. Sector probabilities are not independently normalized, and sector entropies are not averaged.

For an equal cut, reshaping the complete state gives a matrix M with singular values s_j, and the observable remains

\[
S_A(t)=-\sum_{s_j^2>10^{-14}}s_j^2\log(s_j^2).
\]

The logarithm is natural, so the output is in nats. Machine precision, the probability cutoff, and the SVD method are unchanged. There is no Trotter approximation or state-space truncation.

## 5. Timing and evolution return values

The evolution helpers now return

```mathematica
<|
    "entropy" -> entropyVector,
    "maxNormError" -> maximumNormSquaredError
|>
```

This lets the actual batch norm travel with the entropy result. `runOneHx` validates the association and extracts the vector before indexing it. **The saved trajectory still uses the familiar `"entropy"` key containing only the vector.** Existing entropy-data consumers need not treat it as an association.

The new saved field `"maxEvolutionNormError"` is distinct from `"maxValidationNormError"`:

- `maxEvolutionNormError`: largest norm-squared error in the reconstructed production batches, including recurrence drift.
- `maxValidationNormError`: largest norm-squared error among the independent direct-phase validation states.

The latter alone cannot reveal recurrence drift. Both are useful and are kept separate.

The `runSummary` remains a list of result associations. Non-association or unsuccessful entries are caught before summary-key access. A failed numerical result therefore produces a diagnostic rather than a cascade of `Null[[...]]` messages.

## 6. Time grids and reset boundaries

### Linear grid

The retained grid is

\[
t_i=t_{\min}+i\,dt,\qquad i=0,\ldots,nSteps.
\]

Each batch is limited by the next reset boundary. When a boundary is reached, coefficients are regenerated from the original spectral coefficients at that absolute time. This preserves the scalar recurrence's reset schedule. A final partial batch is allowed and both endpoints are retained.

Fractional positive dt and nonzero tmin were explicitly tested. A one-sample trajectory (`nSteps=0`) is supported by the numerical routine. At the restored defaults there are 100,001 samples per field, or 2,100,021 samples across 21 fields.

### Logarithmic grid

The existing logarithmic selector remains. Its arbitrary-time helper computes every phase directly at the supplied time, in batches. It does not reuse a uniform-step recurrence. The exact requested grid is stored with logarithmic outputs.

Direct phases avoid recurrence accumulation, but do not eliminate the long-time effect of eigenvalue errors. Neither linear nor logarithmic mode constitutes a guarantee of arbitrary-time machine-precision accuracy.

## 7. Plotting conventions restored from the working file

The working file applies `MovingAverage` to the pairs `{t,S}`. The corrected file does the same. With window w:

\[
\bar t_i=\frac1w\sum_{j=i}^{i+w-1}t_j,
\qquad
\bar S_i=\frac1w\sum_{j=i}^{i+w-1}S_j.
\]

For N input samples there are N-w+1 complete-window averages. At the default N=100001 and w=200, there are 99,802 averaged points per curve. On a uniform grid the averaged coordinate is the window midpoint. The previous causal version retained N points, used short startup windows, and placed averages at the last sample's time; it therefore did not reproduce the working plot convention. [Wolfram MovingAverage documentation](https://reference.wolfram.com/language/ref/MovingAverage.html)

Raw data remain untouched. Raw and averaged plots remain point plots, with export-only behavior. A logarithmic horizontal axis omits nonpositive time points from that visualization. If fewer than 200 samples exist, the default moving-average export is skipped rather than attempting an oversized window.

On a logarithmic sampling grid, this remains a sample average of complete windows, not a physical-time-weighted integral. That distinction belongs to later interpretation; no new smoothing method has been introduced.

The existing large-plot guard is retained. At the restored defaults, the 2,100,021 raw samples are below its five-million-point threshold. Increasing tmax substantially may cause plot export to be skipped unless `forceHugePlotExport` is enabled. All numerical trajectories are still saved.

## 8. Wolfram verification

### Test environment and model

Tests ran through the selected Wolfram plugin. Native workers were unavailable: `LaunchKernels[2]` returned no kernels. Numerical routines were evaluated with `CompilationTarget -> "WVM"` in the test copy; the delivered script retains `CompilationTarget -> "C"`.

The independent test Hamiltonian was

\[
H=J\sum_{j=1}^{L-1}Z_jZ_{j+1}+\sum_{j=1}^L(h_zZ_j+h_xX_j),
\]

with open boundaries and Pauli operators. This is a test fixture, not a substitution in the delivered script or a claim about your QMB conventions. Product states were constructed explicitly from normalized random complex two-component site vectors using seed 42 for the targeted tests.

A two-field test also executed the setup, grid creation, unique run directory, manifest, per-field exports, and CSV/WXF summaries with `Map` substituted for `ParallelMap` in the test copy. Both runs succeeded, each saved 21 samples, and all five expected data/summary files existed. This verifies serial orchestration; it does not test inter-kernel transfer.

An attempted combined test including PNG rendering returned an internal error from the remote evaluator. PNG export was therefore not verified end to end in this environment. The data/summaries portion was rerun separately and succeeded.

### Observed results

| Test | Result |
|---|---|
| Original attached `runOneHx`, serial L=4 fixture | Success; 21 entropy samples exported |
| Timing with an internal trailing semicolon and outer result assignment | Reproduced `Null` |
| Corrected per-run linear workflow, L=4, tmin=0.37, dt=0.3, 20 steps | Success; 21 samples |
| Same linear run: selected entropy differences versus direct phases | Maximum 2.14e-16 nats |
| Same linear run: actual evolved-batch norm error | Maximum 2.22e-15 |
| Corrected logarithmic/arbitrary-time branch at six irregular times | Success; six samples exported |
| Same arbitrary-time run: selected direct-phase entropy difference | Maximum 4.44e-16 nats |
| Guard applied to `Null` | Rejected |
| Guard applied to a vector containing `Indeterminate` | Rejected |
| Boundary suite: nSteps in {0,1,13}, batch sizes {1,2,7,64}, resets {1,5,100000} | All 36 cases returned the expected sample count |
| Boundary suite versus independent `MatrixExp` state evolution | Maximum entropy difference 6.88e-15 nats |
| Boundary suite actual norm error | Maximum 1.11e-15 |
| L=4 trajectory with 100001 steps, batch=64, reset=100000 | 100,002 samples returned |
| Long trajectory: samples immediately around the reset and endpoint | Maximum entropy difference from direct spectral phases 2.62e-12 nats |
| Long trajectory actual norm error | Maximum 1.01e-12 |

The long-time comparison uses the same computed eigensystem for both routes. It checks batching, recurrence, and reset consistency, not eigenvalue error against an exact high-precision solution.

These numerical differences are comfortably below the delivered entropy validation threshold of 1e-8 nats and norm threshold of 1e-10 for the tested fixtures. They are not a universal bound for larger systems or later times.

## 9. Updated performance assessment

Parallel independent hx runs remain the appropriate first strategy for this 21-field task. Every worker holds one Hamiltonian-dependent calculation, and each trajectory is saved independently. The worker returns a small summary rather than all trajectories to the master.

Sector reconstruction retains the algebraic saving identified in the original assessment. For D=2^L and reflection dimensions dplus,dminus, dense reconstruction uses dplus^2+dminus^2 entries instead of D^2. At L=8 these dimensions are 136 and 120, so the dense-entry ratio is about 0.502. Sparse lifting then adds the sector contributions in the full basis.

Batching uses matrix-matrix products for several times at once. It may improve matrix reuse and reduce interpreter overhead while keeping every SVD and every time sample. This is a plausible performance improvement, not a measured universal speedup.

The new batch norm check adds O(D times batchSize) work and temporary storage, compared with O(D^2 times batchSize) dense reconstruction. It is lower-order work but not free. The finite-vector check adds a scan of the returned entropy vector; it happens once per trajectory, with smaller checks per batch.

Startup checks are tiny and run once per worker. They verify basic numerical usability, not successful native-C compilation or optimal BLAS threading. A machine may return functioning compiled objects while falling back from the requested C target; inspect compiler messages on the actual installation.

Eight workers are retained as the requested setting, not asserted to be optimal. RAM, native linear-algebra threads, memory bandwidth, and licensed kernel availability determine the practical best count. A serially successful fixture is not evidence of an eight-worker speedup. No timing numbers from the remote fixture should be used to predict Windows production runtime.

## 10. Remaining limits

- Your private QMB implementation was not attached and has not been independently audited here. The delivered file calls it unchanged.
- Native-C compilation and inter-kernel definition/library transfer must work on your Mathematica installation. The corrected startup checks detect basic initialization failures before a long sweep.
- The original Null symptom was not reproduced from the attached source alone; the revision removes the demonstrated timing hazard and prevents malformed results from reaching index operations.
- Sampled eigensystem residuals remain sampled, not a proof that every eigenpair meets a global error bound.
- Long-time direct spectral checks share the same machine eigenvalues and cannot certify those eigenvalues' accuracy.
- The code still generates data for prethermalization analysis; it does not extract plateau lifetimes, fit scaling laws, or establish a mechanism from spectral weights alone.
- Unique run folders prevent ordinary sequential reruns from overwriting each other. Directory allocation is not designed as a synchronization mechanism for two independently launched notebook sessions starting at exactly the same time.
- Plotting loads trajectories together after the sweep; its memory demand is separate from the per-worker evolution footprint.

## 11. How to use the corrected file

1. Open it in your saved Mathematica notebook workflow and verify `qmbInitPath`.
2. Evaluate from the beginning, including definitions and worker initialization. Do not evaluate only the final `ParallelMap` cell in a session retaining older helper definitions.
3. Change parameters in the Parameters section before running. The supplied default is the working reference's L=8, J=hz=1, tmin=0, tmax=10^5, dt=1, with the same 21 hx values and seed 12345.
4. Inspect worker startup results and the per-run validation summaries. Successful runs save under a new `run_###` folder.
5. Use `"entropy"` for the raw trajectory, `"maxEvolutionNormError"` for actual production drift, and `"validationEntropyDifferences"` for the sampled direct-phase comparison.

The corrected workflow preserves the intended physics and optimized computational strategy while making Wolfram result handling and worker setup explicit and testable.
