# Rerunnable entropy sweeps and RPS ensemble averages

## What changed

The updated Mathematica file fixes output folders first and adds a separate RPS ensemble execution section. It keeps the supplied parameters: L=8, J=hz=1, hx exponents from -2 to 0 in increments of 0.2, and a linear grid from 0 to 10^6 with dt=1. The default execution mode remains `"Single"`.

### Folder failure

The original `tagNumber` allowed Mathematica's `*^` scientific notation into a Windows directory name. The `*` character is invalid there. Your screenshot shows such a path containing `tmax_1p*^6`. Directory creation failed, and the later worker exports consequently targeted a directory that did not exist.

Times, time steps, and seeds are now omitted from folder names. A typical layout is:

```text
<notebook directory>/
  entropy_sweeps/
    L8_J1_hz1_hxm2to0/
      run_<unique identifier>/
        description.txt
        metadata.wxf
        run_01.wxf
        run_02.wxf
        ...
        run_summary.wxf
        plots/
```

`hxm2to0` indicates the exponent range: hx=10^alpha with alpha from -2 to 0. The exact exponent list, including its spacing, is recorded in the metadata and description. Each execution gets a new run identifier, even if all parameters are unchanged. Shortened parameter tags are organizational labels; the saved metadata is authoritative.

Numeric tags allow only safe characters, including when J or hz uses scientific notation. Folder creation and an actual write test happen before diagonalization or evolution. In parallel mode, each worker also writes and removes its own probe file before expensive work starts.

The new layout shortens the generated path substantially. If the notebook itself lives under an exceptionally long or unwritable path, set `baseOutputDirectory` to a shorter existing directory. Ordinary filesystem permission failures are reported; the code cannot override them.

### Why plots were also missing

At the attached defaults there are 11 fields times 1,000,001 time samples: 11,000,011 points. The old five-million-point guard therefore skipped plotting automatically, independently of the folder failure. That guard is removed. All saved samples are supplied to plotting; no thinning is introduced. Rendering this many points can still take substantial time and memory.

## Rerun without exiting kernels

The script contains no `CloseKernels[]` call. It reuses the live pool, launches additional kernels only when fewer than requested exist, reloads the QMB package on workers, and redistributes the numerical definitions.

Every call to `runEntropySweep` builds a fresh configuration from the current parameters. The Hamiltonian parameters, dimensions, time grid, initial states, tolerances, and output directory travel to each worker in that configuration. A worker does not rely on an old global `L`, `dt`, `psi0`, or `outputDir`.

To use it:

1. Evaluate setup and definitions once in the saved notebook.
2. Edit/evaluate the parameter assignments you want to change.
3. Evaluate the relevant execution section again.

For example, after the definitions have been loaded:

```mathematica
L = 10;
J = 1.;
hz = 1.;
tmin = 0.;
tmax = 10^5;
dt = 1.;
hxExponents = N[Range[-2, 0, 1/10]];

singleResult = runEntropySweep["Single"];
```

`runEntropySweep` rebuilds derived dimensions, bases, time counts, and field values automatically. `batchSize` is a tunable setting; its existing value remains valid after a size change, and you can reevaluate its `Which[...]` assignment to restore the suggested size-dependent default.

A direct `runEntropySweep[...]` call saves numerical results. The complete execution sections additionally call the plotting function when `exportPlots=True`.

When evaluating the entire file, change its Parameters section first. Evaluating it from the beginning naturally reevaluates those assignments. An external one-line parameter change followed by evaluating an unchanged Parameters section would restore the values written in that section.

An already-running sweep uses its starting parameter snapshot. Change parameters and launch the next sweep after the current one completes. The code does not mutate an in-progress computation when globals change.

If more kernels are already open than `requestedKernels`, the existing pool remains open and is used; lowering this setting does not terminate kernels. `useParallel=False` provides a serial route for testing or installations without worker licenses.

## New complete RPS ensemble section

In Parameters, set:

```mathematica
calculationMode = "Ensemble";
ensembleSize = 20;
ensembleSeed = 12345;
```

Then evaluate the section titled:

```text
NEW: Ensemble of RPS, average first, save only mean trajectories, then plot
```

You can also call the underlying workflow directly:

```mathematica
ensembleResult = runEntropySweep["Ensemble"];
```

The ensemble consists of QMB `RandomChainProductState[L]` outputs generated using seeds

```mathematica
ensembleSeed + Range[0, ensembleSize - 1]
```

Each generation is wrapped in `BlockRandom`, so creating the ensemble does not advance the surrounding random-generator state. The same generated list is used for every hx. Changing `ensembleSeed` changes the realization; changing only hx or the time grid preserves the state list for fixed L, ensemble size, and QMB implementation. With `ensembleSize=1` and `ensembleSeed=seed`, the ensemble and single-state calculations coincide.

The library's distribution of random product states is preserved. No unverified assumption that it is a Haar distribution is added.

### Quantity averaged

For member r, the code first computes its pure-state half-chain entanglement entropy:

\[
S_A^{(r)}(t;h_x)=-\operatorname{Tr}\rho_A^{(r)}(t;h_x)\log\rho_A^{(r)}(t;h_x).
\]

It then saves

\[
\overline S_A(t;h_x)=\frac{1}{R}\sum_{r=1}^R S_A^{(r)}(t;h_x).
\]

This is the mean of individual entropies. It is not the entropy of an averaged density matrix or of an averaged state vector.

Every field's two reflection sectors are diagonalized once. The eigensystem is reused for all RPS members at that field. Sector coherence, batched reconstruction, the linear recurrence/reset schedule, direct logarithmic-grid phases, machine arithmetic, and the 10^-14 probability cutoff are retained.

The running mean is updated after each complete member trajectory:

```mathematica
meanEntropy += (entropy - meanEntropy)/member;
```

Only the mean and current member trajectory need to be retained, rather than a stack of R trajectories. The initial state vectors are held separately and are much shorter than long time trajectories. Rough memory scaling per worker is O(Ntime + R 2^L), in addition to eigensystems and reconstruction batches. Evolution time still grows approximately with the number of members; diagonalization is not repeated for each member.

Individual entropy trajectories are never exported. Each `run_XX.wxf` stores one final, unsmoothed averaged trajectory for that hx, along with its grid, field, state count, seeds, spectra, timings, and aggregate validation errors. A failed member fails the field run; the code does not silently average only successful members or export a partial mean as complete.

The mean is formed before plotting. No sample variance or standard-error trajectory is computed or saved in this version.

## Description and reproducibility

Each run has one human-readable `description.txt` containing:

- System size and Hamiltonian parameters, including the exact field/exponent lists.
- Single-state or ensemble mode, generator name, member count, and seeds.
- Linear-grid start/end/step and sample count, or the logarithmic-grid construction and zero-point setting.
- Batch size, reseed interval, entropy definition, and Wolfram/QMB identification information available to the workflow.

`metadata.wxf` additionally retains the initial state vectors for exact reuse without relying solely on seeds. These are initial vectors, not individual time trajectories. The exact logarithmic time grid is stored in the result files. All times remain reconstructible from the saved data.

## Replot later without recomputing

Use the saved run directory:

```mathematica
savedRunDirectory = "C:\\your\\saved\\run_directory";
replotEntropySweep[savedRunDirectory, 200];
```

Or, immediately after a successful ensemble run:

```mathematica
replotEntropySweep[ensembleResult["outputDir"], movingAverageWindow];
```

The function reads the saved `"entropy"` vectors and their metadata. It does not use the current global Hamiltonian parameters to reconstruct old coordinates or labels, and it does not generate states, diagonalize, or evolve.

It generates raw linear/log-time plots and moving-average linear/log-time plots. Here “raw” means the unsmoothed saved mean in ensemble mode. Curves use points without connecting lines. A blue-purple-red legend identifies the actual hx values.

The moving average is computed from the saved mean only after ensemble averaging. It starts with the original first sample, uses growing windows at startup, then uses the requested fixed window. The positive-time subset is constructed after the averaged data, fixing the stale/undefined-variable ordering in the attached plotting section. The linear plot includes t=0; the log-time plot necessarily excludes it. The saved data are never forced to zero or replaced by the smoothed values.

## Verification and limits

Wolfram tests completed for:

| Test | Result |
|---|---|
| Safe numeric folder tags, including 10^6 and 10^-12 | Valid tags and successful file writes |
| Repeating the same parameters | Distinct writable run folders |
| Single, ensemble, then changed L/J/start/end/dt in one session | All three data workflows succeeded |
| Three-member ensemble versus an independently formed mean using MatrixExp | Maximum difference about 4.94e-15 nats |
| One-member ensemble with the single-state seed | Exactly matched the single-state output in the test |
| Logarithmic ensemble branch including t=0 | Successful saved trajectory with the expected 10 samples |
| Changing unrelated live globals after building a configuration | The configuration's original size, fields, grid, and output folder were used |
| Deliberately invalid second member | Field failed; no partial mean file was created |
| Replot function on a saved fixture | Four valid Legended/Graphics structures were constructed |
| Startup average of {0,2,4,6}, window 3 | {0,1,2,4} |

The invalid-member test exposed a Wolfram control-flow pitfall: a `Return` inside the member loop did not exit the outer workflow as intended. The final worker implementation uses a tagged `Catch`/`Throw`, and the repeated failure test confirmed that incomplete means are not saved.

Tests used an explicit small open-chain Ising fixture and Wolfram's WVM compilation target. Your private QMB package was not available remotely, and native-C compilation and actual parallel worker reuse could not be exercised in that environment. No parallel speedup is claimed. The delivered file retains your QMB calls and C compilation target.

The remote evaluator returned an internal error when attempting actual PNG export. Graph construction was then verified separately with export intercepted, so the four plot structures are checked, but full PNG rendering on your installation remains unverified. Export failures are returned explicitly by the file rather than silently reported as success.

Long-time phase accuracy still depends on computed eigenvalues, not only reseeding. Ensemble averaging reduces sampling fluctuations; it does not repair numerical bias shared by the trajectories or independently establish prethermalization.
