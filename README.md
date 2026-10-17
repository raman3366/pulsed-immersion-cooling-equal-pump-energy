# Conjugate CFD Study of Pulsed Dielectric-Fluid Flow for Thermal Management of an 18650 Battery Module at Equal Pump Energy — data and code

Data and code for the conference paper of the same title (AISEE 2027). The
study simulates a six-cell 18650 lithium-ion module (2 × 3 cells, 20 mm
pitch) standing in HFE-7200 (Novec 7200) inside a closed housing, during a 3C
discharge of 900 s. The model is laminar, variable-property conjugate heat
transfer, solved with `chtMultiRegionFoam` of OpenFOAM v2606 (ESI) on a
structured multi-block mesh (0.57, 1.38 and 3.48 million cells on the
coarse, medium and fine levels). Steady coolant flow from 90 to
360 mL min⁻¹ is compared with two square waves (1 Hz, duty 0.25; 0.5 Hz,
duty 0.5) and a sine wave (2 Hz, amplitude 0.6). This package holds the 14
scored cases: each case's OpenFOAM setup and every monitor series it wrote,
the results table of the paper (Table I), the per-case numerical checks, the
derived temperature histories, and the scripts that regenerate Table I, the
checks and Figs. 1–3 from the case archives.

## Comparison basis

In laminar flow the hydraulic work needed to drive a flow Q(t) through a
fixed passage scales with ⟨Q²⟩, and ⟨Q²⟩ ≥ ⟨Q⟩² for every waveform, with
equality only for steady flow. **At the same mean flow a pulsed schedule
therefore always costs more pumping energy than steady flow.** Pulsation can
only be worth its cost if it removes heat better per unit of pumping energy,
so the primary comparison is at **matched pumping energy**: pulsed points are
read against the steady curve on fronts of pumping energy against the peak
cell temperature and against the cell-to-cell spread (Fig. 2). The runs at
the same mean flow (90 mL min⁻¹) are the mechanism study only (Fig. 3), not
an energy comparison.

The pumping energy is measured in each simulation, on the case's own
isothermal pre-circulation window (t = 2–30 s, whole periods), as the port
work less the change of the liquid's kinetic energy (and, for square waves,
less the time-discretisation loss at the valve edges), scaled to the 900 s
discharge. Once the cells heat up, buoyancy drives the liquid as well, so the
port work of the heated window is not a pump energy (paper, Section II,
"Pumping energy").

## Directory layout

```
README.md                  this file
LICENSE                    MIT (code: scripts/, geometry/)
LICENSE-DATA               CC BY 4.0 (data/ and documentation)
CITATION.cff, .zenodo.json citation metadata
requirements.txt           Python versions used
data/
  aisee_results.csv        Table I and the checks, one row per case (14 rows)
  gates/<key>.json         the scorer's numerical checks of each case
  timeseries/<key>.csv.gz  derived temperature histories of each case
  cases/<case>.tar.xz      each case's OpenFOAM setup + monitor output
  reference/               properties.yaml (cell and coolant constants),
                           heatgen.yaml (3C heat polynomial),
                           novec7200_properties.yaml (coolant property fits)
scripts/
  reproduce.sh             regenerate Table I, the checks and Figs. 1-3, and compare
  make_timeseries.py       writes data/timeseries/ from the case archives
  paper_numbers.py         derived numbers: Table I ratios, Table II, Section III
  check_reproduction.py    compares a regeneration with data/
  post/aisee_results.py    the scorer: gates, Table I quantities, figures
  post/aisee_figures.py    Figs. 1-3
  post/housing_metrics.py  the gates (energy ledger, mean flow, time step, ...)
  post/housing_io.py, v1_io.py      readers of the OpenFOAM monitor files and logs
  post/gates.py, pump_energy.py, precirc_check.py, campaign_stats.py, gci2.py
  model/housing_physics.py          case model (waveforms, heat model, regions)
  model/init_from.py                shared pre-circulated state (INIT_FROM; unused here)
  mesh/k0_blockmesh.py, k0_targets.py  structured mesh generator and mesh targets
  mesh/k0_measure.py                mesh statistics and the thin-head cell set
  hydrostatic_init.py               hydrostatic initial pressure
geometry/
  make_geometry.py         cell geometry manifest (STL) read by makePhysics.py
```

## Cases

All cases start from rest at 25 °C, ramp the inlet waveform up over 2 s,
circulate isothermally until t_pre = 30 s and then switch on the heat source
(Liu et al.'s 3C heat-generation polynomial, 1.73–2.22 W per cell). "Heated
time" is the simulated time after t_pre (t_end = 30 s + heated time). The
peak flow Q̂ is Q̄ for steady flow, Q̄(1 + A) for the sine and Q̄/D for the
square waves. The arm labelled C2 in the paper is `C2p` in the data.

| key | role | waveform | Q̄ (mL/min) | Q̂ (mL/min) | grid | Δt (ms) | heated time (s) | case directory (`data/cases/<dir>.tar.xz`) |
|---|---|---|---:|---:|---|---:|---:|---|
| steady_Q90 | curve | steady | 90 | 90 | medium | 12.5 | 900 | h_dtladder_steady_Q90_medium_dt12.5 |
| steady_Q135 | curve | steady | 135 | 135 | medium | 12.5 | 900 | h_aisee_steady_Q135_medium_dt12.5 |
| steady_Q180 | curve | steady | 180 | 180 | medium | 12.5 | 900 | h_aisee_steady_Q180_medium_dt12.5 |
| steady_Q360 | curve | steady | 360 | 360 | medium | 12.5 | 900 | h_aisee_steady_Q360_medium_dt12.5 |
| C1_Q90 | mean-flow | square 1 Hz, D 0.25 | 90 | 360 | medium | 12.5 | 900 | h_dtladder_sq1D0.25_Q90_medium_dt12.5 |
| C3_Q90 | mean-flow | sine 2 Hz, A 0.6 | 90 | 144 | medium | 12.5 | 900 | h_dtladder_sin2A0.6_Q90_medium_dt12.5 |
| C1_Q46.6 | matched-energy | square 1 Hz, D 0.25 | 46.6 | 186.4 | medium | 12.5 | 900 | h_aisee_sq1D0.25_Q46.6_medium_dt12.5 |
| C2p_Q64.2 | matched-energy | square 0.5 Hz, D 0.5 | 64.2 | 128.4 | medium | 12.5 | 900 | h_aisee_sq0.5D0.5_Q64.2_medium_dt12.5 |
| C3_Q82.9 | matched-energy | sine 2 Hz, A 0.6 | 82.9 | 132.6 | medium | 12.5 | 900 | h_aisee_sin2A0.6_Q82.9_medium_dt12.5 |
| steady_Q90_coarse_h270 | grid | steady | 90 | 90 | coarse | 12.5 | 270 | h_aisee_steady_Q90_coarse_dt12.5_h270 |
| steady_Q90_fine_h270 | grid | steady | 90 | 90 | fine | 12.5 | 270 | h_aisee_steady_Q90_fine_dt12.5_h270 |
| C1_Q90_dt6.25_h300 | dt | square 1 Hz, D 0.25 | 90 | 360 | medium | 6.25 | 300 | h_aisee_sq1D0.25_Q90_medium_dt6.25_h300 |
| steady_Q360_dt6.25 | dt | steady | 360 | 360 | medium | 6.25 | 900 | h_aisee_steady_Q360_medium_dt6.25 |
| steady_Q360_cand2 | scheme | steady | 360 | 360 | medium | 12.5 | 900 | h_aisee_steady_Q360_medium_dt12.5_cand2 |

Roles:
`curve` — the steady curve of Fig. 2;
`mean-flow` — pulsed at the mean flow of the steady reference, 90 mL/min
(mechanism study, Fig. 3; costs more energy than steady flow by construction);
`matched-energy` — pulsed at a reduced mean flow, 90/√(⟨Q²⟩/Q̄²) mL/min,
which a viscous-limit estimate puts at the pumping energy of steady 90 mL/min
(the comparison itself uses the measured energy);
`grid` — mesh check on the coarse and fine levels (first 270 s of heating);
`dt` — time-step check at Δt/2;
`scheme` — steady 360 mL/min with limited gradients (`cellLimited Gauss linear 1`,
the generator's scheme `cand2`) instead of `Gauss linear`.
Every case is a single realisation.

## Results table: `data/aisee_results.csv`

One row per case, written by `scripts/post/aisee_results.py` (Table I of the
paper is taken from it). Floats are printed with 6 significant digits.
Temperatures in °C are absolute; differences in K. "Heated samples" are the
monitor samples with t > t_pre.

| column | definition |
|---|---|
| `key` | case key (table above) |
| `role` | `curve`, `mean-flow`, `matched-energy`, `grid`, `dt` or `scheme` (above) |
| `arm` | `steady`, `C1` (square 1 Hz, D 0.25), `C2p` (square 0.5 Hz, D 0.5; C2 in the paper), `C3` (sine 2 Hz, A 0.6) |
| `waveform` | the arm in words |
| `Qbar_mL_min` | mean flow rate Q̄ (mL/min) |
| `grid` | mesh level: `coarse`, `medium`, `fine` |
| `dt_ms` | fixed time step (ms) |
| `source` | label of the computing batch that produced the run (`reused`: the three 90 mL/min runs taken over from the time-step study that preceded this paper, same model, mesh and protocol) |
| `t_heat_s` | heated time t_end − t_pre (s) |
| `Tmax_peak_C` | peak cell temperature: the largest value over the heated samples of max over the six cells of the cell maximum temperature (°C) |
| `Tmax_rise_K` | `Tmax_peak_C` above the inlet temperature, 25 °C (K) |
| `dTcell_end_K` | cell-to-cell spread at the end: maximum minus minimum of the six cell volume-mean temperatures at the last sample, t_end (K) |
| `dTcell_mean_K` | time mean of the spread over the heated samples (trapezoidal rule) (K) |
| `Tmodule_end_C` | module mean temperature at t_end: the six cell volume means weighted by the meshed cell volumes (°C) |
| `P_pump_uW` | pumping power: (E_port − ΔE_k − E_edge)/(t_b − t_a) on the isothermal window [t_a, t_b] (µW); E_port is the port work, ΔE_k the change of the liquid's kinetic energy, E_edge the valve-edge loss (square waves; 0 otherwise) |
| `E_pump_mJ` | pumping energy over the 900 s discharge, `P_pump_uW` × 900 s (mJ); the energy axis of Fig. 2 |
| `E_pump_raw_mJ` | the same before the valve-edge correction, (E_port − ΔE_k)/(t_b − t_a) × 900 s (mJ) |
| `E_edge_removed_mJ` | the valve-edge loss removed, E_edge/(t_b − t_a) × 900 s (mJ); E_edge = E_k,exch/n_e², with E_k,exch the total variation of the liquid's kinetic energy over the window and n_e the number of steps per 50 ms valve edge (4 at 12.5 ms, 8 at 6.25 ms); 0 for steady and sine flow |
| `pump_window_s` | the window [t_a, t_b] (s) on which the pumping energy is measured |
| `ledger_rel_pct` | residual of the enthalpy ledger of the heated stage, in the solver's discrete (BDF2) form, as a percentage of the injected energy |
| `gates_passed` | gates passed / gates evaluated |
| `gates_failed` | names of the failed gates, separated by `;` (empty: all passed) |

## Time series: `data/timeseries/<key>.csv.gz`

Written by `scripts/make_timeseries.py` with the scorer's own readers
(`aisee_results.series`, `housing_io.column`), one file per case,
comma-separated with a header row, gzip-compressed, one row per monitor
sample from the start of the run: every 0.05 s for steady flow, 0.0625 s for
the 2 Hz sine and the 1 Hz square wave, and 0.025 s for the 0.5 Hz square
wave. The heat source is on for t > 30 s.

| column | definition |
|---|---|
| `t_s` | simulation time (s) |
| `Tmax_cell0_C` … `Tmax_cell5_C` | largest temperature in cell i (°C) |
| `TvolAverage_cell0_C` … `TvolAverage_cell5_C` | volume-mean temperature of cell i (°C) |
| `Tmax_module_C` | max over the six cells of `Tmax_cell<i>_C` (°C); its maximum over t > 30 s is `Tmax_peak_C` |
| `dTcell_K` | max minus min of the six cell volume means (K); its last value is `dTcell_end_K` |
| `Tmodule_mean_C` | volume-weighted mean of the six cell volume means (°C); its last value is `Tmodule_end_C` |

Cells 0–2 and 3–5 are the two lines of three cells along the flow direction
x, upstream first.

## Numerical checks: `data/gates/<key>.json`

Every case is scored by `scripts/post/housing_metrics.py`, whose hard checks
("gates") a case must pass before it enters a figure; a case that fails one
is kept, flagged in `gates_failed` and marked in the paper (Table I, Fig. 2),
never plotted silently. Each JSON file holds `key`, `case_dir`, the list
`gates` (`name`, `ok`, `detail`) and the scorer's `info` block (energy-ledger
terms, pumping-energy measurement, Courant maxima and their locations, outlet
backflow, pre-circulation check and more). In the published files the
case path in `info.precirc.case` is reduced to the case directory name and
non-finite numbers are written as `null`.

| gate | check |
|---|---|
| `lattice` | the fixed Δt is commensurate with the waveform: integer steps per period, ON time and valve edge, and monitor sampling |
| `fixed_dt[log.chtMultiRegionFoam#k]` | every step of solver run k has the fixed Δt |
| `courant` | Courant number by region: claim-region bulk, the thin cells at the cap plane and in the headspace, and the outflow column at the riser entrance. All runs of this study carry the case generator's `DIAGNOSTIC` flag, and for such runs the Courant maxima are **reported, not gated** (the gate passes and its detail gives the maxima against the nominal limits; the paper reports them in Table II). The three `reused` runs predate the region monitors and report the global maximum |
| `max_principle` | the liquid stays between T_in − 0.05 K and the largest cell temperature at every heated sample |
| `mean_flow` | mean inlet flow over whole periods within 0.5 % of the intended Q̄ |
| `ledger` | enthalpy ledger of the heated stage closes within 2 % of the injected energy |
| `source_power` | the injected power of every cell equals the model power (volume normalisation) |
| `hydraulic_closure[pre]` | steady flow only: the mechanical-energy budget of the liquid on the isothermal window (port work against viscous dissipation, buoyancy work and kinetic-energy change) closes within a provisional 10 % |
| `closure_sign[pre]`, `closure_stationarity[pre]` | sine only: necessary checks of the same budget (sign of the residual by phase bin; stationarity over sub-windows) |
| `continuity[<monitor>:run]` | six key monitor series are continuous and consistent with the solver log (start and end times, sampling interval) |

The stationarity of the 30 s pre-circulation (`info.precirc`) is reported,
not gated: in 30 s it is not stationary, which the paper discusses under
Limitations.

### Flagged cases

Six cases fail one or two gates (column `gates_failed`) and are flagged in
Table I and Fig. 2:

| key | failed gate(s) | detail (from `data/gates/<key>.json`) |
|---|---|---|
| C1_Q90 | `max_principle` | liquid minimum 0.0512 K below T_in (bound 0.05 K), for 0.3 s (6 samples, t = 81.5–81.8 s) |
| steady_Q135 | `hydraulic_closure[pre]` | budget residual +10.99 % (provisional tolerance 10 %) |
| steady_Q180 | `hydraulic_closure[pre]` | residual +15.08 % |
| steady_Q360 | `max_principle`; `hydraulic_closure[pre]` | liquid minimum 0.245 K below T_in; residual +32.20 % |
| steady_Q360_dt6.25 | `max_principle`; `hydraulic_closure[pre]` | liquid minimum 0.245 K below T_in; residual +32.31 % |
| steady_Q360_cand2 | `hydraulic_closure[pre]` | residual +40.20 % |

The steady budget residual is measured against the cell-centred dissipation
estimate; it grows with the flow rate and falls with mesh refinement while
the port work changes by less than 0.6 % (paper, Section II, "Pumping energy").

## Case archives: `data/cases/<case>.tar.xz`

One archive per case, with the case directory name as the top folder (the
original layout is kept). Each holds:

- the case setup: `system/` (including the rendered `blockMeshDict`,
  `caseMeta.yaml` with every run parameter, `monitors` and `samplers`),
  `constant/` (dictionaries; no mesh), `0.orig/`, `Allrun`, `Allrun.pre`,
  `Allrun.heat`, `makeMesh.py`, `makePhysics.py`, `case.env` (the generator
  switches of the case), `k0_design.json` and `k0_mesh.json` (mesh design and
  measured mesh statistics) and `thin_head.json` (the thin-cell set of the
  Courant monitor);
- every monitor file the solver wrote, `postProcessing/<region>/<monitor>/<start time>/*.dat`
  (volume and surface integrals, extrema and fluxes of the fluid and the six
  cells, sampled every step or every few steps);
- `log.chtMultiRegionFoam`: a **slim extract** of the solver log, not the full
  log. It keeps, unchanged and in order, only the lines the scorer reads
  (`Time =`, the fluid region's `Courant Number` line, `deltaT =`,
  `Create time`, `Date`, `End`, and the heat sources'
  `- selected N cell(s) with volume V`), under a five-line comment header.
  It gives exactly the same gate results as the full log.

`steady_Q90` (`h_dtladder_steady_Q90_medium_dt12.5`) was archived with a
smaller set of setup files: `system/caseMeta.yaml`, `case.env` and the two
`k0_*.json` files, but not its dictionaries. `caseMeta.yaml` records every
parameter of the run; the same generator wrote the other two `h_dtladder_*`
cases, whose dictionaries are complete.

The setup files are the files the runs used, with comment-only edits: two
comment lines that referred to the authors' internal tooling were reworded
(`Allrun.pre`, `Allrun.heat`, `makeMesh.py`, `makePhysics.py`, `system/cell*/fvSolution`).
Comments in the case files also name scripts of the authors' working
repository that generated or drove the cases (for example
`scripts/mesh/make_k0housing.py`, `scripts/resume_case.sh`,
`generate_cases.py`) and internal decision labels (`D-…`, `NEW-…`, `WP …`).
Those scripts are not needed to run a case and are not part of this
package; every script that the case files call is (see "Re-running a case").

## Reproducing Table I, the checks and Figs. 1–3

Requirements: Linux or macOS with bash, GNU tar with xz support, and Python 3
with numpy, scipy, matplotlib, PyYAML (and pandas, see `requirements.txt`).
OpenFOAM is not needed.

```bash
git clone <this repository> && cd <repository>
bash scripts/reproduce.sh            # or: bash scripts/reproduce.sh /path/to/workdir
```

`reproduce.sh`

1. extracts every `data/cases/*.tar.xz` into `work/cases/` (about 3 GB);
2. runs the scorer on every case, recomputing every gate:

   ```bash
   AISEE_CASE_ROOT=work/cases AISEE_OUT=work/results \
   python3 scripts/post/aisee_results.py --refresh --allow-fail \
       C1_Q90 steady_Q360 steady_Q360_dt6.25 steady_Q135 steady_Q360_cand2 steady_Q180
   ```

   which writes `work/results/aisee_results.csv` (Table I), the gate records
   `work/results/aisee_cache/<case>.json` and Figs. 1–3
   `work/results/figures/fig_waveforms.pdf`, `fig_pareto.pdf` and
   `fig_histories.pdf`. The cases named in `--allow-fail` are the six flagged
   cases above; any other failing gate stops the scorer;
3. writes the time series to `work/results/timeseries/`;
4. prints the derived numbers of the paper with `scripts/paper_numbers.py`
   (also saved to `work/results/paper_numbers.txt`): the energy ratios r and
   r_QS of Table I, the gate summaries and the matched-time Δt/2,
   limited-gradient and three-mesh comparisons of Table II (with the grid
   convergence index of the spread), and the matched-energy and
   matched-mean-flow comparisons of Section III;
5. runs `scripts/check_reproduction.py work/results data`, which requires the
   regenerated table to equal `data/aisee_results.csv` cell by cell, every gate
   record to equal `data/gates/<key>.json`, every time series to equal
   `data/timeseries/<key>.csv.gz`, and the three figures to exist, and exits
   non-zero otherwise. (The scorer writes the CSV with CRLF line endings,
   Python's `csv` default; the published file has LF endings. The comparison
   is made cell by cell.)

The whole run takes a few minutes on a desktop machine. The environment
variables `AISEE_CASE_ROOT` (extracted case directories) and `AISEE_OUT`
(output directory) can also be set by hand to run the scorer,
`make_timeseries.py` or `paper_numbers.py` directly. The noise scale quoted
in the paper (0.004 K in the peak temperature, 0.4 mK in the end spread)
comes from a time-step pair and a seed pair of the steady 90 mL/min case run
before this study; those two runs are not part of this package, and the
values are constants in `scripts/post/aisee_results.py`.

## Re-running a case with OpenFOAM v2606

Each archive holds the case as it was run. The scripts in it do the following:

- **`Allrun`** (serial run): sources the OpenFOAM `RunFunctions`, runs
  `./Allrun.pre` and then the solver named in `system/controlDict`
  (`chtMultiRegionFoam`) in serial.
- **`Allrun.pre`** (pre-processing, no solver; it deletes `log.*`, time
  directories, `processor*`, `postProcessing/` and the mesh, so run it on a
  fresh extraction and never to resume a run):
  1. finds the code repository (the first parent directory of the case that
     holds `data/reference/properties.yaml`, else `$PULSED_REPO`) and sources
     `case.env`;
  2. writes the cell geometry manifest into `constant/triSurface` with
     `geometry/make_geometry.py`;
  3. runs `makePhysics.py`, which re-renders every property-dependent
     dictionary from `scripts/model/housing_physics.py` and the YAML data of
     the repository: fluid and cell thermophysical properties, the heat
     source (`constant/cell*/fvOptions`), the inlet waveform
     (`system/fluid/changeDictionaryDict`), schemes, `controlDict`, monitors,
     samplers and `system/caseMeta.yaml`; then sets the time scheme
     (`backward`, BDF2);
  4. runs `makeMesh.py <GRID>`, which writes `system/blockMeshDict` with
     `scripts/mesh/k0_blockmesh.py`, then `blockMesh`;
  5. copies `0.orig` to `0`, runs `splitMeshRegions -cellZones -overwrite`
     (one fluid and six cell regions), `changeDictionary -region <r>` for every
     region, and removes `U` and `p_rgh` from the solid regions;
  6. initialises a hydrostatic pressure in the fluid (`scripts/hydrostatic_init.py`
     writes `system/fluid/setExprFieldsDict`; `setExprFields -region fluid` at
     `writePrecision 17`);
  7. runs `checkMesh -allRegions -constant` and records mesh statistics in
     `k0_mesh.json` (`scripts/mesh/k0_measure.py`);
  8. writes the cell set of the thin headspace cells (`k0_measure.py
     --thin-head`, also `thin_head.json`) and turns it into the fluid cellZone
     `thinHead` with `topoSet -region fluid -dict system/topoSetDict.thinHead`;
     the Courant-number monitor in `system/monitors` reads that zone (this
     step is absent from the older `Allrun.pre` of the two reused
     `h_dtladder_*` cases, see below);
  9. with `INIT_FROM` set (not used by these cases) copies a shared
     pre-circulated state with `scripts/model/init_from.py` and runs
     `Allrun.heat`.
- **`Allrun.heat`**: the heat-on stage of the two-stage path (a separate
  isothermal pre-circulation run, or `INIT_FROM`): re-renders the physics with
  `RUN_STAGE=heat`, checks the stationarity of the pre-circulation
  (`scripts/post/precirc_check.py`) and optionally applies a temperature tilt
  seed. **None of the 14 cases used it**: all ran as one solver run from
  t = 0 to t_end (`run_stage: full` in `caseMeta.yaml`), with the heat source
  switched on at t_pre = 30 s by the time window of `constant/cell*/fvOptions`.

The production runs used `Allrun.pre`, then `decomposePar -allRegions`
(scotch, one MPI rank per physical core) and
`mpirun -np <N> chtMultiRegionFoam -parallel > log.chtMultiRegionFoam`.

### Running the setup from this repository

This repository holds every module and data file that `Allrun.pre`,
`makePhysics.py` and `makeMesh.py` read for these cases, at the paths the
scripts expect: `geometry/make_geometry.py`, `scripts/model/housing_physics.py`
(with `scripts/post/gates.py`), `scripts/mesh/k0_blockmesh.py`,
`scripts/mesh/k0_targets.py`, `scripts/hydrostatic_init.py`,
`scripts/mesh/k0_measure.py`, `scripts/model/init_from.py` (`INIT_FROM`
only), `scripts/post/precirc_check.py` (`Allrun.heat` only), and
`data/reference/properties.yaml`, `heatgen.yaml` and
`novec7200_properties.yaml`. The generator options these cases use are
covered (Liu's 3C heat polynomial with time in seconds, temperature-dependent
coolant properties); other C-rates, the depth-of-discharge time base and the
alternative heat model read data files that are not included.

A case extracted anywhere under the repository finds it on its own; anywhere
else, set `PULSED_REPO` to the repository root. The scripts call `python3`
(with numpy and PyYAML) unless `<repository>/.venv/bin/python` exists.

```bash
# from the repository root, in the OpenFOAM v2606 environment
mkdir -p work/rerun
tar -xJf data/cases/h_aisee_steady_Q90_coarse_dt12.5_h270.tar.xz -C work/rerun
cd work/rerun/h_aisee_steady_Q90_coarse_dt12.5_h270
# export PULSED_REPO=/path/to/repository   # only for a case outside the repository
./Allrun.pre          # geometry, physics, mesh, regions, initial fields, checkMesh, thinHead zone
# edit numberOfSubdomains in system/decomposeParDict (scotch), then:
decomposePar -allRegions
mpirun -np <N> chtMultiRegionFoam -parallel > log.chtMultiRegionFoam 2>&1
```

**What was tested for this package.** `Allrun.pre` was run on a fresh
extraction of the coarse grid-check case
(`h_aisee_steady_Q90_coarse_dt12.5_h270`), outside the repository with
`PULSED_REPO` set to a clean clone of it, using OpenFOAM v2606 (the same
version and build as the production runs, on an arm64 Linux machine). It
completed in about one minute: all seven regions passed `checkMesh`
(567,756 cells), the fluid cellZone `thinHead` was created (72,896 cells, as
in the archived `thin_head.json`), every re-rendered file in `system/`,
`constant/` and `0.orig/` is identical to the archived one, `k0_design.json`
is identical, and `k0_mesh.json` agrees with the archived one to round-off
(relative differences up to 1e-11; only the reported locations of tied
extrema differ). The decomposition and the solver were not run.

For every case, the Python steps of `Allrun.pre` (geometry, `makePhysics.py`,
`makeMesh.py`) were also run against this repository and the rendered files
compared with the archive:

- the other `h_aisee_*` cases: identical (after the time-scheme substitution
  of step 3, which normalises the whitespace of one line of
  `system/fluid/fvSchemes`), with one exception:
- `h_aisee_steady_Q360_medium_dt6.25` (`steady_Q360_dt6.25`): the archived
  `system/monitors` and `system/samplers` carry a 4-step monitor cadence and a
  sampler start of 1010 s, which match neither the case's own
  `caseMeta.yaml` (sampler interval 10 steps, t_end = 930 s) nor its monitor
  files (cell monitors every 0.05 s, i.e. 8 steps of 6.25 ms). The re-rendered
  files have the 8-step cadence and the 920 s start, consistent with both;
  every other file is identical. Treat the archived copies of these two files
  as not authoritative for this case;
- `h_dtladder_sq1D0.25_Q90_medium_dt12.5` and
  `h_dtladder_sin2A0.6_Q90_medium_dt12.5` (`C1_Q90`, `C3_Q90`, reused from an
  earlier study) were generated by an earlier version of the generator.
  Re-rendering them adds the later Courant-region monitors (`CoOut`,
  `CoBulkSleeves`, and a mask that reads the cellZone `thinHead`) and new
  `caseMeta.yaml` entries, and their `Allrun.pre` predates the thin-head step,
  so after their own `Allrun.pre` the solver would stop at the first Courant
  monitor. To re-run them as they were run, keep the archived dictionaries
  and run the OpenFOAM steps directly (below; their archived monitors do not
  read `thinHead`);
- `h_dtladder_steady_Q90_medium_dt12.5` (`steady_Q90`) has no archived
  dictionaries (see "Case archives"); its parameters are in `caseMeta.yaml`
  and `case.env`.

### Running from the archived dictionaries

The rendered dictionaries are in every archive except `steady_Q90`:
`system/blockMeshDict` (the mesh of the case's level), every dictionary of
`system/` and `constant/`, and `system/fluid/setExprFieldsDict` (the
hydrostatic initialisation). The OpenFOAM steps of `Allrun.pre` can be run on
them directly, without re-rendering (this was also run up to `checkMesh` on
the coarse case, with the same result):

```bash
cd <extracted case directory>
blockMesh
cp -r 0.orig 0
splitMeshRegions -cellZones -overwrite
for r in $(foamListRegions); do changeDictionary -region "$r"; done
for r in $(foamListRegions solid); do rm -f 0/"$r"/U 0/"$r"/p_rgh; done
# Allrun.pre runs the next step with writePrecision 17 in system/controlDict
setExprFields -region fluid
checkMesh -allRegions -constant
# h_aisee_* cases only: the fluid cellZone thinHead read by the Courant monitor
python3 "$PULSED_REPO"/scripts/mesh/k0_measure.py . --thin-head --level <GRID>
topoSet -region fluid -dict system/topoSetDict.thinHead
```

## Software environment

Simulations: OpenFOAM v2606 (ESI), Ubuntu package `openfoam2606`
2606.0~rc2-1 (build `_481094f-20260618`, `linux64GccDPInt32Opt`, GCC 13.3.0),
Open MPI 4.1.6, Ubuntu 24.04 LTS, on x86-64 servers with two 44-core Intel
Xeon (Sapphire Rapids) processors; the same platform for every case.

Scoring and figures (this package): Python 3.12.3 with numpy 1.26.4, pandas
2.1.4, matplotlib 3.6.3, scipy 1.11.4 and PyYAML 6.0.1 (the Ubuntu 24.04
system packages `python3-numpy`, `python3-pandas`, `python3-matplotlib`,
`python3-scipy`, `python3-yaml`; `requirements.txt` pins the same versions
for pip). Other recent versions are expected to work but were not tested;
the scorer uses `numpy.trapz`, so NumPy 1.x is the tested series.

## Not included

- field data (time directories), meshes (`polyMesh`) and decomposed
  `processor*` directories;
- the full solver logs (about 75 MB compressed per case) and the logs of the
  pre-processing utilities;
- the sampled field planes (`postProcessing/phasePlanes`, phase-locked slices
  of T, U and p_rgh).

The sampled field planes and the full solver logs are kept by the authors and
are available from the corresponding author on request. Field data and meshes
were not retained; the meshes are regenerated by each case's `Allrun.pre`.

## License

- Code (`scripts/`, `geometry/`): MIT License, see `LICENSE`.
- Data (`data/`, including the case archives) and documentation: Creative
  Commons Attribution 4.0 International (CC BY 4.0), see `LICENSE-DATA`.

## How to cite

If you use this data or code, please cite the paper and this dataset:

A. R. Lakshmipathi, D. Balakrishnan, AUTHOR 3 TO BE ADDED and AUTHOR 4 TO BE
ADDED, *Conjugate CFD Study of Pulsed Dielectric-Fluid Flow for Thermal
Management of an 18650 Battery Module at Equal Pump Energy — data and code*,
Zenodo, 2026. DOI: to be added.

Machine-readable metadata: `CITATION.cff` and `.zenodo.json`.

## Contact

Corresponding author: Anantha Raman Lakshmipathi, Department of Mechanical
Engineering, Rajalakshmi Institute of Technology, Chennai, India
(raman3366@gmail.com, ORCID 0000-0003-3465-286X).

## Provenance

The simulation set-up, the post-processing and figure code and this
documentation were prepared with the assistance of an AI system (Claude,
Anthropic) under the authors' direction; the authors checked the results.
