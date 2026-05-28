# Agent Notes

## Task submission conventions for Li8 NO2B runs

Use this convention when preparing or submitting BIGSTICK Li8 NCSM jobs with NO2B `.bin` interactions.

- Place each run under:
  - `runs/Li8/<force>_emax<emax>_e3max<e3max>/hw<hw>_Nmax<Nmax>_beta<beta>/`
- Inside each run directory, use:
  - `input/` for BIGSTICK input files and Slurm scripts
  - `output/` as the working directory for the submitted BIGSTICK job
- Do not put `auto12` or `no2b` in submitted input, output, or log file names unless the user explicitly asks for it. Those are input-method details, not part of the preferred naming convention.
- Name BIGSTICK input files as:
  - `input/Li8_<force>_emax<emax>_e3max<e3max>_hw<hw>_Nmax<Nmax>_beta<beta>_<nstates>states.in`
- Use the same base name as the BIGSTICK output prefix, i.e. the second line of the `.in` file:
  - `Li8_<force>_emax<emax>_e3max<e3max>_hw<hw>_Nmax<Nmax>_beta<beta>_<nstates>states`
- BIGSTICK itself appends output suffixes such as `.res`, `.log`, `.wfn`, `.lcoef`, `.dres`, and `.occres` to that output prefix in the current working directory.
- Name Slurm scripts as either:
  - `input/run_Li8_<force>_hw<hw>_Nmax<Nmax>_beta<beta>_<nstates>states.slurm`
  - or a shorter equivalent when unambiguous, e.g. `input/run_li8_hw16_Nmax10_beta0.slurm`
- Set Slurm stdout/stderr logs under `output/`, for example:
  - `output/log_Li8_<force>_emax<emax>_e3max<e3max>_hw<hw>_Nmax<Nmax>_beta<beta>_<nstates>states_%j.out`
  - `%j` is the Slurm job id.
- Before submitting a Slurm job, check current node/partition availability and choose resources based on `Nmax`:
  - For `Nmax < 8`, choose a node with `256G` memory or less when available, and use about `32` OpenMP threads.
  - For `Nmax >= 8`, choose a node with `512G` memory or more, and use `64` or more OpenMP threads.
  - If needed for larger runs, request enough resources to occupy the whole selected node.
  - Do not submit immediately after choosing resources. Present the selected partition/node class, memory class, thread count, and Slurm resource settings to the user and wait for confirmation before `sbatch`.
- For the NO2B interaction, create a short symlink or copy inside `output/`, for example:
  - `output/li8_hw<hw>_no2b.bin`
  - pointing to `/lustre/home/2401110128/Forces/no2b/<force>/Li8/no2b_Li8_<force>_hw<hw>_emax<emax>_e3max<e3max>.bin`
- Run BIGSTICK from the `output/` directory so the code writes all generated files there and can open the short NO2B filename directly.

## NO2B binary interaction input

- Ported the NO2B binary reader from the sibling `BigstickPublick-3N` code into this repository.
- The implementation is contained in `src/binput_tbme.f90`.
- New interaction format code: `no2`.
- Intended input: the `no2b_*.bin` minipack binary written by the sibling `normal-order/src/no2b` program.
- The reader opens the file as unformatted stream input and reads:
  - header: oscillator frequency, `Emax`, orbit count, one-body count, two-body count
  - orbit table: `n`, `l`, `2j`, `2tz`
  - zero-body term
  - one-body matrix elements, folded into effective two-body pieces
  - normalized two-body matrix elements, placed into `ppme`, `pnme`, and `nnme`
- Usage in BIGSTICK interaction input:
  1. enter `no2`
  2. enter the `.bin` filename
  3. enter oscillator frequency and center-of-mass strength when prompted

## Mg24/USDB density matrix comparison notes

- Reference comparison setup used Mg24 with the USDB interaction.
- KSHELL was compiled from `/home/mengziyan/kshell`; the interactive KSHELL input script is in its `bin` directory.
- BIGSTICK and KSHELL spectra were first checked against each other before comparing densities.
- Documentation checked:
  - BIGSTICK manual: `docs/BIGSTICK_Manual.pdf`, especially Applications / Two-body densities.
  - KSHELL Brown density notes: `/home/mengziyan/kshell/Brown et al. - Lecture Notes in Nuclear Structure Physics_density.pdf`, especially the OBTD and TBTD formula sections.
- The BIGSTICK and KSHELL OBTD/TBTD definitions appear conventionally compatible after accounting for reduced-matrix-element normalization, pair normalization, and the time-reversal phase in the annihilation pair.
- Diagonal same-state OBTD for the lowest three Mg24 states agrees between BIGSTICK and KSHELL.
- Diagonal same-state TBTD:
  - An initial BIGSTICK run with OpenMP threading showed pp and nn agreement but pn disagreement.
  - A single-thread rerun with `OMP_NUM_THREADS=1` made pp, nn, and pn all agree with KSHELL at roughly `1e-7` to `1e-5`.
  - Confirmed cause: in `src/bdenslib4.f90`, the PN backward branch had an OpenMP reduction on `dmatpnhc` while the loop body updates `dmatpn`. This created a data race for threaded runs.
  - Source fix retained per user request: the PN backward branch now uses `reduction(+:dmatpn)`.
- Historical transition-density issue:
  - Before the retained full-density indexing changes, transition OBTD agreed but full transition TBTD disagreed.
  - The mismatch was concentrated in full transition TBTD elements with non-scalar rank, `Jab != Jcd`, and/or off-diagonal pair blocks.
  - Direct scalar diagonal-style entries were already consistent once the OpenMP PN race was removed.
  - With the current retained fixes, Mg24/USDB full TBTD agrees with KSHELL at about `8.1e-06` max absolute difference when the correct `M` projection is used and invalid `-999` placeholder values are excluded.

## M-projection handling for BIGSTICK density benchmarks

- BIGSTICK density output can contain `-999.0000`; the BIGSTICK manual says this means the Clebsch-Gordan coefficient used for reduction vanishes and the calculation should be rerun with a different `M`.
- The root PDF `打印邮件.pdf` is an email exchange with Calvin W. Johnson. It confirms the same point: when density matrix elements are missing because CG coefficients vanish, rerun BIGSTICK with `M=1`, i.e. `2M=2`, instead of `M=0`.
- In BIGSTICK `src/bdenslib5.f90`, the reduction coefficient is hard-wired to the same left/right projection and `q=0`:
  - `cleb4reduce = cleb(jtot*2,0,JJi,Jz,JJf,Jz)`
- Therefore BIGSTICK currently evaluates `<Jf M|T^K_0|Ji M>` and divides by the corresponding CG coefficient. If that coefficient is zero, the density element is not usable.
- KSHELL handles this more automatically in `src/transit.F90`: it checks for a vanishing reduction CG, constructs an `M+1` right-state with the angular-momentum raising operator, and recomputes the missing entries. The Mg24 KSHELL log prints `*** called mup ***`.
- BIGSTICK has no equivalent density-stage `m-up` workflow in the current code. The `.wfn` header stores a fixed `Jz`; the two-body density routine reads left and right states from one `.wfn` file and assumes the same `Jz` for both.
- Benchmark rule:
  - Run BIGSTICK with `2Jz=0` for states involving `J=0`.
  - Run BIGSTICK again with `2Jz=2` for nonzero-`J` states and use that run as a complete rerun for those states.
  - Match states between `2Jz=0`, `2Jz=2`, and KSHELL by energy and angular momentum, not by raw state index.
  - Do not compare or treat `-999`, `-499.5000171`, `-706.3996865`, or other scaled `-999` placeholder/invalid values as physical densities.
  - Do not merge `2Jz=2` entries into the old `2Jz=0` output for benchmarking; compare the appropriate complete run directly.
- Important interpretation:
  - A properly reduced density matrix element should be independent of the chosen `M` projection, as long as the left/right vectors are the same physical states, the CG coefficient used for reduction is nonzero, and the phase/reduced-matrix-element convention is consistent.
  - Therefore, after reduction, `2Jz=0` and `2Jz=2` BIGSTICK outputs should agree on all common non-`-999` entries for the same physical `2+` states.
  - The observed `M0/M1` agreement for diagonal rank-0 entries but disagreement for common nonzero-rank entries is evidence of a non-scalar TBTD reduction/coupling/convention problem, not a physical `M` dependence.
  - CG zeros explain missing/invalid `-999` placeholder entries only; they do not explain nonzero common entries that differ between `M` runs or disagree with KSHELL.
- For Mg24/USDB:
  - `0+ -> 2+` with rank 2 is safe in `M=0` because the relevant `(0 2 2; 0 0 0)` 3j coefficient is nonzero.
  - `2+ -> 2+` with odd tensor rank fails in `M=0` because `(2 K 2; 0 0 0)` vanishes for odd `K`; use the `2Jz=2` run for these entries.
  - With the current retained full-density fixes, the `2Jz=2` rerun fills the missing nonzero-`J` entries and agrees with KSHELL within the current benchmark tolerance.

## Latest TBTD testing notes

- Retained source fixes/changes:
  - `src/bdenslib4.f90`: PN backward density branch uses `reduction(+:dmatpn)`, fixing the threaded PN data race.
  - `src/bdenslib4.f90`, `src/bdenslib5.f90`, `src/bjumplib_master.f90`, `src/bjumplib_weld.f90`, and `src/bparallel_opbundles.f90`: current full-density code keeps the `dens2bflag` full non-Hermitian same-species density path rather than folding everything through the triangular Hermitian-conjugate indexing.
  - Commit `5670a1c` records this as a tentative fix; keep judging correctness only by direct agreement with KSHELL or another independent reference, not by whether a difference merely gets smaller.

### Mg24/USDB after the retained fix

- KSHELL reference was the existing `/home/mengziyan/kshell/run_mg24_usdb_density_compare/log_Mg24_usdb_all_density.txt`.
- BIGSTICK energy matching:
  - `2Jz=0`: state 1 `0+` at `-87.10445`, state 2 `2+` at `-85.60215`, state 3 `2+` at `-82.98830`.
  - `2Jz=2`: state 1 `2+` at `-85.60215`, state 2 `2+` at `-82.98830`; the `0+` state is absent, as expected.
- Preferred comparison script:
  - `runs/mg24_usdb_density_compare/my_independent_compare.py`
  - Route state pairs involving `0+` through the complete `2Jz=0` run.
  - Route `2+ <-> 2+` pairs through the complete `2Jz=2` run.
  - Do not merge individual `2Jz=2` entries into the old `2Jz=0` output; choose the complete run appropriate to the state pair.
- Current comparison result:
  - BIGSTICK routed entries: `2Jz=0` used `3704`, invalid `-999` placeholders skipped `0`; `2Jz=2` used `12824`, invalid `-999` placeholders skipped `0`.
  - Common BIGSTICK/KSHELL entries: `16528`.
  - Max absolute difference: about `8.1e-06`.
  - Mean absolute difference: about `3.84e-07`.
- Threaded check:
  - `OMP_NUM_THREADS=4` BIGSTICK rerun files: `mg24_usdb_m0_3_2full_omp4.den2b`, `mg24_usdb_m1_2_2full_omp4.den2b`.
  - `OMP_NUM_THREADS=4` and `OMP_NUM_THREADS=1` density outputs matched exactly for the compared entries.
  - `OMP_NUM_THREADS=4` versus KSHELL still has max absolute difference about `8.1e-06`.

### Mg24/IMSRG five-state threaded check

- Interaction:
  - IMSRG source `.snt`: `/home/mengziyan/shell-model-obs/temp/EM1.8_2.0_sd-shell_o18_hw16_emax4_e3max4.snt`.
  - KSHELL run copy: `/home/mengziyan/kshell/run_mg24_imsrg/H.snt`.
  - SHA256 matched for the two files: `68da04899df3915978ed8def8f79a1cc1f885fd02fb96d2b036c67e58aec7429`.
- KSHELL reference:
  - Run directory: `/home/mengziyan/kshell/run_mg24_imsrg`.
  - Spectrum log: `log_Mg24_m0_5states.txt`.
  - Density log: `log_Mg24_density_5states.txt`.
  - KSHELL transit log prints `*** called mup ***`, so KSHELL filled missing `M=0` reduction cases internally.
- BIGSTICK threaded runs stored under `runs/mg24_imsrg`:
  - `input.bigstick_mg24_m0_5_omp4` -> `mg24_imsrg_m0_5_omp4.*`, `2Jz=0`, `OMP_NUM_THREADS=4`.
  - `input.bigstick_mg24_m0_5_2full_omp4` -> `mg24_imsrg_m0_5_2full_omp4.den2b`.
  - `input.bigstick_mg24_m1_5_omp4` -> `mg24_imsrg_m1_5_omp4.*`, `2Jz=2`, `OMP_NUM_THREADS=4`.
  - `input.bigstick_mg24_m1_5_2full_omp4` -> `mg24_imsrg_m1_5_2full_omp4.den2b`.
- BIGSTICK/KSHELL energy matching:
  - `2Jz=0`: state 1 `0+` at `-74.32028`, state 2 `2+` at `-72.48098`, state 3 `2+` at `-69.63942`, state 4 `4+` at `-68.84279`, state 5 `3+` at `-68.33667`.
  - `2Jz=2`: state 1 `2+` at `-72.48098`, state 2 `2+` at `-69.63942`, state 3 `4+` at `-68.84279`, state 4 `3+` at `-68.33667`, state 5 `4+` at `-66.88278`.
- Projection choice for this five-state benchmark:
  - Use the `2Jz=0` BIGSTICK file only for state pairs involving the `0+` state.
  - Use the `2Jz=2` BIGSTICK file for all nonzero-`J` pairs among `2+`, `2+`, `4+`, and `3+`.
  - A CG check showed that `M=1` has no reduction CG zeros for the allowed ranks among `J=2,2,4,3`; `M=0` has many zeros for those nonzero-`J` pairs.
  - The complete `2Jz=0` full TBTD file still contains global `-999`/scaled `-999` invalid placeholders, but none are used in the routed comparison; the `2Jz=2` routed entries also contain none.
- Comparison script:
  - `runs/mg24_imsrg/compare_tbtd_imsrg_5_omp4.py`
  - It canonicalizes same-species pair ordering, routes `M=0`/`M=1` by state pair, and enumerates Lanczos eigenvector phase signs before comparing.
- Current comparison result:
  - Routed invalid `-999` placeholders skipped: `0` for `M=0`, `0` for `M=1`.
  - Common BIGSTICK/KSHELL entries: `63042`.
  - Max absolute difference: `1.54e-05`.
  - Mean absolute difference: `2.84e-07`.
  - KSHELL-only entries: `14`; largest absolute KSHELL-only value about `9.2e-06`.

## Current interpretation

- The retained density changes make Mg24/USDB and Mg24/IMSRG full TBTD agree with KSHELL at about `1e-5` absolute precision when the correct `M` projection is used and `-999` placeholder values are excluded.
- `-999` and scaled `-999` values are invalid placeholders caused by vanishing reduction CG coefficients. They are not physical densities and must not be included in comparisons.
- Use complete reruns at the appropriate `M` projection and match physical states by energy and angular momentum. Do not judge an edit by whether it merely reduces a discrepancy; judge it by direct agreement with KSHELL or another independent reference.

## Verification

- `git diff --check`
- `make -C make gfortran-openmp`

The default `make openmp` target was not used for verification because this environment does not have `ifort`.
