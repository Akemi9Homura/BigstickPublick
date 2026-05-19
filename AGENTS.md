# Agent Notes

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
- Transition density comparison:
  - Transition OBTD agrees.
  - Full transition TBTD still disagrees even in the single-thread BIGSTICK rerun.
  - The remaining mismatch is concentrated in full transition TBTD elements with non-scalar rank, `Jab != Jcd`, and/or off-diagonal pair blocks.
  - Direct scalar diagonal-style entries are consistent once the OpenMP PN race is removed.
  - Current suspect area is `src/bdenslib5.f90`, especially `couple_2bdensXX` and `couple_2bdensPN` handling of direct versus Hermitian-conjugate blocks, exchanged pair ordering, and associated phases for full transition TBTD.
  - The exact line-level fix for transition TBTD has not been confirmed yet.
  - Trial restores of the old full non-Hermitian XX indexing / density jump-folding code paths were tested and rejected because they broke already-correct rank-0 and `Jab == Jcd` diagonal checks.

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
  - Do not compare or treat `-999`, `-499.5000171`, `-706.3996865`, or other scaled sentinel values as physical densities.
  - Do not merge `2Jz=2` entries into the old `2Jz=0` output for benchmarking; compare the appropriate complete run directly.
- Important interpretation:
  - A properly reduced density matrix element should be independent of the chosen `M` projection, as long as the left/right vectors are the same physical states, the CG coefficient used for reduction is nonzero, and the phase/reduced-matrix-element convention is consistent.
  - Therefore, after reduction, `2Jz=0` and `2Jz=2` BIGSTICK outputs should agree on all common non-sentinel entries for the same physical `2+` states.
  - The observed `M0/M1` agreement for diagonal rank-0 entries but disagreement for common nonzero-rank entries is evidence of a non-scalar TBTD reduction/coupling/convention problem, not a physical `M` dependence.
  - CG zeros explain missing/sentinel entries only; they do not explain nonzero common entries that differ between `M` runs or disagree with KSHELL.
- For Mg24/USDB:
  - `0+ -> 2+` with rank 2 is safe in `M=0` because the relevant `(0 2 2; 0 0 0)` 3j coefficient is nonzero.
  - `2+ -> 2+` with odd tensor rank fails in `M=0` because `(2 K 2; 0 0 0)` vanishes for odd `K`; use the `2Jz=2` run for these entries.
  - However, the same-state nonzero-rank TBTD mismatch is not solved just by changing `M`: a `2Jz=2` rerun fills the missing odd-rank entries but those entries still disagree with KSHELL.

## Latest TBTD testing notes

- The PN OpenMP reduction fix is retained in `src/bdenslib4.f90`.
- For benchmark reproducibility, BIGSTICK density runs are still being made with `OMP_NUM_THREADS=1` unless threaded behavior is the thing being tested.
- New benchmark runs stored under `runs/mg24_usdb_density_compare`:
  - `input.bigstick_mg24_usdb_m0_3` and `mg24_usdb_m0.*`: `2Jz=0`, lowest three states.
  - `input.bigstick_mg24_usdb_m1_2` and `mg24_usdb_m1.*`: `2Jz=2`, lowest two states.
  - `input.bigstick_mg24_usdb_m0_3_2full` and `mg24_usdb_m0_3_2full.den2b`: full TBTD for the three `2Jz=0` states.
  - `input.bigstick_mg24_usdb_m1_2_2full` and `mg24_usdb_m1_2_2full.den2b`: full TBTD for the two `2Jz=2` states.
- Energy matching:
  - `2Jz=0`: state 1 `0+` at `-87.10445`, state 2 `2+` at `-85.60215`, state 3 `2+` at `-82.98830`.
  - `2Jz=2`: state 1 `2+` at `-85.60215`, state 2 `2+` at `-82.98830`; the `0+` state is absent, as expected.
- Comparison scripts added:
  - `runs/mg24_usdb_density_compare/compare_tbtd_mproj.py` maps and merges `2Jz=0`/`2Jz=2`; this was useful diagnostically but should not be used as the final benchmark procedure.
  - `runs/mg24_usdb_density_compare/compare_tbtd_m1_only.py` compares only the complete `2Jz=2` BIGSTICK rerun for the two `2+` states against KSHELL, mapping BIGSTICK states `{1: 2, 2: 3}`.
- Preferred current benchmark comparison for the two `2+` states:
  - Use `python3 compare_tbtd_m1_only.py`.
  - `2Jz=2` BIGSTICK entries after state mapping: `12824`, skipped sentinels: `0`.
  - Common entries with KSHELL: `9864`.
  - diagonal rank-0 remains good: `n=520`, max difference about `5.6e-6`.
  - diagonal rank>0 still disagrees: max difference about `0.3182`.
  - transition TBTD between the two `2+` states still disagrees: max difference about `0.6097`.
- Diagnostic combined BIGSTICK-vs-KSHELL TBTD result from `compare_tbtd_mproj.py`:
  - `2Jz=0` entries: `10892`, skipped sentinels: `5252`.
  - mapped `2Jz=2` entries: `12824`, skipped sentinels: `0`, filled into combined set: `5636`.
  - combined common entries with KSHELL: `12720`.
  - diagonal rank-0 remains good: `n=780`, max difference about `5.5e-6`.
  - filled odd-rank diagonal entries still disagree: max difference about `0.2507`.
  - overall same-state nonzero-rank TBTD still disagrees: max difference about `0.2614`.
  - transition TBTD still disagrees: max difference about `0.2575`.
  - This combined result is not the benchmark procedure; it is only evidence that the `M=0` sentinels are not the remaining mismatch.
- Interpretation of the new benchmark:
  - `M=0` CG zeros explain the sentinel values only.
  - They do not explain the remaining nonzero-rank TBTD mismatch.
  - The next debugging target remains BIGSTICK's coupling from raw `dmatpp/dmatnn/dmatpn` to reduced nonzero-rank TBTD, especially `src/bdenslib5.f90`.
- Experimental changes tested and reverted:
  - Restoring `dens2bflag` branches in `src/bigstick_main.f90`, `src/bjumplib_master.f90`, and `src/bparallel_opbundles.f90`.
  - Switching same-species density jumps/coupling in `src/bjumplib_weld.f90` and `src/bdenslib5.f90` from triangular Hermitian-style indexing to full non-Hermitian block indexing.
  - Disabling Hermitian-conjugate XX coupling in `src/bdenslib5.f90`.
  - These variants either doubled/broke rank-0 diagonal entries or damaged entries that already agreed with KSHELL, so they were not kept.
- Remaining known problem:
  - Full TBTD with nonzero tensor rank and especially `Jab != Jcd` still disagrees with KSHELL.
  - PN off-diagonal coupled blocks still show large differences; transition TBTD still has large differences.
- Next steps:
  - Compare BIGSTICK and KSHELL at the uncoupled m-scheme density level for a small set of failing keys before changing more source.
  - Start with representative same-state failures such as `(state 2 -> 2, a b Jab c d Jcd rank) = (2 2 2 2 2 3 3 2)` and PN failures such as `(3 5 2 2 5 3 2)`.
  - Trace those keys through `src/bdenslib4.f90` raw `dmatpp/dmatnn/dmatpn` accumulation and `src/bdenslib5.f90` coupling, then compare against KSHELL `operator_mscheme.f90:get_cpld_tbtd`.
  - Do not re-enable broad `dens2bflag` jump/opbundle branches unless the uncoupled comparison shows missing raw density blocks rather than a coupling/phase convention problem.

## Verification

- `git diff --check`
- `make -C make gfortran-openmp`

The default `make openmp` target was not used for verification because this environment does not have `ifort`.
