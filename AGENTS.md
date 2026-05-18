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

## Verification

- `git diff --check`
- `make -C make gfortran-openmp`

The default `make openmp` target was not used for verification because this environment does not have `ifort`.
