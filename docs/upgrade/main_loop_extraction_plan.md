# Phase 3 Main Loop Extraction Plan

This is the first safe slice for extracting the serial NERDSS simulation loop.
It intentionally changes no runtime code. The current loop interleaves reaction
ordering, RNG consumption, propagation, compaction, and restart/output writes in
one large block, so the first extraction should preserve call order byte-for-byte
before moving any computation.

## Current Serial Loop Boundary

The serial main loop lives in `EXEs/nerdss.cpp`.

- Restart pre-write for restart simulations: lines 777-795.
- Iteration loop header and per-step timer setup: lines 807-815.
- Population, zeroth-order, and unimolecular-state reactions: lines 824-865.
- Implicit-lipid dissociation: lines 867-883.
- Separation measurement and possible reaction collection: lines 885-953.
- Bimolecular/transmission reaction decision and execution: lines 959-1221.
- Overlap checking and propagation: lines 1223-1391.
- Per-step trajectory, PDB, and transition writes: lines 1397-1424.
- Empty complex compaction: lines 1427-1475.
- Empty molecule compaction: lines 1477-1575.
- Per-step cleanup of reweighting/crossing state: lines 1577-1595.
- Restart, bonded-complex, checkpoint, and time/status output: lines 1617-1755.
- Final timestep output block after the loop: lines 1758-1845 and following.

The MPI loop in `EXEs/nerdss_mpi.cpp` is not the Phase 3 extraction target, but
it is useful as a reference because several concerns have already been separated
there:

- `check_perform_zeroth_first_order_reactions(...)`: lines 536-540.
- `measure_separations_to_identify_possible_reactions(...)`: lines 561-566.
- `perform_bimolecular_reactions(...)`: lines 623-626 and 799-802.
- `check_overlap(...)`: lines 644-646 and 820-822.
- `remove_empty_slots(...)`: lines 978-981.
- `write_output(...)`: lines 1002-1007.

## Non-Negotiable Preservation Rules

- Do not change the order of calls that can consume RNG. In the serial loop this
  includes reaction checks, `determine_if_reaction_occurs(...)`, transmission
  decisions using `rand_gsl()`, association/dissociation handlers, overlap
  resampling, propagation, and restart RNG writes.
- Do not change molecule or complex iteration order, including the order of
  subcell members, neighboring cells, `moleculeList`, `complexList`,
  `Molecule::emptyMolList`, or `Complex::emptyComList`.
- Do not change when temporary cross-reaction vectors, reweighting vectors,
  `trajStatus`, `isDissociated`, `justBoundThisStep`, or `justUnboundThisStep`
  are cleared.
- Do not move restart writes across cleanup, compaction, or RNG-state writes.
  Restart files are part of the behavioral surface, not just diagnostics.
- Keep serial and MPI changes separate unless a later slice explicitly validates
  both execution modes.

## Proposed Extraction Order

1. Extract a `SimulationLoopState` or similarly named context object that only
   groups existing references needed by the loop. The first code slice should
   compile with the context created but still call the original inline code.
2. Extract output-only helpers whose call sites are already single-purpose:
   trajectory/PDB/transition writes, bonded-complex writes, checkpoint writes,
   and time/status writes. Validate restart files and normal output file names.
3. Extract per-step cleanup helpers for clearing molecule cross-reaction state
   and complex `ncross`/`trajStatus`, preserving the exact serial cleanup order.
4. Extract empty-slot compaction. Prefer reusing or adapting the existing MPI
   `remove_empty_slots(...)` shape only after serial parity is proven on cases
   with destroyed molecules and complexes.
5. Extract reaction-discovery and reaction-execution helpers only after the
   earlier helpers are validated. These are the highest RNG-risk seams.
6. Extract propagation/overlap helpers last. These touch resampling,
   `propCalled`, `trajStatus`, and complex propagation order.

## First Code Slice Candidate

The safest no-behavior code candidate is an output helper for the repeated
restart-write block:

- restart pre-write: `EXEs/nerdss.cpp` lines 783-794;
- periodic restart write: `EXEs/nerdss.cpp` lines 1621-1635;
- checkpoint restart write: `EXEs/nerdss.cpp` lines 1647-1655;
- final restart write: `EXEs/nerdss.cpp` lines 1762-1770.

This still needs careful validation because `write_rng_state()` and
`write_rng_state_simItr(...)` must remain at the same logical point relative to
`write_restart(...)`.

## Validation Checklist For Any Code Slice

- Build and run the existing serial smoke runner:
  `python3 tools/run_smoke_tests.py --artifact-dir /tmp/nerdss-main-loop-smoke`.
- Run a restart regression with `rngwrite = true` that compares a continuous run
  against a restart-resumed run from the same seed. Compare at least `restart.dat`,
  `DATA/trajectory.xyz`, `DATA/species_time.dat`, and final coordinates.
- Include a destroyed-molecule/destroyed-complex case before extracting
  compaction, because empty-slot ordering is observable through molecule and
  complex indices.
- Include an association-capable case before extracting reaction discovery or
  execution, because this validates `crossbase`, `mycrossint`, `crossrxn`,
  `probvec`, and RNG order.
- Include an overlap/resampling case before extracting propagation, because
  propagation can consume RNG and mutates `trajStatus`.
- For documentation-only slices, validate by confirming the diff contains no
  source or build-system changes.

## Follow-Up Checklist

- Add a small restart-regression runner next to `tools/run_smoke_tests.py`.
- Capture baseline artifacts for the smoke input before any loop code moves.
- Add one minimal input with a bimolecular association path and fixed seed.
- Add one minimal input with creation/destruction so compaction can be tested.
- Add one minimal input that forces restart writing with `rngwrite = true`.
- Only then begin the first no-behavior helper extraction.
