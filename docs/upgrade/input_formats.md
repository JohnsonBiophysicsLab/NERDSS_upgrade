# Phase 4 input format modernization

This first Phase 4 slice adds a JSON schema for a future NERDSS input format
without changing the runtime parser or deprecating legacy inputs. Existing
`.inp` and `.mol` files remain the compatibility contract.

## Goals

- Define a PR-reviewable JSON shape for simulation parameters, boundaries,
  molecule references, reaction declarations, and observables.
- Preserve legacy reaction equations as strings. The BNGL-like reaction grammar
  is intentionally not reimplemented in this slice.
- Keep molecule templates compatible with existing `.mol` files by allowing a
  JSON molecule entry to reference a legacy molecule file or embed inspected
  molecule metadata.
- Provide a small standard-library inspection helper for validation samples and
  migration discussions.

## Files

- `schemas/nerdss-input.schema.json` is the draft JSON Schema 2020-12 contract.
- `tools/legacy_input_to_json.py` inspects a legacy `.inp` file and emits JSON
  that matches the schema's compatibility-oriented shape for common samples.

## JSON shape

The top-level document is versioned:

```json
{
  "format": "nerdss-input-json",
  "schemaVersion": "0.1.0",
  "source": {
    "kind": "legacy-inspection",
    "path": "sample_inputs/VALIDATE_SUITE/clock_model/clock_model.inp"
  },
  "compatibility": {
    "legacyInputPath": "sample_inputs/VALIDATE_SUITE/clock_model/clock_model.inp",
    "legacyPreserved": true
  },
  "simulation": {
    "parameters": {},
    "boundaries": {},
    "molecules": [],
    "reactions": []
  }
}
```

`simulation.parameters` and `simulation.boundaries` map directly to the current
`start parameters` and `start boundaries` blocks. Known keys such as `nItr`,
`timeStep`, `waterBox`, and `sphereR` are typed in the schema; additional keys
are allowed as compatibility values so legacy extension points are not lost.

`simulation.molecules` records the molecule names and copy numbers declared in
the `.inp` file. Each molecule can keep a `moleculeFile` reference to the legacy
`Name.mol` file. Inline `template` content is optional and intended for
inspection, documentation, or future native JSON inputs.

`simulation.reactions` keeps each reaction equation verbatim and stores the
following key/value lines in `parameters`. This avoids a partial reaction parser
rewrite while still making rates, geometry parameters, labels, and coupled
reaction metadata visible to tooling.

## Legacy compatibility design

Legacy compatibility is explicit:

- JSON-native readers should accept `moleculeFile` references and may delegate
  `.mol` parsing to the existing parser.
- The schema allows additional typed legacy values in parameter objects.
- `nan` and `null` from legacy angle arrays are represented as JSON `null`.
- Symbolic values such as `M_PI` are preserved as strings.
- The helper does not rewrite reaction equations, molecule state syntax, or
  advanced copy-number/state expressions.

This keeps the modernization path additive. A future Phase 4 slice can add a C++
JSON reader that translates the schema into the existing parser data structures
while retaining the legacy `.inp` path as a supported frontend.

## Helper usage

Inspect a sample and write JSON:

```sh
python3 tools/legacy_input_to_json.py \
  sample_inputs/VALIDATE_SUITE/clock_model/clock_model.inp \
  --include-mol-files \
  --output /tmp/clock_model.input.json
```

List discovered sections only:

```sh
python3 tools/legacy_input_to_json.py \
  sample_inputs/VALIDATE_SUITE/clock_model/clock_model.inp \
  --list
```

The helper is intentionally limited. It understands block boundaries, comments,
simple `key = value` lines, molecule declarations of the common
`Name : copies` form, and reaction records formed by an equation followed by
key/value parameters. Unsupported legacy details are preserved as raw strings.
