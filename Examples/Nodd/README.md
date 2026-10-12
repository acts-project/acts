# nODD examples

This directory contains ACTS adapters and scripts for the external
[nODD repository](https://github.com/asalzburger/nodd). Detector definitions,
materials, module placements and services remain owned by that repository.
The maintained entry point is its complete BeamPipe + Tracker compact.

`DD4hepNoddDetector` loads DD4hep for Geant4 and optionally constructs ACTS Gen3
geometry. `getNoddDetector()` selects the complete tracker with Gen3 enabled;
`gen3=False` supports DD4hep-only recording and explicit legacy pixel controls.
Sensor planes use native DD4hep transforms, bounds and XYZ axes. Every chip and
both paired strip faces remain separate. Root system IDs in the assembly are
pixel1–3, short strips4–5, long strips6–7; all other readout fields are retained.
A native Be slab at the wall mid-radius decorates its finite cylindrical portal;
native silicon slabs decorate sensors. Passive supports/services still need
a dedicated mapping campaign. No upstream ODD material map is applied.

Gen3 containers remain coaxial about the beam while sensor planes preserve their
actual tilts/stereo. A correctness-first try-all navigation policy avoids binning
assumptions for staggered/mixed modules. Performance binning is future work.
Navigation envelopes (default1mm) are software margins, not physical material.
`geometryReport()` returns JSON source IDs, converted centres/axes/bounds and
volume/portal counts for independent validation. DRAFT engineering/coverage
limitations remain.

## Build and configure

Build the nODD repository's DD4hep factory and maintained combined compact first,
using its documented compatible DD4hep/ROOT environment; build all default targets. The default inputs are:

- `<nodd>/build/dd4hep/detector/tracker/tracker.xml` and its XML includes.
- Factory libraries `nODDPixelBarrel`, `nODDBeamPipe`, `nODDShortStripBarrel`,
  `nODDShortStripEndcap`, `nODDLongStrip` and their `.components` registries.
- `<nodd>/build/dd4hep/detector/libnODDPixelBarrel.components`.

Enable `ACTS_BUILD_EXAMPLES_DD4HEP` and `ACTS_BUILD_PYTHON_BINDINGS` when building
ACTS. `ACTS_BUILD_EXAMPLES_GEANT4` and ROOT examples are also required for the
material recorder. The library and Python bindings follow the existing DD4hep
build option; a separate nODD dependency is not downloaded or linked at build
time. All new detector sources, binding source, Python helpers and scripts live
here; the top-level Examples and Python/Examples CMake files register them.

After activating the built ACTS runtime:

```sh
export NODD_PATH=/path/to/nodd
python - <<'PY'
from acts.examples.nodd import getNoddDetector
detector = getNoddDetector()
print(detector.config.xmlFileNames)
PY
```

`getNoddDetector(nodd_dir=..., build_dir=..., compact_file=...)` accepts explicit
overrides. Relative build directories are relative to the nODD checkout;
explicit compact paths are relative to the current directory. `NODD_PATH` is
used only when `nodd_dir` is omitted. The helper checks the input files, loads
the factory by absolute path and adds its component directory to
`DD4HEP_LIBRARY_PATH`, preserving existing entries. Keep the build and compact
compatible; changing the XML path does not rebuild the factory.

## Record tracker material

```sh
python Examples/Nodd/Scripts/nodd_material_recording.py \
  --nodd-dir "$NODD_PATH" \
  -n 10 -t 100 --seed 228 \
  --eta-range -4 4 --phi-range 0 360 \
  -o material/nodd_tracker
```

The installed command is also `nodd_material_recording.py`. Use `--build-dir`
for a nondefault nODD build, or `--compact` for an explicitly selected compact.
The output is `material/nodd_tracker.root`, tree `material_tracks`; the output
directory is created automatically. `--material-track-collection` changes both
the event-store collection and ROOT tree name. The output argument is a stem,
without the `.root` extension.

The pipeline mirrors `Examples/Scripts/Python/material_recording.py`: seed228,
point origin, uniform eta/phi generation, 1–10 GeV neutral geantinos, HepMC3
conversion, single-thread Geant4 material recording and ROOT material-track
writing with pre/post-step information and recalculated totals. Defaults are
1000 events and100 tracks per event. These are straight material probes, not
charged-track propagation, hit simulation or the luminous-region acceptance
study. Recording runs through the current external tracker model, including its
passive supports and effective services. No new geometry values are selected.

## Tests

From the ACTS checkout root in the installed ACTS Python environment:

```sh
pytest Examples/Nodd/Tests/test_nodd.py -v
NODD_PATH=/path/to/nodd pytest Examples/Nodd/Tests/test_nodd.py -v
```

Path/error tests run with the DD4hep nODD bindings. The external-model tests are
skipped when `NODD_PATH` is absent; with it they load the actual factory and run
eight material probes in a subprocess, checking the ROOT tree and nonzero finite
material totals. Geant4 is required for that recording test. Each subprocess
has its own DD4hep/Geant4 state.

## Validate Gen3 geometry and navigation

```sh
python Examples/Nodd/Scripts/nodd_geometry_validation.py \
  --nodd-dir "$NODD_PATH" --output validation/nodd_gen3 --tracks 128
```

This saves the actual surface/volume/pipe-portal report and ROOT propagation
steps and summaries for 128 seeded 10 GeV straight probes over eta [-4,4].
The companion nODD tools independently compare every sensitive surface with
the source inventory and check saved crossings. This is a software navigation
check; field propagation and detector acceptance require separate studies.
