## Detector Metadata

The `detray::detectors` component provides the detector metadata headers
(`detray/detectors/<name>_metadata.hpp`). Only hand-written headers, currently
`telescope_metadata.hpp`, live in `include/`. All others are generated at build
time by the `detray_metadata` package in `codegen/`; see
[`codegen/README.md`](codegen/README.md) for the vocabulary, layout and how to
add shipped or external configurations.
