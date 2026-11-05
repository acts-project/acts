# Generate a configuration file for the smearing digitizer

## Optional geometric charge quantization

A geometric digitization entry may contain:

```json
"charge-quantization": {
  "enabled": true,
  "bits": 4,
  "offset": 0.01,
  "step": 0.011105528959000351
}
```

Missing configuration or `enabled: false` preserves native digitization.
The C++ configuration is a plain value, disabled by default. JSON output always
includes its settings, including disabled defaults when the input block was absent.
Enabled quantization requires `doMerge=true`, geometric indices and
`digital=false`. ACTS preserves the existing native smearing and suppression,
sums accepted contributions on identical channels, and quantizes once before
final measurement weighting. No noise or suppression is added. Neighboring
cell connectivity and cluster membership are unchanged.

For `M = 2**bits - 1`, the response is
`clamp(1 + floor((activation-offset)/step), 1, M)`, with boundary ties assigned
upwards. Bit depths 1–16 are supported. Zero is reserved for absent channels.
Invalid or non-representable binning is rejected at initialization. Changing
`bits` alone preserves `step`, increasing dynamic range; adjusting `step` can
instead hold a chosen activation range fixed.

The existing cell activation field stores codes on enabled geometry entries;
codes weight both positions and the existing uncertainty formulas. Field types
and writer schemas are unchanged; four-bit information is not packed storage.
Native values are available by running the disabled configuration with the same
inputs and seed. Consult the configuration for mixed native/code detectors.

The last code begins at `offset + (M-1)*step`; clipping beyond its nominal bin
starts at `offset + M*step`. Detector parameters belong in detector-specific
configuration, not in the generic response implementation.

## Smearing configuration generator

As a convenience, a simple helper script is provided to help producing JSON configuration files for smearing digitization.
The script is located in the source tree at `Examples/Algorithms/Digitization/scripts/smearing-config.py`.

Each volume configuration is one logical block consisting of the following options:

* `--digi-smear-volume=X`: Specifies the numeric identifier of the
  volume configured in the present block. This option also acts as a
  delimiter for different blocks.
* `--digi-smear-indices=X[:Y...]`: Specifies the dimensions to be configured in a comma-delimited list
* `--digi-smear-types=X[:Y...]`: Specifies the smearer type for each dimension. The following values are allowed:
  * 0: Gaussian
  * 1: Truncated gaussian
  * 2: Clipped gaussian
  * 3: Uniform
  * 4: Digital
* `--digi-smear-parameters=X[:Y...]`: Configuration parameters, 1 for Gaussian, 3 for all others (1 parameter, 2 range values)

Such blocks may be repeated as often as needed to configure many volumes, but each block has to begin with `--digi-smear-volume` and must contain all of the above arguments.

## Example

The following snippet will print a JSON string containing the
configuration for two different volumes. The first volume applies
gaussian smearing to the first two spatial dimensions and uniform
smearing to the time dimension. Correspondingly, the parameter list
for this volume has five entries: two for the simple Gaussians and
three for the uniform smearer. After this, another volume gets only
one spatial dimension smeared by a simple gaussian.

```bash
Examples/Algorithms/Digitization/scripts/smearing-config.py \
   --digi-smear-volume=8 \
   --digi-smear-indices=0:1:5 \
   --digi-smear-types=0:0:3 \
   --digi-smear-parameters=10:20:2.5:-25:25 \
   --digi-smear-volume=11 \
   --digi-smear-indices=1 \
   --digi-smear-types=0 \
   --digi-smear-parameters=12.5
```
