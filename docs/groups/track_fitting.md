@defgroup track_fitting Track Fitting
@brief Algorithms and data structures for track fitting

The track fitters estimate the track parameters from a fixed set of
measurements. Material effects (multiple scattering and energy loss) are applied
on each surface using a shared pointwise formalism; see
@ref kf_material_effects for the details of that formalism and its
implementation.

For electrons, whose Bremsstrahlung energy loss is strongly non-Gaussian, the
@ref Acts::GaussianSumFitter "Gaussian Sum Filter" models the track state as a
mixture of Gaussians instead; see @ref gsf.

The @ref Acts::TrackletJoiner joins a second tracklet, filtered from a state of
a first one, and combines their smoothed estimates without a new propagation;
this is the two-way track finding.
