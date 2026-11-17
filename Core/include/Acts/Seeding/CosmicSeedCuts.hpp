// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <cstdint>

namespace Acts {

/// Detector region of a space point, as seen by the cosmic seed finder.
///
/// It selects the covariance scale factor and tells pixel space points from
/// strip space points for the subdetector cuts.
enum class CosmicSpacePointClass : std::uint8_t {
  BarrelPixel = 0,  ///< Pixel space point in the barrel
  EndcapPixel = 1,  ///< Pixel space point in an endcap
  BarrelStrip = 2,  ///< Strip space point in the barrel
  EndcapStrip = 3,  ///< Strip space point in an endcap
};

/// @brief Cuts applied by CosmicTripletSeedFinder when it forms
/// top-middle-bottom triplets from cosmic-ray space points.
///
/// Cosmic muons cross the detector from top to bottom, so the seed finder
/// orders and separates space points along the global y axis instead of the
/// radius. Lengths are in mm, momenta in MeV and the magnetic field in T.
struct CosmicSeedCuts {
  /// Require the vertical order top.y > middle.y > bottom.y.
  bool requireTopMiddleBottom = true;
  /// Maximum kink angle between the two segments of a triplet [rad]. With
  /// useSubdetectorAngularCuts it applies only to triplets with a strip middle
  /// space point; pixel ones are cut by pixelMinCosDeltaPhi instead.
  double maxKinkAngle = 0.05;
  /// Minimum vertical separation of each segment [mm]: |y_top - y_middle| and
  /// |y_middle - y_bottom| must lie in [minDeltaY, maxDeltaY]. Together with
  /// maxDeltaY it plays the role of mindRadius/maxdRadius in Athena's
  /// SiSpacePointsSeedMaker_Cosmic, where dR is dy.
  double minDeltaY = 5.0;
  /// Maximum vertical separation of each segment [mm], see minDeltaY.
  double maxDeltaY = 1000.0;
  /// Minimum alignment of each segment with the y axis, |dy| / |dr| with
  /// |dr| = sqrt(dx^2 + dy^2 + dz^2).
  double minCosThetaY = 0.3;
  /// Maximum distance of the middle space point to the top-bottom line [mm].
  /// With usePtDependentKinkResidual it is the lower bound of the tolerance on
  /// the component out of the transverse plane. A negative value disables the
  /// cut.
  double maxResidual = 5.0;
  /// Minimum z of the space points used for seeding [mm].
  double zMin = -3000.0;
  /// Maximum z of the space points used for seeding [mm].
  double zMax = 3000.0;
  /// Acceptance |eta| < maxAbsEta, applied as a window in dz/dr.
  double maxAbsEta = 2.3;
  /// Minimum transverse momentum [MeV], used only with usePtCut and a
  /// non-zero magneticFieldBz.
  double minPt = 300.0;
  /// Apply the minPt cut, with pT estimated from the triplet curvature.
  bool usePtCut = false;
  /// Magnetic field along z used to estimate pT [T].
  double magneticFieldBz = 2.0;
  /// Multiple scattering coefficient used in the longitudinal compatibility
  /// and in the pT-dependent tolerances.
  double multipleScatteringCof = 60.3;
  /// Reject triplets whose longitudinal discriminant dT is positive, i.e.
  /// incompatible with a straight (or curved) trajectory within the errors.
  bool useLongitudinalCompatibility = false;
  /// Cut on the polar kink and on the residual out of the transverse plane,
  /// with tolerances widened by multiple scattering at the estimated pT.
  /// Requires usePtCut, minPt > 0 and magneticFieldBz != 0.
  bool usePtDependentKinkResidual = false;
  /// Largest pT used in the scattering term [MeV]: a higher estimate is
  /// replaced by this value, as in Athena.
  double maxPtForScattering = 5000.0;

  /// Scale the space point variances with covarianceScale. On by default, as
  /// in Athena; it can be disabled for tests.
  bool useCovarianceScale = true;
  /// Variance scale factor per CosmicSpacePointClass, indexed by its value.
  /// Only the dominant coordinate is scaled: z in the barrel, r in the
  /// endcaps.
  std::array<double, 4> covarianceScale = {16., 16., 100., 16.};

  /// Require all space points of a seed to be of the same type, all pixel or
  /// all strip. Needs the space point classes.
  bool requireSameSubdetector = false;

  /// Reject seeds in which top and middle, or middle and bottom, lie on the
  /// same surface (module). Needs the surface of each space point in
  /// CosmicTripletSeedFinder::createSeeds().
  bool rejectSameSurface = true;

  /// Choose the angular cut between the two segments by the type of the
  /// middle space point: maxKinkAngle for strips, pixelMinCosDeltaPhi for
  /// pixels. Needs the space point classes.
  bool useSubdetectorAngularCuts = false;
  /// Minimum cosine of the transverse angle between the two segments of a
  /// triplet with a pixel middle space point, see useSubdetectorAngularCuts.
  double pixelMinCosDeltaPhi = 0.96;

  /// Recompute the position of the strip space points in the seeds with the
  /// direction of their best seed. Needs a copy of the space points with the
  /// same indices in CosmicTripletSeedFinder::createSeeds().
  bool useStripCalibration = false;
  /// Tolerance on the position along the strips when a strip space point is
  /// recomputed, relative to the strip half length (1 means inside the strip).
  double stripCalibrationTolerance = 1.1;

  /// Compute the seed quality from dT and the segment lengths. Otherwise every
  /// seed gets seedQuality. The quality is used to sort and filter seeds in
  /// the CKF.
  bool useSeedQuality = true;
  /// Quality of the seeds when useSeedQuality is false, and of all two space
  /// point seeds.
  float seedQuality = 0.0f;
};

}  // namespace Acts
