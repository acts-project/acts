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
#include <utility>
#include <vector>

namespace Acts {

// class of a space point for the purpose of covariance scaling
enum class CosmicSpacePointClass : std::uint8_t {
  BarrelPixel = 0,
  EndcapPixel = 1,
  BarrelStrip = 2,
  EndcapStrip = 3,
};

// cuts applied by CosmicTripletSeedFinder when forming top-middle-bottom
// triplets from cosmic ray space points
struct CosmicSeedCuts {
  // order: top.y() > middle.y() > bottom.y()
  bool requireTopMiddleBottom = true;
  // with useSubdetectorAngularCuts, it is the maximum kink for strip SPs
  // and pixel SPs are cut by pixelMinCosDeltaPhi instead
  double maxKinkAngle = 0.05;
  // vertical separation for each segment [mm], for any middle SP (pixel or
  // strip): |y_top - y_middle| and |y_middle - y_bottom| in [minDeltaY,
  // maxDeltaY] these are mindRadius/maxdRadius of SiSpacePointsSeedMaker_Cosmic
  // (m_drmin = 5, m_drmax = 1000), where dR is dy, and the deltaRMin/deltaRMax
  // of the Acts::DoubletSeedFinder that makes the doublets (its r is y for
  // cosmics, see CosmicTripletSeedFinder)
  double minDeltaY = 5.0;
  double maxDeltaY = 1000.0;
  // minimal alignment with the y direction: |dy| / |dr|, where dr = sqrt(dx^2 +
  // dy^2 + dz^2)
  double minCosThetaY = 0.3;
  // maximal distance of the middle SP to the line top-bottom [mm]
  // with usePtDependentKinkResidual it is the floor of the tolerance of the
  // component out of the transverse plane (see passGeometricCuts)
  // negative = no cut
  double maxResidual = 5.0;
  // zMin/zMax [mm] of the space points to be considered for seeding
  double zMin = -3000.0;
  double zMax = 3000.0;
  // Athena's binning, used only with useCartesianGridBinning = false:
  // bins in z [mm]: the middle SP is taken from one bin and its top and bottom
  // candidates from the neighbouring bins (buildFrameWork() / production3Sp of
  // SiSpacePointsSeedMaker_Cosmic). zBinEdges are the boundaries between bins,
  // strictly increasing: N edges give N + 1 bins, the first and last open and
  // bounded only by zMin/zMax. A z exactly on an edge belongs to the bin
  // farther from z = 0, as in Athena. The default is Athena's binning.
  std::vector<double> zBinEdges = {-2500., -1400., -925., -450., -250.,
                                   250.,   450.,   925.,  1400., 2500.};
  // neighbouring bins of each bin
  std::vector<std::pair<int, int>> zBinNeighbors = {
      {0, 1},  {-1, 1}, {-1, 1}, {-1, 2}, {-1, 1}, {-1, 1},
      {-1, 1}, {-2, 1}, {-1, 1}, {-1, 1}, {-1, 0}};
  // binning of the space points (zBinEdges and zBinNeighbors are then not
  // used): equidistant bins in x, y and z [mm] of an
  // Acts::CartesianSpacePointGrid, gridBins = {nx, ny, nz} bins over
  // [gridXMin, gridXMax], [gridYMin, gridYMax] and [zMin, zMax]. The middle SP
  // is taken from one bin and its top and bottom candidates from the bins up
  // to gridNeighbors = {x, y, z} bins away on each side, per axis. A space
  // point outside the x or y range takes no part. Neighbours that do not reach
  // maxDeltaY in y, or the x and z separation of the doublets of an inclined
  // muon, lose seeds. The default x/y range covers the ODD (r < 1080 mm).
  // The default bins only z, 12 bins with 3 neighbours on each side: it finds
  // every muon that Athena's z bins find, and more of the muons inclined in z
  // (Athena's central bin reaches only 200 mm in z)
  bool useCartesianGridBinning = true;
  std::array<std::uint32_t, 3> gridBins = {1, 1, 12};
  double gridXMin = -1100.0;
  double gridXMax = 1100.0;
  double gridYMin = -1100.0;
  double gridYMax = 1100.0;
  std::array<int, 3> gridNeighbors = {0, 0, 3};
  // acceptance in |eta| < maxAbsEta
  double maxAbsEta = 2.3;
  // momentum cut: pT > minPtFactor * minPt [MeV] (requires
  // magneticFieldBz != 0), done by Acts::TripletSeedFinder (minHelixDiameter2)
  double minPt = 300.0;
  bool usePtCut = false;
  double magneticFieldBz = 2.0;  // T
  // the pT cut is pT > minPtFactor * minPt, and without pT cut the scattering
  // term uses pT = minPtFactor * minPt: Athena's 0.9 (m_ipt = 1 / (0.9 ptmin)
  // in SiSpacePointsSeedMaker_Cosmic)
  double minPtFactor = 0.9;
  // pT [MeV] = ptPerFieldRadius * Bz [T] * R [mm], used for the pT cut, the
  // pT-dependent kink/residual and the field given to Acts::TripletSeedFinder
  // (bFieldInZ). Athena's 0.3 (K = 2 / (300 Bz[kT])) instead of the exact
  // 0.299792458, so that the pT cut is the same as in Athena
  double ptPerFieldRadius = 0.3;  // MeV / (T mm)
  // the scattering term of dT is COF (1 + cot^2) / pT^2
  double multipleScatteringCof = 60.3;
  // if true, rejects triplets whose longitudinal discriminant dT is positive,
  // i.e. incompatible with a straight (or curved) trajectory within the errors.
  // The cut is the one of Acts::TripletSeedFinder; if false, its scattering
  // terms are opened so that it does not cut on dT
  bool useLongitudinalCompatibility = false;
  bool usePtDependentKinkResidual = false;
  // maximum pT for scattering: if the estimated pT is above this value,
  // it is replaced by this value in the scattering term (same as in Athena)
  double maxPtForScattering = 5000.0;

  // active by default, but can be disabled for testing
  bool useCovarianceScale = true;
  std::array<double, 4> covarianceScale = {16., 16., 100., 16.};

  bool requireSameSubdetector = false;

  // if true, the top-middle and middle-bottom SPs cannot be on the same surface
  // (module) this requires passing the surface of each SP to createSeeds()
  bool rejectSameSurface = true;

  // if true, the cuts on the angular separation of the segments are applied
  // depending on the subdetector type (pixel or strip)
  bool useSubdetectorAngularCuts = false;
  double pixelMinCosDeltaPhi = 0.96;

  // if true, the positions of strip space points in the seeds are recalculated
  // with the direction of their best seed this requires passing a copy of the
  // space points with the same indices to createSeeds()
  bool useStripCalibration = false;
  double stripCalibrationTolerance = 1.1;
  // the recalculation is iterative: crossing with the triplet direction, new
  // direction from the new crossings, and so on, at most
  // stripCalibrationMaxIterations times, until no space point moves more than
  // stripCalibrationConvergence [mm]; without convergence the triplet is
  // rejected. It is also rejected when |cos| between the triplet direction and
  // the normal of a strip module is below stripCalibrationMinCosToNormal (a
  // direction almost parallel to the module crosses the strips far away)
  int stripCalibrationMaxIterations = 50;
  double stripCalibrationConvergence = 1e-3;
  double stripCalibrationMinCosToNormal = 1e-6;

  // quality assigned to the seeds, used for sorting and filtering in the CKF
  bool useSeedQuality = true;
  float seedQuality = 0.0f;
};

}  // namespace Acts
