// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/CosmicSeedCuts.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace Acts {

/// @brief Seed finder for cosmic-ray tracks.
///
/// Builds top-middle-bottom triplets ordered along the global y axis and keeps
/// those that pass the configured CosmicSeedCuts. It can also build two space
/// point seeds (bottom, top) for tracks with too few hits for a triplet.
class CosmicTripletSeedFinder {
 public:
  /// Check the cuts and precompute the quantities derived from them.
  /// @param cuts Cuts applied to the seeds
  /// @throws std::invalid_argument if the cuts are inconsistent
  explicit CosmicTripletSeedFinder(const CosmicSeedCuts& cuts = {});

  /// Create triplet seeds from the space points of one event.
  /// @param spacePoints Input space points
  /// @param outputSeeds Container the seeds are appended to. If it has no
  ///        space point container yet, @p spacePoints is assigned to it.
  /// @param spClasses One class per space point, same indices as
  ///        @p spacePoints. Required with useCovarianceScale,
  ///        requireSameSubdetector or useSubdetectorAngularCuts.
  /// @param calibratedSpacePoints Copy of @p spacePoints with the same
  ///        indices. With useStripCalibration, each strip space point used in a
  ///        seed is moved to the position computed with the direction of its
  ///        best seed. May be null.
  /// @param spSurfaces Sensitive surface of each space point, same indices as
  ///        @p spacePoints. Required with rejectSameSurface.
  void createSeeds(
      const SpacePointContainer& spacePoints, SeedContainer& outputSeeds,
      const std::vector<CosmicSpacePointClass>* spClasses = nullptr,
      SpacePointContainer* calibratedSpacePoints = nullptr,
      const std::vector<GeometryIdentifier>* spSurfaces = nullptr) const;

  /// Create two space point seeds (bottom, top), like the find2Sp strategy of
  /// Athena. The kink cut is not applied, since there is no third point.
  /// @param spacePoints Input space points
  /// @param outputSeeds Container the seeds are appended to. If it has no
  ///        space point container yet, @p spacePoints is assigned to it.
  /// @param spClasses One class per space point, same indices as
  ///        @p spacePoints. Required with requireSameSubdetector.
  /// @param skipSpacePoints One flag per space point; space points with a
  ///        non-zero flag are not used (for example those already in a
  ///        triplet seed). May be null.
  /// @param spSurfaces Sensitive surface of each space point, same indices as
  ///        @p spacePoints. Required with rejectSameSurface.
  void createTwoSpacePointSeeds(
      const SpacePointContainer& spacePoints, SeedContainer& outputSeeds,
      const std::vector<CosmicSpacePointClass>* spClasses = nullptr,
      const std::vector<std::uint8_t>* skipSpacePoints = nullptr,
      const std::vector<GeometryIdentifier>* spSurfaces = nullptr) const;

 private:
  // position of a SP in a triplet: the one from the container,
  // or the recalculated one with the direction of the triplet if it is a strip
  // SP
  struct Point {
    double x = 0., y = 0., z = 0.;
    std::size_t index = 0;
  };

  // information about the event, derived from the space points and their
  // classes
  struct EventData {
    const SpacePointContainer* spacePoints = nullptr;
    const double* covR = nullptr;
    const double* covZ = nullptr;
    const std::uint8_t* isStrip = nullptr;
    const std::uint8_t* isCalibratable = nullptr;
    const OuterStripSpacePointCalibrationDetailsDerived* stripDetails = nullptr;
    const GeometryIdentifier* surfaces = nullptr;
  };

  struct EventBuffers {
    std::vector<double> covR;
    std::vector<double> covZ;
    std::vector<std::uint8_t> isStrip;
    std::vector<std::uint8_t> isCalibratable;
    std::vector<OuterStripSpacePointCalibrationDetailsDerived> stripDetails;
    // the extent to which the y-coordinate of a strip SP changes upon
    // recalculation [mm]
    double stripShiftY = 0.0;
  };

  // fill buffers() with the necessary information for each space point,
  // depending on the cuts that are enabled
  EventData prepareEvent(const SpacePointContainer& spacePoints,
                         const std::vector<CosmicSpacePointClass>* spClasses,
                         const std::vector<GeometryIdentifier>* spSurfaces,
                         bool needCovariance, bool needType,
                         EventBuffers& buffers) const;

  bool calibrateStripSpacePoints(Point& top, Point* middle, Point& bottom,
                                 const EventData& ev) const;

  bool passDoubletCuts(const Point& top, const Point& bottom) const;

  bool tripletDirections(const Point& top, const Point& middle,
                         const Point& bottom,
                         std::array<std::array<float, 3>, 3>& dirs) const;

  bool passCosmicCuts(const Point& top, const Point& middle,
                      const Point& bottom, float& quality,
                      const EventData& ev) const;

  bool passCurvatureAndQuality(const Point& top, const Point& middle,
                               const Point& bottom, float& quality,
                               const EventData& ev) const;

  struct PtDependentTerms {
    double polarKink = 0.;         // |theta1 - theta2| [rad]
    double kinkSlack = 0.;         // multiple scattering in angle [rad]
    double parallelResidual = 0.;  // residual out of the transverse plane [mm]
    double residualSlack = 0.;     // multiple scattering in residual [mm]
  };
  bool ptDependentTerms(const Point& top, const Point& middle,
                        const Point& bottom, PtDependentTerms& out) const;

  static double residualToTopBottomLine(const Point& top, const Point& middle,
                                        const Point& bottom);

  static float estimateVertexZ(const Point& top, const Point& middle,
                               const Point& bottom);

  CosmicSeedCuts m_cuts;
  double m_dzdrMax = 0.0;

  bool m_ptDependent = false;
  bool m_ptPrefilter = false;
  // minimal radius of curvature for the prefilter: 0.9 minPt / (0.3 |Bz|) [mm]
  double m_rhoMin = 0.0;
  // cosine of the maximum kink angle for the prefilter: 1 - k^2/2 <= cos(k)
  // for every k in [0, maxKinkAngle]
  double m_kinkCosThreshold = -1.0;
  // lower limit of the cosine of the angle between each segment and the y-axis
  // for the prefilter
  double m_segCosThetaYMin = 0.0;
};

}  // namespace Acts
