// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTS/Digitization/CartesianSegmentation.hpp"
#include "ACTS/Digitization/DigitizationModule.hpp"
#include "ACTS/Tools/ISpacePointBuilder.hpp"
#include "ACTS/Tools/SingleHitSpacePointBuilder.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

/// @brief Structure for easier bookkeeping of hits.
struct DoubleHitSpacePoint
{
  /// Storage of the hit cluster on a surface
  std::pair<PlanarModuleCluster const*, PlanarModuleCluster const*> hitModuleFront;
  /// Storage of the hit cluster on another surface
  std::pair<PlanarModuleCluster const*, PlanarModuleCluster const*>  hitModuleBack;
  /// Combined coordinate of a cluster of hits on the front side
  Vector3D clusterPointFront = {0., 0., 0.};
  /// Combined coordinate of a cluster of hits on the back side
  Vector3D clusterPointBack = {0., 0., 0.};
  /// Storage of a space point. Zero vector indicates unset point
  Vector3D spacePoint = {0., 0., 0.};
};

/// @brief Configuration of the class to steer its behaviour
struct DoubleHitSpacePointConfig
{
  /// Accepted difference in eta for two hits
  double diffTheta2 = 1.;
  /// Accepted difference in phi for two hits
  double diffPhi2 = 1.;
  /// Accepted distance between two hits
  double diffDist = 100. * units::_mm;
  /// Allowed increase of strip length
  double stripLengthTolerance = 0.01;
  /// Allowed increase of strip length wrt gaps between strips
  double stripLengthGapTolerance = 0.01;
  /// Assumed position of the vertex
  Vector3D vertex = {0., 0., 0.};
  /// Perform the perpendicular projection for space point finding
  bool usePerpProj = false;
  /// Cluster the hits on the front side strip module
  bool clusterFrontHits = true;
  /// Cluster the hits on the back side strip module
  bool clusterBackHits = true;
};

/// @class TwoHitsSpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the hits strip detectors need further treatment. This class takes
/// the digitized hits and combines them on two different detector elements to a
/// result of the combined detector element. The class is intended to handle
/// strip detector elements in particular.
///
/// @note Used abbreviation: "Strip Detector Element" -> SDE
///
template <>
class SpacePointBuilder<DoubleHitSpacePoint, DoubleHitSpacePointConfig>
    : public SpacePointBuilder<SingleHitSpacePoint, void>
{
public:
  /// Default constructor
  SpacePointBuilder<DoubleHitSpacePoint, DoubleHitSpacePointConfig>() = delete;

  /// @brief This function is intended to use a single hit for the formation of
  /// a space point. Since this is not needed for this class this function is
  /// deleted.
  static void
  addHits(std::vector<DoubleHitSpacePoint>&              spacePointStorage,
          const std::vector<PlanarModuleCluster const*>& hits)
      = delete;

  /// @brief Searches possible combinations of two hits on different surfaces
  /// that may come from the same particles
  /// @param spacePointStorage storage of the space points
  /// @param hits1 vector of hits on a surface
  /// @param hits2 vector of hits on another surface
  /// @param cfg optional configuration to steer the combination process of @p
  /// hits1 and @p hits2
  /// @note The structure of @p hits is meant to be hits[Surfaces][Hits on a
  /// surface]
  static void
  addHits(std::vector<DoubleHitSpacePoint>&                spacePointStorage,
          const std::vector<PlanarModuleCluster const*>&   hits1,
          const std::vector<PlanarModuleCluster const*>&   hits2,
          const std::shared_ptr<DoubleHitSpacePointConfig> cfg = nullptr);

  /// @brief Calculates the space points out of a given collection of hits
  /// on several strip detectors and stores the data
  /// @param spacePointStorage storage of the data
  /// @param cfg optional configuration to steer the calculation of the space
  /// points
  /// @note If no configuration is set, the default values will be used
  static void
  calculateSpacePoints(std::vector<DoubleHitSpacePoint>& spacePoints,
                       const std::shared_ptr<DoubleHitSpacePointConfig> cfg
                       = nullptr);

private:
  /// @brief Storage container for variables related to the calculation of space
  /// points
  struct SpacePointParameters
  {
    /// Vector pointing from bottom to top end of first SDE
    Vector3D q;
    /// Vector pointing from bottom to top end of second SDE
    Vector3D r;
    /// Twice the vector pointing from vertex to to midpoint of first SDE
    Vector3D s;
    /// Twice the vector pointing from vertex to to midpoint of second SDE
    Vector3D t;
    /// Cross product between SpacePointParameters::q and
    /// SpacePointParameters::s
    Vector3D qs;
    /// Cross product between SpacePointParameters::r and
    /// SpacePointParameters::t
    Vector3D rt;
    /// Magnitude of SpacePointParameters::q
    double qmag = 0.;
    /// Parameter that determines the hit position on the first SDE
    double m = 0.;
    /// Parameter that determines the hit position on the second SDE
    double n = 0.;
    /// Regular limit of the absolut values of SpacePointParameters::m and
    /// SpacePointParameters::n
    double limit = 1.;
    /// Limit of SpacePointParameters::m and SpacePointParameters::n in case of
    /// variable vertex
    double limitExtended = 1.;

    /// @brief reset structure and allows to reuse the same structure
    void
    reset()
    {
      // Set every vector to nullvector. This allows checks, if a variable was
      // already set.
      q  = {0., 0., 0.};
      r  = {0., 0., 0.};
      s  = {0., 0., 0.};
      t  = {0., 0., 0.};
      qs = {0., 0., 0.};
      rt = {0., 0., 0.};
      // Set every double to default values
      qmag          = 0;
      m             = 0;
      n             = 0;
      limit         = 1.;
      limitExtended = 1.;
    }
  };

  /// @brief Calculates (Delta theta)^2 + (Delta phi)^2 between two hits
  /// @param pos1 the first hit
  /// @param pos2 the second hit
  /// @param cfg optional configuration to steer the combination process of @p
  /// hit1 and @p hit2
  /// @return the squared sum in case of success, otherwise -1
  static double
  differenceOfHits(const Vector3D&                       pos1,
                   const Vector3D&                       pos2,
                   const std::shared_ptr<DoubleHitSpacePointConfig> cfg);

  static std::vector<BinningData>
  binningData(PlanarModuleCluster const* hit);

  /// @brief Calculates the bin of a hit
  /// @param hit recorded hit
  /// @return channel 0 and 1 of the hit
  static std::pair<size_t, size_t>
  binOfHit(PlanarModuleCluster const* hit);
  
  /// @brief Create a sparse matrix of hits based on the bin numbers as look up
  /// @param hits list of hits
  /// @return matrix with hits
  /// @note If the hits are given from more than one surface the matrix will be empty. This prevents false listing since it is based on bin indices.
  static std::vector<std::vector<PlanarModuleCluster const*>>
  sortHits(const std::vector<PlanarModuleCluster const*>& hits);
                        

  /// @brief Build pair of hits on neighboring bins
  /// @param hits collection of hits on a single surface
  /// @param peformClustering configuration that steers if a clustering should be performed
  /// @return collection of found pairs
  static const std::vector<std::pair<PlanarModuleCluster const*,
                                     PlanarModuleCluster const*>>
  clusterSpacePoints(
      const std::vector<PlanarModuleCluster const*>& hits, bool performClustering);

	/// @brief Calculate the mean of the coordinates of two hits
	/// @param cluster pair of neighbouring hits
	/// @return vector with the combined coordinates
  static const Vector3D
  clusterPoint(std::pair<PlanarModuleCluster const*,
                                     PlanarModuleCluster const*> cluster);

  /// @brief Calculates the top and bottom ends of a SDE
  /// that corresponds to a given hit
  /// @param hit object that stores the information about the hit
  /// @return vectors to the top and bottom end of the SDE
  static std::pair<Vector3D, Vector3D>
  endsOfStrip(const PlanarModuleCluster& hit);

  /// @brief Calculates a space point whithout using the vertex
  /// @note This is mostly to resolve space points from cosmic data
  /// @param a vector to the top end of the first SDE
  /// @param c vector to the top end of the second SDE
  /// @param q vector from the bottom to the top end of the first SDE
  /// @param r vector from the bottom to the top end of the second SDE
  /// @return parameter that indicates the location of the space point; returns
  /// 1. if it failed
  /// @note The meaning of the parameter is explained in more detail in the
  /// function body
  static double
  calcPerpProj(const Vector3D& a,
               const Vector3D& c,
               const Vector3D& q,
               const Vector3D& r);

  /// @brief This function tests if a space point can be estimated by a more
  /// tolerant treatment of construction. In fact, this function indirectly
  /// allows shifts of the vertex.
  /// @param spaPoPa container that stores geometric parameters and rules of the
  /// space point formation
  /// @param cfg optional configuration to steer the recovery of space points
  /// @return indicator if the test was successful
  static bool
  recoverSpacePoint(SpacePointParameters&                            spaPoPa,
                    const std::shared_ptr<DoubleHitSpacePointConfig> cfg);
};

}  // namespace Acts
