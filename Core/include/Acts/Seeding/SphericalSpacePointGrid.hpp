// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <numbers>
#include <vector>

namespace Acts::Experimental {

/// A spherical space point grid for seeding, binned in (phi, eta, r).
///
/// A space point's eta bin is found from cot(theta) = z / r (a single division),
/// matched against bin edges that are stored as sinh(eta). Since
/// cot(theta) = sinh(eta), comparing z / r against those edges
/// places the point into exactly its eta bin, with no atan2 / tan / log / asinh.
class SphericalSpacePointGrid {
 public:
  /// Space point index type used in the grid.
  using SpacePointIndex = std::uint32_t;
  /// Type alias for bin container holding space point indices
  using BinType = std::vector<SpacePointIndex>;
  /// Type alias for phi axis with equidistant binning and closed boundaries
  using PhiAxisType = Axis<AxisType::Equidistant, AxisBoundaryType::Closed>;
  /// Type alias for the middle (cot(theta)) axis with variable binning and open
  /// boundaries. Edges are stored as sinh(eta) = cot(theta), so the binning is
  /// configured in eta while space points are placed by z / r (see class doc).
  using CotThetaAxisType = Axis<AxisType::Variable, AxisBoundaryType::Open>;
  /// Type alias for r axis with variable binning and open boundaries
  using RAxisType = Axis<AxisType::Variable, AxisBoundaryType::Open>;
  /// Spherical space point grid type (phi, eta, r)
  using GridType = Grid<BinType, PhiAxisType, CotThetaAxisType, RAxisType>;
  /// Type alias for binned group over the spherical grid
  using BinnedGroupType = BinnedGroup<GridType>;

  /// Configuration parameters for the spherical space point grid.
  struct Config {
    /// minimum pT
    float minPt = 0 * UnitConstants::MeV;
    /// minimum extension of sensitive detector layer relevant for seeding as
    /// distance from x=y=0 (i.e. in r)
    /// WARNING: if rMin is smaller than impactMax, the bin size will be 2*pi,
    /// which will make seeding very slow!
    float rMin = 0 * UnitConstants::mm;
    /// maximum extension of sensitive detector layer relevant for seeding as
    /// distance from x=y=0 (i.e. in r)
    float rMax = 600 * UnitConstants::mm;
    /// minimum pseudorapidity relevant for seeding
    float etaMin = -4;
    /// maximum pseudorapidity relevant for seeding
    float etaMax = 4;
    /// maximum distance in r from middle space point to bottom or top
    /// space point
    float deltaRMax = 270 * UnitConstants::mm;
    /// maximum impact parameter in mm
    float impactMax = 0 * UnitConstants::mm;
    /// minimum phi value for phiAxis construction
    float phiMin = -std::numbers::pi_v<float>;
    /// maximum phi value for phiAxis construction
    float phiMax = std::numbers::pi_v<float>;
    /// Multiplicator for the number of phi-bins. The minimum number of phi-bins
    /// depends on min_pt, magnetic field: 2*pi/(minPT particle phi-deflection).
    /// phiBinDeflectionCoverage is a multiplier for this number. If
    /// numPhiNeighbors (in the configuration of the BinFinders) is configured
    /// to return 1 neighbor on either side of the current phi-bin (and you want
    /// to cover the full phi-range of minPT), leave this at 1.
    int phiBinDeflectionCoverage = 1;
    /// maximum number of phi bins
    int maxPhiBins = 10000;
    /// width of the equidistant eta bins used when `etaBinEdges` is empty
    float deltaEtaMax = 0.8;
    /// enable non equidistant binning in eta (edges given in pseudorapidity;
    /// mapped to cot(theta) = sinh(eta) internally)
    std::vector<float> etaBinEdges{};
    /// enable non equidistant binning in r
    std::vector<float> rBinEdges{};

    /// magnetic field
    float bFieldInZ = 0 * UnitConstants::T;

    /// bin finder for bottom space points
    std::optional<GridBinFinder<3ul>> bottomBinFinder;
    /// bin finder for top space points
    std::optional<GridBinFinder<3ul>> topBinFinder;
    /// navigation structure for the grid
    std::array<std::vector<std::size_t>, 3ul> navigation;
  };

  /// Construct a spherical space point grid with the given configuration and
  /// an optional logger.
  /// @param config Configuration for the spherical grid
  /// @param logger Optional logger instance for debugging output
  explicit SphericalSpacePointGrid(
      const Config& config,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("SphericalSpacePointGrid", Logging::Level::INFO));

  /// Clear the grid and drop all state. The object will behave like a newly
  /// constructed one.
  void clear();

  /// Get the bin index for a space point given its azimuthal angle,
  /// cot(theta) = z / r, and radial distance.
  /// @param position The position of the space point in (phi, cotTheta, r)
  /// coordinates
  /// @return The index of the bin in which the space point is located, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> binIndex(const Vector3& position) const {
    if (!grid().isInside(position)) {
      return std::nullopt;
    }
    return grid().globalBinFromPosition(position);
  }
  /// Get the bin index for a space point.
  /// @param phi The azimuthal angle of the space point in radians
  /// @param cotTheta cot(theta) = z / r, used to place the point in its eta bin
  /// via cot(theta) = z / r = sinh(eta)
  /// @param r The radial distance of the space point from the origin
  /// @return The index of the bin in which the space point is located, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> binIndex(float phi, float cotTheta, float r) const {
    return binIndex(Vector3(phi, cotTheta, r));
  }

  /// Insert a space point into the grid.
  /// @param index The index of the space point to insert
  /// @param phi The azimuthal angle of the space point in radians
  /// @param cotTheta cot(theta) = z / r, used to place the point in its eta bin
  /// via cot(theta) = z / r = sinh(eta)
  /// @param r The radial distance of the space point from the origin
  /// @return The index of the bin in which the space point was inserted, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> insert(SpacePointIndex index, float phi,
                                    float cotTheta, float r);
  /// Insert a space point into the grid, binning it by cot(theta) = z / r
  /// (a single division; assumes r > 0); see the class documentation.
  /// @param sp The space point to insert
  /// @return The index of the bin in which the space point was inserted, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> insert(const ConstSpacePointProxy& sp) {
    return insert(sp.index(), sp.phi(), sp.z() / sp.r(), sp.r());
  }

  /// Fill the grid with space points from the container.
  /// @param spacePoints The space point container to fill the grid with
  void extend(const SpacePointContainer::ConstRange& spacePoints);

  /// Sort the bins in the grid by the space point radius, which is required by
  /// some algorithms that operate on the grid.
  /// @param spacePoints The space point container to sort the bins by radius
  void sortBinsByR(const SpacePointContainer& spacePoints);

  /// Compute the range of radii in the grid. This requires the grid to be
  /// filled with space points and sorted by radius. The sorting can be done
  /// with the `sortBinsByR` method.
  /// @param spacePoints The space point container to compute the radius range
  /// @return The range of radii in the grid
  Range1D<float> computeRadiusRange(
      const SpacePointContainer& spacePoints) const;

  /// Mutable bin access by index.
  /// @param index The index of the bin to access
  /// @return Mutable reference to the bin at the specified index
  BinType& at(std::size_t index) { return grid().at(index); }
  /// Const bin access by index.
  /// @param index The index of the bin to access
  /// @return Const reference to the bin at the specified index
  const BinType& at(std::size_t index) const { return grid().at(index); }

  /// Mutable grid access.
  /// @return Mutable reference to the grid
  GridType& grid() { return *m_grid; }
  /// Const grid access.
  /// @return Const reference to the grid
  const GridType& grid() const { return *m_grid; }

  /// Access to the binned group.
  /// @return Reference to the binned group
  const BinnedGroupType& binnedGroup() const { return *m_binnedGroup; }

  /// Get the number of space points in the grid.
  /// @return The number of space points in the grid
  std::size_t numberOfSpacePoints() const { return m_counter; }

  /// Get the number of bins in the grid.
  /// @return The number of bins in the grid
  std::size_t numberOfBins() const { return grid().size(); }

 private:
  Config m_cfg;

  std::unique_ptr<const Logger> m_logger;

  GridType* m_grid{};
  std::optional<BinnedGroupType> m_binnedGroup;

  std::size_t m_counter{};

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts::Experimental
