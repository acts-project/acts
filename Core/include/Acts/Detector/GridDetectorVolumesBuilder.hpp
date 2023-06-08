// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Detector/detail/CylindricalGridVolumesHelper.hpp"
#include "Acts/Detector/detail/IndexedGridFiller.hpp"
#include "Acts/Detector/interface/IDetectorComponentBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

#include <array>
#include <vector>

namespace Acts {

class Surface;

namespace Experimental {

class Detector;

struct NoSurfaces {
  std::vector<std::shared_ptr<Surface>> operator()() const { return {}; }
};

/// @brief A dedicated container builder for detectors based on a grid structure
/// It supports the constructruction of:
///
/// - a cylindrical grid in z/r/phi (and subsets) - with trapezoidal
/// approximation
/// - a cartesian grid in x/y/z (and subsets) (not yet implemented)
///
class GridDetectorVolumesBuilder : public IDetectorComponentBuilder {
 public:
  /// @brief Nested configuration struct
  struct Config {
    // Name of the detector builder
    std::string name = "";
    /// Connection point for a function to provide surfaces
    std::function<std::vector<std::shared_ptr<Surface>>()> surfaces =
        NoSurfaces{};
    /// Binning description
    std::vector<ProtoBinning> binning = {};
    /// The segments for the polyhedron reference generator
    unsigned int polyhedronSegements = 4;
    /// Use a polygon approximation for the cylinders
    bool approximateCylinders = false;
    /// The innermost is a full bin
    bool innerMostSingleBin = false;
    /// A gobal transform
    Transform3 transform = Transform3::Identity();
    /// Auxiliary information
    std::string auxillary = "";
  };

  /// Constructor
  ///
  /// @param cfg is the configuration struct
  /// @param logger logging instance for screen output
  GridDetectorVolumesBuilder(const Config& cfg,
                             std::unique_ptr<const Logger> logger =
                                 getDefaultLogger("GridDetectorVolumesBuilder",
                                                  Logging::INFO));

  /// Final implementation of a volume builder that is purely defined
  /// by an internal and external structure builder
  ///
  /// @param gctx The geometry context for this call
  ///
  /// @return an outgoing detector component
  DetectorComponent construct(const GeometryContext& gctx) const final;

 private:
  /// @brief templated construction method depending on the axis generator
  ///
  /// @tparam axis_generator_t the axis generator type
  //
  /// @param gctx The geometry context for this call
  /// @param ag the axis generator instance
  ///
  /// @return the detector compontns object
  template <typename axis_generator_t>
  DetectorComponent constructT(const GeometryContext& gctx,
                               const axis_generator_t& ag) const {
    // Define a matching grid type dependent on the axis generator
    using GridType =
        typename axis_generator_t::template grid_type<std::vector<std::size_t>>;
    GridType grid(ag());

    auto surfaces = m_cfg.surfaces();
    ACTS_DEBUG("Retreived " << surfaces.size()
                            << " surfaces to assign to the grid.");

    IndexedSurfacesImpl<GridType> isGrid(std::move(grid), m_binningValues,
                                         m_cfg.transform);
    detail::IndexedGridFiller isGridFiller;
    detail::PolyhedronReferenceGenerator polyGen{true,
                                                 m_cfg.polyhedronSegements};
    isGridFiller.fill(gctx, isGrid, surfaces, polyGen);

    detail::CylindricalGridVolumesHelper::Options cgOptions;
    cgOptions.polygonApproximation = m_cfg.approximateCylinders;

    // Construct the grid volumes, also fills the root volumes grid
    auto [volumes, portalContainer, rootVolumesGrid] =
        detail::CylindricalGridVolumesHelper::buildVolumes<GridType,
                                                           axis_generator_t>(
            gctx, isGrid.grid, ag, surfaces, cgOptions);

    ACTS_DEBUG("Retreived " << volumes.size() << " grid volumes.");

    // The detector component is ready and can be returned
    return DetectorComponent{
        volumes, portalContainer,
        RootDetectorVolumes{
            volumes, rootVolumesFromGrid(std::move(rootVolumesGrid),
                                         m_binningValues, m_cfg.transform)}};
  }

  /// Filter-check 3 x equidistant
  static constexpr std::array<Acts::detail::AxisType, 3u> eeeType = {
      Acts::detail::AxisType::Equidistant, Acts::detail::AxisType::Equidistant,
      Acts::detail::AxisType::Equidistant};

  /// Filter-check 2 x arbitrary, 1 x equidistant
  static constexpr std::array<Acts::detail::AxisType, 3u> vveType = {
      Acts::detail::AxisType::Variable, Acts::detail::AxisType::Variable,
      Acts::detail::AxisType::Equidistant};

  /// Filter-check 3 x variable
  static constexpr std::array<Acts::detail::AxisType, 3u> vvvType = {
      Acts::detail::AxisType::Variable, Acts::detail::AxisType::Variable,
      Acts::detail::AxisType::Variable};

  /// Filter-check cylindrical binning
  static constexpr std::array<BinningValue, 3u> cylindricalBinning = {
      BinningValue::binZ, BinningValue::binR, BinningValue::binPhi};

  /// Filter-check cartesian binning
  static constexpr std::array<BinningValue, 3u> cartesianBinning = {
      BinningValue::binX, BinningValue::binY, BinningValue::binZ};

  /// Configuration object
  Config m_cfg;

  /// Binning values from parsing
  std::array<BinningValue, 3u> m_binningValues = {};

  /// Private acces method to the logger
  const Logger& logger() const { return *m_logger; }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Experimental
}  // namespace Acts
