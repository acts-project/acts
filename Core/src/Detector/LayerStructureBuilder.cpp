// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/LayerStructureBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Detector/detail/IndexedSurfacesGenerator.hpp"
#include "Acts/Detector/detail/ReferenceGenerators.hpp"
#include "Acts/Detector/detail/SupportHelper.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <cmath>
#include <cstddef>
#include <ostream>
#include <set>
#include <stdexcept>
#include <utility>

namespace Acts {
namespace Experimental {
class DetectorVolume;
}  // namespace Experimental
}  // namespace Acts

namespace {

/// Helper for 1-dimensional generators
///
/// @tparam aType is the axis boundary type: closed or bound
///
/// @param gctx the geometry context
/// @param lSurfaces the surfaces of the layer
/// @param assignToAll the indices assigned to all
/// @param binning the binning struct
///
/// @return a configured surface candidate updators
template <Acts::detail::AxisBoundaryType aType>
Acts::Experimental::SurfaceCandidatesUpdator createUpdator(
    const Acts::GeometryContext& gctx,
    const std::vector<std::shared_ptr<Acts::Surface>>& lSurfaces,
    const std::vector<size_t>& assignToAll,
    const Acts::Experimental::ProtoBinning& binning) {
  // The surface candidate updator & a generator for polyhedrons
  Acts::Experimental::SurfaceCandidatesUpdator sfCandidates;
  Acts::Experimental::detail::PolyhedronReferenceGenerator rGenerator;
  // Indexed Surface generator for this case
  Acts::Experimental::detail::IndexedSurfacesGenerator<decltype(lSurfaces)> isg{
      lSurfaces, assignToAll, {binning.binValue}, {binning.expansion}};
  if (binning.axisType == Acts::detail::AxisType::Equidistant) {
    // Equidistant
    Acts::Experimental::detail::GridAxisGenerators::Eq<aType> aGenerator{
        {binning.edges.front(), binning.edges.back()}, binning.bins()};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  } else {
    // Variable
    Acts::Experimental::detail::GridAxisGenerators::Var<aType> aGenerator{
        binning.edges};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  }
  return sfCandidates;
}

/// Helper for 2-dimensional generators
///
/// @tparam aType is the axis boundary type, axis a: closed or bound
/// @tparam bType is the axis boundary type, axis b: closed or bound
///
/// @param gctx the geometry context
/// @param lSurfaces the surfaces of the layer
/// @param assignToAll the indices assigned to all
/// @param aBinning the binning struct of axis a
/// @param bBinning the binning struct of axis b
///
/// @return a configured surface candidate updators
template <Acts::detail::AxisBoundaryType aType,
          Acts::detail::AxisBoundaryType bType>
Acts::Experimental::SurfaceCandidatesUpdator createUpdator(
    const Acts::GeometryContext& gctx,
    const std::vector<std::shared_ptr<Acts::Surface>>& lSurfaces,
    const std::vector<size_t>& assignToAll,
    const Acts::Experimental::ProtoBinning& aBinning,
    const Acts::Experimental::ProtoBinning& bBinning) {
  // The surface candidate updator & a generator for polyhedrons
  Acts::Experimental::SurfaceCandidatesUpdator sfCandidates;
  Acts::Experimental::detail::PolyhedronReferenceGenerator rGenerator;
  // Indexed Surface generator for this case
  Acts::Experimental::detail::IndexedSurfacesGenerator<decltype(lSurfaces)> isg{
      lSurfaces,
      assignToAll,
      {aBinning.binValue, bBinning.binValue},
      {aBinning.expansion, bBinning.expansion}};
  // Run through the cases
  if (aBinning.axisType == Acts::detail::AxisType::Equidistant and
      bBinning.axisType == Acts::detail::AxisType::Equidistant) {
    // Equidistant-Equidistant
    Acts::Experimental::detail::GridAxisGenerators::EqEq<aType, bType>
        aGenerator{{aBinning.edges.front(), aBinning.edges.back()},
                   aBinning.bins(),
                   {bBinning.edges.front(), bBinning.edges.back()},
                   bBinning.bins()};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  } else if (bBinning.axisType == Acts::detail::AxisType::Equidistant) {
    // Variable-Equidistant
    Acts::Experimental::detail::GridAxisGenerators::VarEq<aType, bType>
        aGenerator{aBinning.edges,
                   {bBinning.edges.front(), bBinning.edges.back()},
                   bBinning.bins()};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  } else if (aBinning.axisType == Acts::detail::AxisType::Equidistant) {
    // Equidistant-Variable
    Acts::Experimental::detail::GridAxisGenerators::EqVar<aType, bType>
        aGenerator{{aBinning.edges.front(), aBinning.edges.back()},
                   aBinning.bins(),
                   bBinning.edges};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  } else {
    // Variable-Variable
    Acts::Experimental::detail::GridAxisGenerators::VarVar<aType, bType>
        aGenerator{aBinning.edges, bBinning.edges};
    sfCandidates = isg(gctx, aGenerator, rGenerator);
  }
  // Return the candidates
  return sfCandidates;
}

}  // namespace

Acts::Experimental::LayerStructureBuilder::LayerStructureBuilder(
    const Acts::Experimental::LayerStructureBuilder::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : IInternalStructureBuilder(), m_cfg(cfg), m_logger(std::move(logger)) {
  if (m_cfg.surfacesProvider == nullptr) {
    throw std::invalid_argument(
        "LayerStructureBuilder: surfaces provider is nullptr.");
  }
}

Acts::Experimental::InternalStructure
Acts::Experimental::LayerStructureBuilder::construct(
    const Acts::GeometryContext& gctx) const {
  // Trivialities first: internal volumes
  std::vector<std::shared_ptr<DetectorVolume>> internalVolumes = {};
  DetectorVolumeUpdator internalVolumeUpdator = tryNoVolumes();

  // Print the auxilliary information
  if (not m_cfg.auxilliary.empty()) {
    ACTS_DEBUG(m_cfg.auxilliary);
  }

  // Retrieve the layer surfaces
  SurfaceCandidatesUpdator internalCandidatesUpdator;
  auto internalSurfaces = m_cfg.surfacesProvider->surfaces(gctx);
  ACTS_DEBUG("Building internal layer structure from "
             << internalSurfaces.size() << " provided surfaces.");

  // Check whether support structure is scheduled to be built, and if so
  // collect those that should be assigned to all bins
  std::vector<size_t> assignToAll = {};
  if (not m_cfg.supports.empty()) {
    ACTS_DEBUG("Adding " << m_cfg.supports.size() << " support structures.")
    // The surface candidate updator
    for (const auto& support : m_cfg.supports) {
      // Throw an excpetion is misconfigured
      if (support.type == Surface::SurfaceType::Other) {
        throw std::invalid_argument(
            "LayerStructureBuilder: support surface type not specified.");
      }
      ACTS_VERBOSE("- Build support of type '"
                   << Acts::Surface::s_surfaceTypeNames[support.type] << "'.");
      if (support.splits > 1u) {
        ACTS_VERBOSE("  Support surface is modelled with " << support.splits
                                                           << " planes.");
      }
      // To correctly attach the support structures, estimate the extent
      Extent internalExtent;
      for (const auto& s : internalSurfaces) {
        auto sPolyhedron = s->polyhedronRepresentation(gctx, m_cfg.nSegments);
        internalExtent.extend(sPolyhedron.extent(), support.constraints);
      }
      // Use the support bulder helper to add support surfaces
      detail::SupportHelper::addSupport(
          internalSurfaces, assignToAll, internalExtent, support.type,
          support.values, support.transform, support.splits);
    }
  }

  // Create the indexed surface grids
  if (m_cfg.binnings.size() == 1u) {
    ACTS_DEBUG("- 1-dimensional surface binning detected.");
    // Capture the binning
    auto binning = m_cfg.binnings[0u];
    if (binning.boundaryType == Acts::detail::AxisBoundaryType::Closed) {
      internalCandidatesUpdator =
          createUpdator<Acts::detail::AxisBoundaryType::Closed>(
              gctx, internalSurfaces, assignToAll, binning);
    } else {
      internalCandidatesUpdator =
          createUpdator<Acts::detail::AxisBoundaryType::Bound>(
              gctx, internalSurfaces, assignToAll, binning);
    }
  } else if (m_cfg.binnings.size() == 2u) {
    ACTS_DEBUG("- 2-dimensional surface binning detected.");
    // Capture the binnings
    const auto& binning0 = m_cfg.binnings[0u];
    const auto& binning1 = m_cfg.binnings[1u];

    if (binning0.boundaryType == Acts::detail::AxisBoundaryType::Closed) {
      internalCandidatesUpdator =
          createUpdator<Acts::detail::AxisBoundaryType::Closed,
                        Acts::detail::AxisBoundaryType::Bound>(
              gctx, internalSurfaces, assignToAll, binning0, binning1);
    } else if (binning1.boundaryType ==
               Acts::detail::AxisBoundaryType::Closed) {
      internalCandidatesUpdator =
          createUpdator<Acts::detail::AxisBoundaryType::Bound,
                        Acts::detail::AxisBoundaryType::Closed>(
              gctx, internalSurfaces, assignToAll, binning0, binning1);
    } else {
      internalCandidatesUpdator =
          createUpdator<Acts::detail::AxisBoundaryType::Bound,
                        Acts::detail::AxisBoundaryType::Bound>(
              gctx, internalSurfaces, assignToAll, binning0, binning1);
    }
  }
  // Check if everything went ok
  if (not internalCandidatesUpdator.connected()) {
    throw std::runtime_error(
        "LayerStructureBuilder: could not connect surface candidate updator.");
  }

  // Return the internal structure
  return InternalStructure{internalSurfaces, internalVolumes,
                           std::move(internalCandidatesUpdator),
                           std::move(internalVolumeUpdator)};
}
