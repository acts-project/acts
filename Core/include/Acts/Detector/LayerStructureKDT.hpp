// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/GridAxisGenerators.hpp"
#include "Acts/Detector/IndexedGridFiller.hpp"
#include "Acts/Detector/IndexedSurfacesGenerator.hpp"
#include "Acts/Detector/SupportBuilder.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/KDTree.hpp"

#include <array>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace Acts {

namespace Experimental {

/// @brief A wrapper class around a KDTree of surfaces
///
/// It also deals with the conversion from global query to
/// KDTree lookup positions
///
template <size_t kDIM = 2u, size_t bSize = 100u,
          typename reference_generator = PolyhedronReferenceGenerator>
class SurfacesKDT {
 public:
  /// Broadcast the surface KDT type
  using KDTS =
      KDTree<kDIM, std::shared_ptr<Surface>, ActsScalar, std::array, bSize>;

  /// Broadcast the query definition
  using Query = std::array<ActsScalar, kDIM>;

  /// Broadcast the entry
  using Entry = std::pair<Query, std::shared_ptr<Surface>>;

  /// Constructor from a vector of surfaces
  ///
  /// @param gctx the geometry context of this call
  /// @param surfaces the surfaces to be filled into the tree
  /// @param casts the cast list from global position into kdtree local
  /// @param rgen the reference point generator
  ///
  SurfacesKDT(const GeometryContext& gctx,
              const std::vector<std::shared_ptr<Surface>>& surfaces,
              const std::array<BinningValue, kDIM>& casts,
              const reference_generator& rgen = PolyhedronReferenceGenerator{})
      : m_kdt(nullptr), m_casts(casts), m_rGenerator(rgen) {
    // Simple check if the dimension is correct
    if (kDIM == 0u or kDIM != m_casts.size()) {
      throw std::invalid_argument(
          "SurfacesKDT: dimension and/or cast rules are incorrect.");
    }
    // Fill the tree from surfaces
    std::vector<Entry> kdtEntries;
    kdtEntries.reserve(surfaces.size());
    for (auto& s : surfaces) {
      // Generate the references and the center of gravity from it
      const auto references = m_rGenerator.references(gctx, *s);
      const auto ref = cog(references);
      //  Now cast into the correct fill position
      std::array<ActsScalar, kDIM> fill;
      fillCasts(ref, fill, std::make_integer_sequence<std::size_t, kDIM>{});
      kdtEntries.push_back({fill, s});
    }
    // Create the KDTree
    m_kdt = std::make_unique<KDTS>(std::move(kdtEntries));
  }

  /// Query with a Range object
  ///
  /// @param range is the range to be queried
  ///
  /// @return the matching surfaces from the KDT structure
  std::vector<std::shared_ptr<Surface>> surfaces(
      const RangeXD<kDIM, ActsScalar>& range) const {
    // Strip the surfaces
    std::vector<std::shared_ptr<Surface>> surfacePtrs;
    auto surfaceQuery = m_kdt->rangeSearchWithKey(range);
    std::for_each(surfaceQuery.begin(), surfaceQuery.end(),
                  [&](auto& s) { surfacePtrs.push_back(s.second); });
    return surfacePtrs;
  }

  /// Query with an Extent object
  ///
  /// @param extent is the range Extent to be queried
  ///
  /// @return the matching surfaces fpulled from the KDT structure
  std::vector<std::shared_ptr<Surface>> surfaces(const Extent& extent) const {
    RangeXD<kDIM, ActsScalar> qRange;
    for (auto [ibv, v] : enumerate(m_casts)) {
      qRange[ibv] = extent.range(v);
    }
    return surfaces(qRange);
  }

 private:
  /// The KDTree as single source for the surfaces (maybe shared)
  std::unique_ptr<KDTS> m_kdt = nullptr;

  /// Cast values that turn a global position to lookup position
  std::array<BinningValue, kDIM> m_casts = {};

  /// Helper to generate refernce points for filling
  reference_generator m_rGenerator;

  /// Unroll the cast loop
  /// @param position is the position of the update call
  /// @param a is the array to be filled
  template <typename Array, std::size_t... idx>
  void fillCasts(const Vector3& position, Array& a,
                 std::index_sequence<idx...> /*indices*/) const {
    ((a[idx] = VectorHelpers::cast(position, m_casts[idx])), ...);
  }

  /// Helper method to calculate the center of gravity
  ///
  /// @param positions are the reference positions that go in
  /// @note will do nothing if vector size is equal to 1
  ///
  /// @return the center of gravity
  Vector3 cog(const std::vector<Vector3>& positions) const {
    // Throw exception if misconfigured, or generation did not work
    if (positions.empty()) {
      throw std::runtime_error("SurfacesKDT: reference points empty.");
    }
    // Return the single one you have
    if (positions.size() == 1u) {
      return positions[0u];
    }
    // Build the center of gravity of the n positions
    Vector3 c(0., 0., 0.);
    ActsScalar weight = 1. / positions.size();
    std::for_each(positions.begin(), positions.end(),
                  [&](const auto& p) { c += weight * p; });
    return c;
  }
};

/// @brief Class to generate indexed surface structures
///
/// It acts as a translator from geometry configuration input,
/// to the block building, and hence dispatches some configuration
/// types into actual c++ gnerators for Axes, IndexGrid and eventually
/// SurfaceCandidateUpdators and surfaces.
///
/// To keep the code simple, this deals only with common
/// layer structures, namely:
///
/// - cylindrical structures
///   binning: in z, phi, z - phi
///   support: in delat(r), with envelope extension in z-,z+, phi-, phi+
///
///   @note if a layer transform is given and the identity transform,
///   then the cylinder is forced to be concentric at max z extent
///
/// - disk like structures
///   binning: in r, phi, r - phi
///   support: in delta(z), with envelope extension in r-,r+, phi-, phi+
///
/// - plane like structures
///   binning: in x, y, x - y
///   support: in delta(local z), envelope extension in x-, x+, y-, y+
///
/// @note this helper is only for detectors of regular orientation, i.e.
/// disc and cylinder surfaces are assumed to be aligned with the z axis
/// per surface frame definition
///
/// @note Axis generators are restricted to Closed and Bound, which are
/// the most inclusinve options one can create for a navigation stream
///
/// @tparam kDIM the dimension of the KDTree (and thus query points)
///
/// @note a shared KDTree of surfaces is used for
template <size_t kDIM = 2u>
struct LayerStructureKDT {
  /// Shorthand of surface collection
  using StructureSurfaces = std::vector<std::shared_ptr<Surface>>;

  // Support parameter defintions
  struct Support {
    /// Define whether you want to build support structures
    std::array<ActsScalar, 5u> values = {};
    /// The (optional) layer transform
    std::optional<Transform3> transform = std::nullopt;
    /// Potential splits into planar approximations
    unsigned int splits = 1u;
  };

  // The surface binning
  struct Binning {
    BinningData data;
    size_t expansion = 0u;
  };

  /// The shared KDTree
  std::shared_ptr<SurfacesKDT<kDIM>> surfacesKDT;
  /// The Extent that represents this Layer
  Extent layerExtent;
  /// The representative type
  Surface::SurfaceType representation = Surface::SurfaceType::Other;
  /// The layer support
  std::vector<Support> layerSupports = {};
  /// The surface binning
  std::vector<Binning> surfaceBinning = {};

  /// @brief  This creates the internal structure of a
  /// DetectorVolume that represents a layer like object
  ///
  /// It's a preconfigured struct that is called in the
  /// DetectorBuilding process
  ///
  /// @note on demand it can create extra surfaces to
  /// represent passive structures
  ///
  /// @param gctx The geometry contrext of this call
  ///
  /// @return a tuple of the contained
  std::tuple<StructureSurfaces, SurfaceCandidatesUpdator> create(
      const GeometryContext& gctx) const {
    // Surface binning is not present - misconfigured
    if (surfaceBinning.empty()) {
      throw std::invalid_argument(
          "LayerStructureKDT: No surface binning provided.");
    }
    // Get the surfaces first from the KDTree
    StructureSurfaces lSurfaces = surfacesKDT->surfaces(layerExtent);
    std::vector<size_t> assignToAll = {};
    // The surface generator
    SurfaceCandidatesUpdator sfCandidates;
    // Add support if needed
    for (const auto& ls : layerSupports) {
      // Use the support bulder helper to add support surfaces
      SupportBuilder::addSupport(lSurfaces, assignToAll, layerExtent,
                                 representation, ls.values, ls.transform,
                                 ls.splits);
    }
    // Create the correct axis generators
    if (surfaceBinning.size() == 1u) {
      // Capture the binning
      auto binning = surfaceBinning[0u];
      if (binning.data.option == closed) {
        sfCandidates = createUpdator<detail::AxisBoundaryType::Closed>(
            gctx, lSurfaces, assignToAll, surfaceBinning[0u]);
      } else {
        sfCandidates = createUpdator<detail::AxisBoundaryType::Bound>(
            gctx, lSurfaces, assignToAll, surfaceBinning[0u]);
      }
    } else if (surfaceBinning.size() == 2u) {
      // Capture the binnings
      const auto& binning0 = surfaceBinning[0u];
      const auto& binning1 = surfaceBinning[1u];
      if (binning0.data.option == closed) {
        sfCandidates = createUpdator<detail::AxisBoundaryType::Closed,
                                     detail::AxisBoundaryType::Bound>(
            gctx, lSurfaces, assignToAll, binning0, binning1);
      } else if (binning1.data.option == closed) {
        sfCandidates = createUpdator<detail::AxisBoundaryType::Bound,
                                     detail::AxisBoundaryType::Closed>(
            gctx, lSurfaces, assignToAll, binning0, binning1);
      } else {
        sfCandidates = createUpdator<detail::AxisBoundaryType::Bound,
                                     detail::AxisBoundaryType::Bound>(
            gctx, lSurfaces, assignToAll, binning0, binning1);
      }
    }
    // Check if everything went ok
    if (not sfCandidates.connected()) {
      throw std::runtime_error(
          "LayerStructureKDT: could not connect surface candidate updator.");
    }
    // Rerturn to the caller to build a DetectorVolume
    return {lSurfaces, std::move(sfCandidates)};
  }

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
  template <detail::AxisBoundaryType aType>
  SurfaceCandidatesUpdator createUpdator(
      const GeometryContext& gctx, const StructureSurfaces& lSurfaces,
      const std::vector<size_t>& assignToAll,
      const LayerStructureKDT<kDIM>::Binning& binning) const {
    // The surface candidate updator & a generator for polyhedrons
    SurfaceCandidatesUpdator sfCandidates;
    PolyhedronReferenceGenerator rGenerator;
    // Indexed Surface generator for this case
    IndexedSurfacesGenerator<StructureSurfaces> isg{
        lSurfaces, assignToAll, {binning.data.binvalue}, {binning.expansion}};
    if (binning.data.type == equidistant) {
      // Equidistant
      GridAxisGenerators::Eq<aType> aGenerator{{-M_PI, M_PI},
                                               binning.data.bins()};
      sfCandidates = isg(gctx, aGenerator, rGenerator);
    } else {
      std::vector<ActsScalar> edges = {binning.data.boundaries().begin(),
                                       binning.data.boundaries().end()};
      // Variable
      GridAxisGenerators::Var<aType> aGenerator{edges};
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
  template <detail::AxisBoundaryType aType, detail::AxisBoundaryType bType>
  SurfaceCandidatesUpdator createUpdator(
      const GeometryContext& gctx, const StructureSurfaces& lSurfaces,
      const std::vector<size_t>& assignToAll,
      const LayerStructureKDT<kDIM>::Binning& aBinning,
      const LayerStructureKDT<kDIM>::Binning& bBinning) const {
    // The surface candidate updator & a generator for polyhedrons
    SurfaceCandidatesUpdator sfCandidates;
    PolyhedronReferenceGenerator rGenerator;
    // Indexed Surface generator for this case
    IndexedSurfacesGenerator<StructureSurfaces> isg{
        lSurfaces,
        assignToAll,
        {aBinning.data.binvalue, bBinning.data.binvalue},
        {aBinning.expansion, bBinning.expansion}};
    // Run through the cases
    if (aBinning.data.type == equidistant and
        bBinning.data.type == equidistant) {
      // Equidistant-Equidistant
      GridAxisGenerators::EqEq<aType, bType> aGenerator{
          {aBinning.data.min, aBinning.data.max},
          aBinning.data.bins(),
          {bBinning.data.min, bBinning.data.max},
          bBinning.data.bins()};
      sfCandidates = isg(gctx, aGenerator, rGenerator);
    } else if (bBinning.data.type == equidistant) {
      // Variable-Equidistant
      std::vector<ActsScalar> edges0 = {bBinning.data.boundaries().begin(),
                                        bBinning.data.boundaries().end()};
      GridAxisGenerators::VarEq<aType, bType> aGenerator{
          edges0, {bBinning.data.min, bBinning.data.max}, bBinning.data.bins()};
      sfCandidates = isg(gctx, aGenerator, rGenerator);
    } else if (aBinning.data.type == equidistant) {
      // Equidistant-Variable
      std::vector<ActsScalar> edges1 = {bBinning.data.boundaries().begin(),
                                        bBinning.data.boundaries().end()};
      GridAxisGenerators::EqVar<aType, bType> aGenerator{
          {aBinning.data.min, aBinning.data.max}, aBinning.data.bins(), edges1};
      sfCandidates = isg(gctx, aGenerator, rGenerator);
    } else {
      // Variable-Variable
      std::vector<ActsScalar> edges0 = {aBinning.data.boundaries().begin(),
                                        aBinning.data.boundaries().end()};
      std::vector<ActsScalar> edges1 = {bBinning.data.boundaries().begin(),
                                        bBinning.data.boundaries().end()};
      GridAxisGenerators::VarVar<aType, bType> aGenerator{edges0, edges1};
      sfCandidates = isg(gctx, aGenerator, rGenerator);
    }
    // return the candidates
    return sfCandidates;
  }
};

}  // namespace Experimental

}  // namespace Acts
