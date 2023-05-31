// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Detector/detail/IndexedGridFiller.hpp"
#include "Acts/Detector/detail/IndexedSurfacesGenerator.hpp"
#include "Acts/Detector/detail/ReferenceGenerators.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/KDTree.hpp"

#include <array>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace Acts {

namespace Experimental {

namespace detail {

/// @brief A wrapper class around a KDTree of surfaces
///
/// It also deals with the conversion from global query to
/// KDTree lookup positions
///
template <size_t kDIM = 2u, size_t bSize = 100u,
          typename reference_generator = PolyhedronReferenceGenerator>
class KdtSurfaces {
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
  KdtSurfaces(const GeometryContext& gctx,
              const std::vector<std::shared_ptr<Surface>>& surfaces,
              const std::array<BinningValue, kDIM>& casts,
              const reference_generator& rgen = PolyhedronReferenceGenerator{})
      : m_kdt(nullptr), m_casts(casts), m_rGenerator(rgen) {
    // Simple check if the dimension is correct
    if (kDIM == 0u) {
      throw std::invalid_argument(
          "KdtSurfaces: dimension and/or cast rules are incorrect.");
    }
    // Fill the tree from surfaces
    std::vector<Entry> kdtEntries;
    kdtEntries.reserve(surfaces.size());
    for (auto& s : surfaces) {
      // Generate the references and the center of gravity from it
      const auto references = m_rGenerator.references(gctx, *s);
      const auto ref = cog(references);
      //  Now cast into the correct fill position
      std::array<ActsScalar, kDIM> fill = {};
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
  /// @note no checking on positions.empty() is done as the
  /// positions are to be provided by a generator which itself
  /// is tested for consistency
  ///
  /// @return the center of gravity
  Vector3 cog(const std::vector<Vector3>& positions) const {
    // Build the center of gravity of the n positions
    Vector3 c(0., 0., 0.);
    ActsScalar weight = 1. / positions.size();
    std::for_each(positions.begin(), positions.end(),
                  [&](const auto& p) { c += weight * p; });
    return c;
  }
};

/// @brief Callable struct wrapper around the KDT surface structure
///
/// This allows to create small region based callable structs at
/// configuration level that are then connected to an InternalStructureBuilder
template <size_t kDIM = 2u, size_t bSize = 100u,
          typename reference_generator = PolyhedronReferenceGenerator>
struct KdtSurfacesProvider {
  /// The prefilled surfaces in a KD tree structure, it is generally shared
  /// amongst different providers
  std::shared_ptr<KdtSurfaces<kDIM, bSize, reference_generator>> kdt = nullptr;
  /// The query region
  Extent region;
  /// The call operator that provides the function call
  std::vector<std::shared_ptr<Surface>> operator()() const {
    return kdt->surfaces(region);
  }
};

}  // namespace detail

}  // namespace Experimental

}  // namespace Acts
