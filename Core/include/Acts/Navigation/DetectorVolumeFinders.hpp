// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/detail/IndexedLookupHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationDelegateHelpers.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/NavigationStateFillers.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <stdexcept>

namespace Acts {
namespace Experimental {

using VariableBoundAxis =
    Acts::detail::Axis<Acts::detail::AxisType::Variable,
                       Acts::detail::AxisBoundaryType::Bound>;
using VariableBoundIndexGrid1 =
    Acts::detail::Grid<std::size_t, VariableBoundAxis>;

struct RootVolumeFinder final : public IDetectorVolumeFinder {
  inline const DetectorVolume* find(const GeometryContext& gctx,
                                    const NavigationState& nState) const final {
    auto currentDetector = nState.currentDetector;

    if (currentDetector == nullptr) {
      throw std::runtime_error(
          "DetectorVolumeFinders: no detector set to navigation state.");
    }

    const auto& volumes = currentDetector->rootVolumes();
    for (const auto v : volumes) {
      if (v->inside(gctx, nState.position)) {
        if (v->detectorVolumeFinder()) {
          return v->detectorVolumeFinder()(gctx, nState);
        }
        return v;
      }
    }

    return nullptr;
  }
};

struct TrialAndErrorVolumeFinder final : public IDetectorVolumeFinder {
  inline const DetectorVolume* find(const GeometryContext& gctx,
                                    const NavigationState& nState) const final {
    auto currentVolume = nState.currentVolume;

    if (currentVolume == nullptr) {
      throw std::runtime_error(
          "DetectorVolumeFinders: no volume set to navigation state.");
    }

    if (!currentVolume->inside(gctx, nState.position)) {
      return nullptr;
    }

    const auto& volumes = currentVolume->volumes();
    for (const auto v : volumes) {
      if (v->inside(gctx, nState.position)) {
        if (v->detectorVolumeFinder()) {
          return v->detectorVolumeFinder()(gctx, nState);
        }
        return v;
      }
    }

    return currentVolume;
  }
};

/// @brief The end of world sets the volume pointer of the
/// navigation state to nullptr, usually indicates the end of
/// the known world, hence the name
struct EndOfWorldVolume final : public IDetectorVolumeFinder {
  /// @brief a null volume link - explicitely
  ///
  /// @note the method parameters are ignored
  inline const DetectorVolume* find(
      const GeometryContext& /*gctx*/,
      const NavigationState& /*nState*/) const final {
    return nullptr;
  }
};

/// @brief Single volume updator, it sets the current navigation
/// volume to the volume in question
///
struct SingleDetectorVolume final : public IDetectorVolumeFinder {
  const DetectorVolume* dVolume = nullptr;

  /// @brief Allowed constructor
  /// @param sVolume the volume to which it points
  SingleDetectorVolume(const DetectorVolume* sVolume) noexcept(false)
      : dVolume(sVolume) {
    if (sVolume == nullptr) {
      throw std::invalid_argument(
          "DetectorVolumeUpdators: nullptr provided, use EndOfWorld instead.");
    }
  }

  /// @brief a null volume link - explicitely
  ///
  /// @note the method parameters are ignored
  ///
  inline const DetectorVolume* find(
      const GeometryContext& /*gctx*/,
      const NavigationState& /*nState*/) const final {
    return dVolume;
  }
};

/// @brief This is used for volumes that are indexed in a bound
/// 1-dimensional grid, e.g. a z-spaced array, or an r-spaced array
/// of volumes.
///
struct BoundVolumesGrid1 final : public IDetectorVolumeFinder {
  using IndexedUpdator =
      IndexedLookupHelper<VariableBoundIndexGrid1, DetectorVolumesCollection,
                          DetectorVolumeFiller>;

  // The indexed updator
  IndexedUpdator indexedUpdator;

  /// Allowed constructor with explicit arguments
  ///
  /// @param gBoundaries the grid boundaries
  /// @param bValue the binning value
  /// @param cVolumes the contained volumes
  /// @param bTransform is the optional transform
  BoundVolumesGrid1(
      const std::vector<ActsScalar>& gBoundaries, BinningValue bValue,
      const std::vector<const DetectorVolume*>& cVolumes,
      const Transform3& bTransform = Transform3::Identity()) noexcept(false)
      : indexedUpdator(IndexedUpdator(VariableBoundIndexGrid1(std::make_tuple(
                                          VariableBoundAxis(gBoundaries))),
                                      {bValue}, bTransform)) {
    indexedUpdator.extractor.dVolumes = cVolumes;

    if (gBoundaries.size() != cVolumes.size() + 1u) {
      throw std::invalid_argument(
          "DetectorVolumeUpdators: mismatching boundaries and volume numbers");
    }
    // Initialize the grid entries
    for (std::size_t ib = 1u; ib < gBoundaries.size(); ++ib) {
      indexedUpdator.grid.at(ib) = ib - 1;
    }
  }

  /// @brief This updator relies on an 1D single index grid
  ///
  /// @param gctx the geometry context
  /// @param nState [in,out] the navigation state to be updated
  inline const DetectorVolume* find(const GeometryContext& gctx,
                                    const NavigationState& nState) const final {
    return indexedUpdator.lookup(gctx, nState);
  }
};

template <typename indexed_lookup_t>
struct GenericIndexedVolumeFinder : public IDetectorVolumeFinder {
  indexed_lookup_t indexedLookup;

  GenericIndexedVolumeFinder(indexed_lookup_t&& _indexedLookup)
      : indexedLookup(std::move(_indexedLookup)) {}

  inline const DetectorVolume* find(const GeometryContext& gctx,
                                    const NavigationState& nState) const final {
    // Extract the index grid entry
    return indexedLookup.lookup(gctx, nState);
  }
};

template <typename grid_type>
using IndexedVolumeLookupHelper =
    IndexedLookupHelper<grid_type, IndexedDetectorVolumeExtractor,
                        DetectorVolumeFiller>;

/// @brief  An indexed volume implementation access
///
/// @tparam grid_type is the grid type used for this
template <typename grid_type>
using IndexedDetectorVolumeFinder =
    GenericIndexedVolumeFinder<IndexedVolumeLookupHelper<grid_type>>;

template <typename grid_type>
inline static IndexedDetectorVolumeFinder<grid_type>
makeIndexedDetectorVolumeFinder(
    grid_type igrid, const std::array<BinningValue, grid_type::DIM>& icasts,
    const Transform3& itr = Transform3::Identity()) {
  auto lookupHelper = IndexedVolumeLookupHelper<grid_type>(
      std::forward<grid_type>(igrid), icasts, itr);
  return IndexedDetectorVolumeFinder<grid_type>(std::move(lookupHelper));
}

}  // namespace Experimental
}  // namespace Acts
