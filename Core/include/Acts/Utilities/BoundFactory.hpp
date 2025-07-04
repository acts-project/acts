// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <type_traits>
#include <vector>

namespace Acts {
namespace detail {

/// @brief Concept to define the minimal requirements on the bounds to apply the deduplication mechanism
///        using the BoundFactory.
template <typename BoundsType_t>
concept ComparableBoundConcept = requires(const BoundsType_t& bounds) {
  /// @brief Getter function to distinguish the various bound types (e.g box vs. cylinder)
  { bounds.type() };
  /// @brief Getter function returning all defining parameters of the bounds as an std::vector
  { bounds.values() } -> std::same_as<std::vector<double>>;
};
}  // namespace detail

///   @brief Factory helper class to construct volume or surface bounds where the constructed bounds
///          are cached inside the factory and if the same bound parameters are
///          requested at a later stage the factory automatically returns the
///          cached bounds. This provides a simple sharing mechanism of the same
///          bounds across multiple surfaces / volumes
template <detail::ComparableBoundConcept BoundsType_t>
class BoundFactory {
 public:
  /// @brief Empty default constructor
  BoundFactory() = default;
  /// @brief Delete the copy constructor
  BoundFactory(const BoundFactory& other) = delete;
  /// @brief Delete copy assignment
  BoundFactory& operator=(const BoundFactory& other) = delete;
  /// @brief Pass externally constructed bounds to the factory and run the deduplication
  ///        mechanism on them
  /// @tparam BoundsImpl_t: Template specification of the bounds to deduplicate
  /// @param bounds: Pointer to the bounds to deduplicated
  /// @return Pointer to an equivalent bound object
  template <typename BoundsImpl_t>
  std::shared_ptr<BoundsImpl_t> insert(
      const std::shared_ptr<BoundsImpl_t>& bounds)
    requires(std::is_base_of_v<BoundsType_t, BoundsImpl_t>)
  {
    assert(bounds);
    return std::dynamic_pointer_cast<BoundsImpl_t>(
        *m_boundSet.insert(bounds).first);
  }
  /// @brief Factory method to construct new bounds from the passed arguments
  /// @tparam BoundsImpl_t: Explicit template specification of the bounds to construct
  /// @param args: List of defining bound parameters
  /// @return A pointer to the newly constructed bounds or to an already existing
  ///          equivalent bound object
  template <typename BoundsImpl_t, class... argList>
  std::shared_ptr<BoundsImpl_t> makeBounds(argList... args)
    requires(std::is_base_of_v<BoundsType_t, BoundsImpl_t>)
  {
    return insert(std::make_shared<BoundsImpl_t>(args...));
  }
  /// @brief Returns the number of cached objects inside the bound factory
  std::size_t size() const { return m_boundSet.size(); }

 private:
  ///  @brief Helper struct to actually perform the deduplication of the passed bound objects.
  ///          For a pair of two bound pointers, the struct defines the <
  ///          operator which is then exploited by a std::set to deduplicate
  ///          equivalent bounds.
  struct BoundComparator {
   public:
    /// @brief Implementation of the comparison operator between two bound pointers
    /// @param a: First bound pointer in the comparison
    /// @param b: Second bound pointer in the comparison
    bool operator()(const std::shared_ptr<BoundsType_t>& a,
                    const std::shared_ptr<BoundsType_t>& b) const {
      /// If we deal with two fundamental different bound sets, then just
      /// cast the type to int and return the comparison
      if (a->type() != b->type()) {
        return static_cast<int>(a->type()) < static_cast<int>(b->type());
      }
      const std::vector<double> avalues{a->values()};
      const std::vector<double> bvalues{b->values()};
      /// In case of polygon shaped bounds, the vectors may differ
      if (avalues.size() != bvalues.size()) {
        return avalues.size() < bvalues.size();
      }
      ///  Loop over the defining parameters of the two bounds and compare them
      ///  pairwise. If a difference is spotted, then use it to assign the <
      ///  ordering
      return std::ranges::lexicographical_compare(
          avalues, bvalues,
          [](double parA, double parB) { return parA < parB; });
    }
  };
  std::set<std::shared_ptr<BoundsType_t>, BoundComparator> m_boundSet{};
};

/// @brief Abrivation for a factory to construct surface bounds
using SurfaceBoundFactory = BoundFactory<SurfaceBounds>;
/// @brief Abrivation for a factory to construct volume bounds
using VolumeBoundFactory = BoundFactory<VolumeBounds>;

}  // namespace Acts
