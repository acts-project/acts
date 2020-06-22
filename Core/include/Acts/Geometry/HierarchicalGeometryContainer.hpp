// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <vector>

#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/detail/DefaultGeometryIdGetter.hpp"

namespace Acts {

/// Store homogeneous elements mapped into the geometry hierarchy.
///
/// @tparam value_t stored element type
/// @tparam identifier_getter_t functor to access element geometry identifier.
/// must be default-constructible and have no internal state.
///
/// The core functionality is to find an equivalent element for a given
/// identifier via
///
///     auto it = container.find(GeometryID(...));
///     if (it != container.end()) {
///         ...
///     }
///
/// Trailing zero levels of stored geometry identifiers are used as
/// broadcast values to refer to higher-level objects within the geometry, e.g.
/// a geometry identifier with vanishing approach and sensitive index identifies
/// a layer. An entry will all geometry identifier levels set to zero acts as
/// the global default value.
///
/// The container also supports range-based iteration over all stored elements
///
///     for (const auto& element : container) {
///         ...
///     }
///
/// and index-based access to stored elements and associated geometry
/// identifiers
///
///     GeometryID id3 = container.idAt(3);
///     const auto& element4 = container.valueAt(4);
///
/// @note No guarantees are given for the element order when using range-based
///   or index-based access. Any apparent ordering must be considered an
///   implementation detail and might change.
///
/// Adding elements is potentially expensive as the internal lookup structure
/// must be updated. In addition, modifying an element in-place could change its
/// identifier which would also break the lookup. Thus, the container can not be
/// modified after construction to prevent misuse.
template <typename value_t,
          typename identifier_getter_t = detail::DefaultGeometryIdGetter>
class HierarchicalGeometryContainer {
 public:
  using Iterator = typename std::vector<value_t>::const_iterator;
  using Size = typename std::vector<value_t>::size_type;
  using Value = value_t;

  /// Construct the container from the given elements.
  ///
  /// @param values input elements (must be unique with respect to identifier)
  HierarchicalGeometryContainer(std::vector<Value>&& values);
  /// Construct the container from an initializer list.
  ///
  /// @param values input initializer list
  HierarchicalGeometryContainer(std::initializer_list<Value> values);

  // defaulted constructors and assignment operators
  HierarchicalGeometryContainer() = default;
  HierarchicalGeometryContainer(const HierarchicalGeometryContainer&) = default;
  HierarchicalGeometryContainer(HierarchicalGeometryContainer&&) = default;
  ~HierarchicalGeometryContainer() = default;
  HierarchicalGeometryContainer& operator=(
      const HierarchicalGeometryContainer&) = default;
  HierarchicalGeometryContainer& operator=(HierarchicalGeometryContainer&&) =
      default;

  /// Return an iterator pointing to the beginning of the stored elements.
  Iterator begin() const { return m_values.begin(); }
  /// Return an iterator pointing to the end of the stored elements.
  Iterator end() const { return m_values.end(); }
  /// Check if any elements are stored.
  bool empty() const { return m_values.empty(); }
  /// Return the number of stored elements.
  Size size() const { return m_values.size(); }

  /// Access the geometry identifier for the i-th element with bounds check.
  ///
  /// @throws std::out_of_range for invalid indices
  GeometryID idAt(Size index) const { return m_ids.at(index); }
  /// Access the i-th element in the container with bounds check.
  ///
  /// @throws std::out_of_range for invalid indices
  const Value& valueAt(Size index) const { return m_values.at(index); }

  /// Find the most specific element for a given geometry identifier.
  ///
  /// This can be either the element matching exactly to the given geometry id,
  /// if it exists, or the element from the next available higher level within
  /// the geometry hierachy.
  ///
  /// @param id geometry identifier for which information is requested
  /// @retval iterator to an existing element
  /// @retval `.end()` iterator if no matching element exists
  Iterator find(GeometryID id) const;

 private:
  // NOTE this class assumes that it knows the ordering of the levels within
  //      the geometry id. if the geometry id changes, this code has to be
  //      adapted too. the asserts ensure that such a change is caught.
  static_assert(GeometryID().setVolume(1).value() <
                    GeometryID().setVolume(1).setBoundary(1).value(),
                "Incompatible GeometryID hierarchy");
  static_assert(GeometryID().setBoundary(1).value() <
                    GeometryID().setBoundary(1).setLayer(1).value(),
                "Incompatible GeometryID hierarchy");
  static_assert(GeometryID().setLayer(1).value() <
                    GeometryID().setLayer(1).setApproach(1).value(),
                "Incompatible GeometryID hierarchy");
  static_assert(GeometryID().setApproach(1).value() <
                    GeometryID().setApproach(1).setSensitive(1).value(),
                "Incompatible GeometryID hierarchy");

  using Identifier = GeometryID::Value;
  using IdentifierGetter = identifier_getter_t;

  std::vector<Value> m_values;
  // encoded ids for all elements for faster lookup.
  std::vector<Identifier> m_ids;
  // validity bit masks for the ids: which parts to use for comparison
  std::vector<Identifier> m_masks;

  /// Construct a mask where all leading non-zero levels are set.
  static constexpr Identifier makeLeadingLevelsMask(GeometryID id) {
    // construct id from encoded value with all bits set
    auto allSet = GeometryID(~GeometryID::Value(0u));
    // manually iterate over identifier levels starting from the lowest
    if (id.sensitive() != 0u) {
      // all levels are valid; keep all bits set.
      return allSet.value();
    }
    if (id.approach() != 0u) {
      return allSet.setSensitive(0u).value();
    }
    if (id.layer() != 0u) {
      return allSet.setSensitive(0u).setApproach(0u).value();
    }
    if (id.boundary() != 0u) {
      return allSet.setSensitive(0u).setApproach(0u).setLayer(0u).value();
    }
    if (id.volume() != 0u) {
      return allSet.setSensitive(0u)
          .setApproach(0u)
          .setLayer(0u)
          .setBoundary(0u)
          .value();
    }
    // no valid levels; all bits are zero.
    return Identifier(0u);
  }
  /// Construct a mask where only the highest level is set.
  static constexpr Identifier makeHighestLevelMask() {
    return makeLeadingLevelsMask(GeometryID(0u).setVolume(1u));
  }
  /// Compare the two identifiers only within the masked bits.
  static constexpr bool equalWithinMask(Identifier lhs, Identifier rhs,
                                        Identifier mask) {
    return (lhs & mask) == (rhs & mask);
  }

  /// Ensure ordering and uniqueness.
  void sortAndCheckDuplicates() {
    auto idLess = [=](const auto& lhs, const auto& rhs) {
      return IdentifierGetter()(lhs) < IdentifierGetter()(rhs);
    };
    auto idEqual = [](const auto& lhs, const auto& rhs) {
      return IdentifierGetter()(lhs) == IdentifierGetter()(rhs);
    };
    // ensure elements are sorted
    std::sort(m_values.begin(), m_values.end(), idLess);
    // check that all elements are unique
    auto duplicate =
        std::adjacent_find(m_values.begin(), m_values.end(), idEqual);
    if (duplicate != m_values.end()) {
      throw std::invalid_argument("Input elements contain duplicates");
    }
  }
  /// Fill the lookup containers from the current values.
  void fillLookup() {
    m_ids.clear();
    m_ids.reserve(m_values.size());
    m_masks.clear();
    m_masks.reserve(m_values.size());
    IdentifierGetter getId;
    for (const auto& element : m_values) {
      m_ids.push_back(getId(element).value());
      m_masks.push_back(makeLeadingLevelsMask(getId(element).value()));
    }
  }
};

// implementations

template <typename value_t, typename identifier_getter_t>
inline HierarchicalGeometryContainer<value_t, identifier_getter_t>::
    HierarchicalGeometryContainer(std::vector<Value>&& values)
    : m_values(std::move(values)) {
  sortAndCheckDuplicates();
  fillLookup();
}

template <typename value_t, typename identifier_getter_t>
inline HierarchicalGeometryContainer<value_t, identifier_getter_t>::
    HierarchicalGeometryContainer(std::initializer_list<Value> values)
    : HierarchicalGeometryContainer(
          std::vector<Value>(values.begin(), values.end())) {}

template <typename value_t, typename identifier_getter_t>
inline auto HierarchicalGeometryContainer<value_t, identifier_getter_t>::find(
    GeometryID id) const -> Iterator {
  assert((m_values.size() == m_ids.size()) and
         "Inconsistent container state: #elements != #ids");
  assert((m_values.size() == m_masks.size()) and
         "Inconsistent container state: #elements != #masks");

  // we can not search for the element directly since the relevant one
  // might be stored at a higher level. ids for higher levels would always
  // be sorted before the requested id. searching for the first element
  // after the requested ensures that we include the full hierarchy.
  const auto it = std::upper_bound(m_ids.begin(), m_ids.end(), id.value());
  auto i = std::distance(m_ids.begin(), it);

  // now go up the hierarchy to find the first matching element.
  // example: the container stores four identifiers
  //
  //     2|x|x (volume-only)
  //     2|2|1 (volume, layer, and sensitive)
  //     2|3|x (volume and layer)
  //     2|3|4 (volume, layer, and sensitive)
  //
  // where | marks level boundaries. searching for either 2|3|4, 2|3|7, or
  // 2|4|x would first point to 2|3|4 and thus needs to go up the hierarchy.
  while (0 < i) {
    // index always starts before item of interest due to upper bound search
    --i;

    // if the input id does not even match at the highest hierarchy level
    // with the current comparison id, then have reached the end of this
    // hierarchy. having a special check for the highest level avoids an
    // unbounded search window all the way to the beginning of the container for
    // the global default entry.
    if (not equalWithinMask(id.value(), m_ids[i], makeHighestLevelMask())) {
      // check if a global default entry exists
      if (m_ids.front() == Identifier(0u)) {
        return begin();
      } else {
        return end();
      }
    }

    // since the search is going backwards in the sorted container, it
    // progresses from more specific to less specific elements. the first
    // match is automatically the appropriate one.
    if (equalWithinMask(id.value(), m_ids[i], m_masks[i])) {
      return std::next(begin(), i);
    }
  }

  // all options are exhausted and not matching element was found.
  return end();
}

}  // namespace Acts
