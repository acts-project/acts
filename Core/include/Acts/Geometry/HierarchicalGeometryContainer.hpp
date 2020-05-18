// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/detail/DefaultGeometryIdGetter.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace Acts {

/// Store homogeneous elements anchored in the geometry hierarchy.
///
/// @tparam value_t stored element type
///
/// The core functionality is to find an equivalent element for a given
/// identifier via
///
///     auto it = container.find(GeometryID(...));
///     if (it != container.end()) {
///         ...
///     }
///
/// Trailing zero components of stored geometry identifiers are used as
/// broadcast values to refer to higher-level objects within the geometry, e.g.
/// a geometry identifier with vanishing approach and sensitive index refers
/// to the full layer. All stored geometry identifiers must set at least
/// the highest hierarchy level.
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
/// must be updated. In addition, modifying an element in-place can
/// potentially change its identifier which would also break the lookup. Thus,
/// the container can not be modified after construction to prevent misuse.
///
template <typename value_t>
class HierarchicalGeometryContainer {
 public:
  using Iterator = typename std::vector<value_t>::const_iterator;
  using Size = typename std::vector<value_t>::size_type;
  using Value = value_t;

  /// Construct the container from the given elements.
  ///
  /// @param elements input elements (must be unique with respect to identifier)
  /// @param getId functor to retrieve the id of an element
  ///
  /// @tparam id_getter_t Functor type to retrieve the id of an element
  template <typename id_getter_t = detail::DefaultGeometryIdGetter>
  HierarchicalGeometryContainer(std::vector<Value>&& elements,
                                id_getter_t getId = id_getter_t());

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
  Iterator begin() const { return m_elements.begin(); }
  /// Return an iterator pointing to the end of the stored elements.
  Iterator end() const { return m_elements.end(); }
  /// Check if any elements are stored.
  bool empty() const { return m_elements.empty(); }
  /// Return the number of stored elements.
  Size size() const { return m_elements.size(); }

  /// Access the geometry identifier for the i-th element with bounds check.
  ///
  /// @throws std::out_of_range for invalid indices
  GeometryID idAt(Size index) const { return m_ids.at(index); }
  /// Access the i-th element in the container with bounds check.
  ///
  /// @throws std::out_of_range for invalid indices
  const Value& valueAt(Size index) const { return m_elements.at(index); }

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
  using Identifier = GeometryID::Value;

  std::vector<Value> m_elements;
  // encoded ids for all elements for faster lookup and to avoid having
  // to store a method to retrieve the ids from the stored elements.
  std::vector<Identifier> m_ids;
  // validity bit masks for the ids: which parts to use for comparison
  std::vector<Identifier> m_masks;

  /// Construct a mask where all leading non-zero levels are set.
  static constexpr Identifier makeLeadingLevelsMask(GeometryID id) {
    // construct id from encoded value with all bits set
    auto mask = GeometryID(~GeometryID::Value(0u));
    // NOTE this code assumes some knowledge about the ordering of the levels
    //      within the geometry id. if the geometry id changes, this code has
    //      to be adapted too.
    if (id.sensitive() != 0u) {
      // all levels are valid; keep all at one.
      return mask.value();
    }
    if (id.approach() != 0u) {
      return mask.setSensitive(0u).value();
    }
    if (id.layer() != 0u) {
      return mask.setSensitive(0u).setApproach(0u).value();
    }
    if (id.boundary() != 0u) {
      return mask.setSensitive(0u).setApproach(0u).setLayer(0u).value();
    }
    // identifier must always have non-zero volume
    return mask.setSensitive(0u)
        .setApproach(0u)
        .setLayer(0u)
        .setBoundary(0u)
        .value();
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

  /// Fill the lookup containers from the current values.
  ///
  /// @tparam id_getter_t Functor type to retrieve the id of an element
  template <typename id_getter_t>
  void fillLookup(id_getter_t getId) {
    m_ids.clear();
    m_ids.reserve(m_elements.size());
    m_masks.clear();
    m_masks.reserve(m_elements.size());
    for (const auto& element : m_elements) {
      m_ids.push_back(getId(element).value());
      m_masks.push_back(makeLeadingLevelsMask(getId(element).value()));
    }
  }
};

// implementations

template <typename value_t>
template <typename id_getter_t>
inline HierarchicalGeometryContainer<value_t>::HierarchicalGeometryContainer(
    std::vector<Value>&& elements, id_getter_t getId)
    : m_elements(std::move(elements)) {
  // ensure values are sorted by geometry id
  std::sort(m_elements.begin(), m_elements.end(),
            [&](const auto& lhs, const auto& rhs) {
              return getId(lhs) < getId(rhs);
            });
  // check that all elements are unique
  auto duplicate = std::adjacent_find(m_elements.begin(), m_elements.end(),
                                      [&](const auto& lhs, const auto& rhs) {
                                        return getId(lhs) == getId(rhs);
                                      });
  if (duplicate != m_elements.end()) {
    throw std::invalid_argument("Input elements contain duplicates");
  }
  fillLookup(getId);
}

template <typename value_t>
inline typename HierarchicalGeometryContainer<value_t>::Iterator
HierarchicalGeometryContainer<value_t>::find(GeometryID id) const {
  assert((m_elements.size() == m_ids.size()) and
         "Inconsistent container state: #elements != #ids");
  assert((m_elements.size() == m_masks.size()) and
         "Inconsistent container state: #elements != #masks");

  // we can not search for the element directly since the relevant one
  // might be stored at a higher level. ids for higher levels would always
  // be sorted before the requested id. searching for the first element
  // after the requested ensures that we include the full hierarchy.
  const auto it = std::upper_bound(m_ids.begin(), m_ids.end(), id.value());
  auto i = std::distance(m_ids.begin(), it);

  // now go up the hierarchy to find the first matching element.
  // example: the container stores three identifiers
  //
  //     2|x|x (volume-only)
  //     2|3|x (volume and layer)
  //     2|3|4 (volume, layer, and sensitive)
  //
  // where | marks level boundaries. searching for either 2|3|4, 2|3|7, or
  // 2|4|x would first point to 2|3|4 and thus needs to go up the hierarchy.
  while (0 < i) {
    // always starts below item of interest due to upper bound search
    --i;

    // if the input id does not even match at the highest hierarchy level
    // with the current comparison id, then have reached the end of this
    // hierarchy without finding a matching entry.
    // NOTE this prevents adding a catch-all entry (all components equal to
    //      zero), but avoids having an unbounded search window all the way to
    //      the beginning of the ids container.
    if (not equalWithinMask(id.value(), m_ids[i], makeHighestLevelMask())) {
      return end();
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
