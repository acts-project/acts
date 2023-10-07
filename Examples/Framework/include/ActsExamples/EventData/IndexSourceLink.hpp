// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"

#include <cassert>

namespace ActsExamples {

/// A source link that stores just an index.
///
/// This is intentionally kept as barebones as possible. The source link
/// is just a reference and will be copied, moved around, etc. often.
/// Keeping it small and separate from the actual, potentially large,
/// measurement data should result in better overall performance.
/// Using an index instead of e.g. a pointer, means source link and
/// measurement are decoupled and the measurement represenation can be
/// easily changed without having to also change the source link.
class IndexSourceLink final {
 public:
  /// Construct from geometry identifier and index.
  constexpr IndexSourceLink(Acts::GeometryIdentifier gid, Index idx)
      : m_geometryId(gid), m_index(idx) {}

  // Construct an invalid source link. Must be default constructible to
  /// satisfy SourceLinkConcept.
  IndexSourceLink() = default;
  IndexSourceLink(const IndexSourceLink&) = default;
  IndexSourceLink(IndexSourceLink&&) = default;
  IndexSourceLink& operator=(const IndexSourceLink&) = default;
  IndexSourceLink& operator=(IndexSourceLink&&) = default;

  /// Access the index.
  constexpr Index index() const { return m_index; }

  Acts::GeometryIdentifier geometryId() const { return m_geometryId; }

 private:
  Acts::GeometryIdentifier m_geometryId;
  Index m_index = 0;

  friend bool operator==(const IndexSourceLink& lhs,
                         const IndexSourceLink& rhs) {
    return (lhs.geometryId() == rhs.geometryId()) and
           (lhs.m_index == rhs.m_index);
  }
  friend bool operator!=(const IndexSourceLink& lhs,
                         const IndexSourceLink& rhs) {
    return not(lhs == rhs);
  }
};

/// Container of index source links.
///
/// Since the source links provide a `.geometryId()` accessor, they can be
/// stored in an ordered geometry container.
using IndexSourceLinkContainer = GeometryIdMultiset<IndexSourceLink>;
/// Accessor for the above source link container
///
/// It wraps up a few lookup methods to be used in the Combinatorial Kalman
/// Filter
struct IndexSourceLinkAccessor : GeometryIdMultisetAccessor<IndexSourceLink> {
  using BaseIterator = GeometryIdMultisetAccessor<IndexSourceLink>::Iterator;

  using Iterator = Acts::SourceLinkAdapterIterator<BaseIterator>;

  // get the range of elements with requested geoId
  std::pair<Iterator, Iterator> range(const Acts::Surface& surface) const {
    assert(container != nullptr);
    auto [begin, end] = container->equal_range(surface.geometryId());
    return {Iterator{begin}, Iterator{end}};
  }
};
}  // namespace ActsExamples
