// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cassert>

namespace ActsExamples {

const Acts::Surface* IndexSourceLinkSurfaceAccessor::operator()(
    const Acts::SourceLink& sourceLink) const {
  const auto& indexSourceLink = sourceLink.get<IndexSourceLink>();
  return geometry.findSurface(indexSourceLink.geometryId());
}

std::pair<IndexSourceLinkAccessor::Iterator, IndexSourceLinkAccessor::Iterator>
IndexSourceLinkAccessor::range(const Acts::Surface& surface) const {
  assert(container != nullptr);
  auto [begin, end] = container->equal_range(surface.geometryId());
  return {Iterator{begin}, Iterator{end}};
}

}  // namespace ActsExamples
