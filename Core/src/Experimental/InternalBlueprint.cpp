// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/InternalBlueprint.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <iostream>

Acts::InternalBlueprint::InternalBlueprint(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<Surface>>& surfaces,
    const SurfaceLinks& surfaceLinks,
    const Extent& restriction,
    const std::vector<std::pair<ActsScalar, ActsScalar>>& envelope,
    const std::string& name)
    : m_restriction(restriction), m_envelope(envelope), m_surfaceLinks(surfaceLinks), m_name(name)
{
  m_surfaces.reserve(surfaces.size());
  measure(gctx, surfaces);
}

Acts::ActsScalar Acts::InternalBlueprint::min(BinningValue bval, bool addenv) const {
  if (addenv) {
    return m_extent.min(bval) - m_envelope[bval].first;
  }
  return m_extent.min(bval);
}

Acts::ActsScalar Acts::InternalBlueprint::max(BinningValue bval, bool addenv) const {
  if (addenv) {
    return m_extent.max(bval) + m_envelope[bval].second;
  }
  return m_extent.max(bval);
}

Acts::ActsScalar Acts::InternalBlueprint::medium(BinningValue bval, bool addenv) const {
  return 0.5 * (min(bval, addenv) + max(bval, addenv));
}

Acts::ActsScalar Acts::InternalBlueprint::range(BinningValue bval, bool addenv) const {
  return std::abs(max(bval, addenv) - min(bval, addenv));
}

std::ostream& Acts::InternalBlueprint::toStream(std::ostream& sl) const {
  sl << "InternalBlueprint with dimensions (min/max)" << std::endl;
  m_extent.toStream(sl);
  return sl;
}

void Acts::InternalBlueprint::measure(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<Surface>>& surfaces) {

  for (auto& sf : surfaces) {
    auto sfPolyhedron = sf->polyhedronRepresentation(gctx, 1);    
    auto sfExtent= sfPolyhedron.extent();
    // Check if the surface is contained within restriction
    if (m_restriction.contains(sfExtent)){
        m_extent.extend(sfExtent);
        m_surfaces.push_back(sf);
    }
  }
}

void Acts::InternalBlueprint::add(const GeometryContext& gctx,
                               std::shared_ptr<Surface> surface) {
  measure(gctx, { surface });
}

// Overload of << operator for std::ostream for debug output
std::ostream& Acts::operator<<(std::ostream& sl, const Acts::InternalBlueprint& lbp) {
  return lbp.toStream(sl);
}
