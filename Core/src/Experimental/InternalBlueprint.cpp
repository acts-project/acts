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
    const SurfaceLinksGenerator& surfaceLinksGenerator,
    const std::vector<std::pair<BinningValue, Envelope>>& binValueEnvelopes,
    const std::string& name)
    : m_surfaces(surfaces),
      m_surfaceLinksGenerator(surfaceLinksGenerator),
      m_name(name) {
  // Create the envelope & store the bin values
  ExtentEnvelope eEnvelope = zeroEnvelopes;
  std::vector<BinningValue> binValues;
  for (const auto& bve : binValueEnvelopes) {
    binValues.push_back(bve.first);
    eEnvelope[bve.first] = bve.second;
  }
  m_extent = GeometricExtent(eEnvelope);
  m_surfaces.reserve(surfaces.size());

  for (auto& sf : m_surfaces) {
    auto sfPolyhedron = sf->polyhedronRepresentation(gctx, 1);
    for (const auto& v : sfPolyhedron.vertices) {
      m_extent.extend(v, binValues);
    }
  }
}

Acts::InternalBlueprint::InternalBlueprint(
    const std::vector<std::pair<std::shared_ptr<Surface>, GeometricExtent>>&
        surfacesExtent,
    const SurfaceLinksGenerator& surfaceLinksGenerator, const std::string& name)
    : m_surfaceLinksGenerator(surfaceLinksGenerator), m_name(name) {
  // Store surfaces and extent
  m_surfaces.reserve(surfacesExtent.size());
  for (auto& sfe : surfacesExtent) {
    m_surfaces.push_back(sfe.first);
    m_extent.extend(sfe.second, s_binningValues, false);
  }
}

std::ostream& Acts::InternalBlueprint::toStream(std::ostream& sl) const {
  sl << "InternalBlueprint with dimensions (min/max)" << std::endl;
  m_extent.toStream(sl);
  return sl;
}

// Overload of << operator for std::ostream for debug output
std::ostream& Acts::operator<<(std::ostream& sl,
                               const Acts::InternalBlueprint& lbp) {
  return lbp.toStream(sl);
}
