// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/JsonDetectorElement.hpp"

#include "ActsPlugins/Json/AlgebraJsonConverter.hpp"
#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"

namespace Acts {

JsonDetectorElement::JsonDetectorElement(const nlohmann::json &jSurface,
                                         double thickness)
    : m_jSurface(jSurface), m_thickness(thickness) {}

const Transform3 &JsonDetectorElement::localToGlobalTransform(
    const GeometryContext & /*gctx*/) const {
  return m_transform;
}

double JsonDetectorElement::thickness() const {
  return m_thickness;
}

std::shared_ptr<Surface> JsonDetectorElement::createSurface() {
  auto surface = Acts::SurfaceJsonConverter::fromJson(m_jSurface);
  m_transform = Transform3JsonConverter::fromJson(m_jSurface["transform"]);

  surface->assignSurfacePlacement(shared_from_this());
  surface->assignThickness(m_thickness);

  return surface;
}

}  // namespace Acts
