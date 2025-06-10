// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/JsonDetectorElement.hpp"

#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"

namespace Acts {

JsonDetectorElement::JsonDetectorElement(const nlohmann::json &jSurface,
                                         double thickness)
    : m_thickness(thickness) {
  m_surface = Acts::SurfaceJsonConverter::fromJson(jSurface);
  m_transform = Transform3JsonConverter::fromJson(jSurface["transform"]);
  m_surface->assignDetectorElement(*this);
}

const Surface &JsonDetectorElement::surface() const {
  return *m_surface;
}

Surface &JsonDetectorElement::surface() {
  return *m_surface;
}

const Transform3 &JsonDetectorElement::transform(
    const GeometryContext & /*gctx*/) const {
  return m_transform;
}

double JsonDetectorElement::thickness() const {
  return m_thickness;
}

}  // namespace Acts
