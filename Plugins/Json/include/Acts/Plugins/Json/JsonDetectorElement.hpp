// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

class JsonDetectorElement : public DetectorElementBase {
  std::shared_ptr<Surface> m_surface;
  Transform3 m_transform{};
  double m_thickness{};

 public:
  JsonDetectorElement(const nlohmann::json &jSurface, double thickness);

  Surface &surface() override;
  const Surface &surface() const override;

  double thickness() const override;

  const Transform3 &transform(const GeometryContext &gctx) const override;
};

}  // namespace Acts
