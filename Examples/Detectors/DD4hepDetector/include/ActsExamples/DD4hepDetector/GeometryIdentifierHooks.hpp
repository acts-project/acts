// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <map>
#include <vector>

/// Geometry identifiers hooks to be used with DD4Hep detectors
/// to add some extra identifier to sensitives surfaces.
namespace det {
namespace GeometryIdentifierHooks {

struct RadiusGeometryIdentifierDecorator final : Acts::GeometryIdentifierHook {
  std::map<int, std::vector<double>> volumeToRadialCuts;

  /// Use the extra identifier in the ODD endcap to separate the two row of
  /// modules in R
  /// @param identifier geometry identifier
  /// @param surface coresponding surface
  virtual Acts::GeometryIdentifier decorateIdentifier(
      Acts::GeometryIdentifier id, const Acts::Surface&) const override;
};

}  // namespace GeometryIdentifierHooks
}  // namespace det
