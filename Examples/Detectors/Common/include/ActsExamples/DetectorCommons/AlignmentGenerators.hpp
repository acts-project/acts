// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <memory>
#include <unordered_map>
#include <vector>

namespace Acts {
class TrackingGeometry;
class ITransformStore;
}  // namespace Acts

namespace ActsExamples {

struct NoRandom {
  /// @brief A no random number generator
  /// @return 1.0
  double operator()() const { return 1.0; }
};

/// @brief A simple alignment generator for the geometry to apply a global shift
struct GlobalShift {
  /// @brief The configuration struct
  GlobalShift() = default;

  // Constructor with arguments
  /// @param gctx The geometry context - for the nominal transforms
  /// @param trackingGeometry The tracking geometry
  /// @param selection The selection of elements to be shifted via GeometryIdentifier
  /// @param ishift The shift vector
  GlobalShift(const Acts::GeometryContext& gctx,
              const Acts::TrackingGeometry& trackingGeometry,
              const std::vector<Acts::GeometryIdentifier>& selection,
              const Acts::Vector3& ishift);

  // The shift to be applied
  Acts::Vector3 shift = Acts::Vector3::Zero();

  /// Internal storage - the nominal transforms of the selected elements
  std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>
      nominalTransforms;

  /// The call operator to generate the transform store
  /// @param rng The random number generator
  std::shared_ptr<Acts::ITransformStore> operator()(
      std::function<double()>& rng);
};

/// @brief A simple alignment generator for the geometry to apply a radial expansion
struct PerpendicularScale {
  /// @brief The configuration struct
  PerpendicularScale() = default;

  // Constructor with arguments
  /// @param gctx The geometry context - for the nominal transforms
  /// @param trackingGeometry The tracking geometry
  /// @param selection The selection of elements to be shifted via GeometryIdentifier
  /// @param iexpansion The radial expansion factor
  PerpendicularScale(const Acts::GeometryContext& gctx,
                     const Acts::TrackingGeometry& trackingGeometry,
                     const std::vector<Acts::GeometryIdentifier>& selection,
                     double iexpansion);

  // The expansion to be applied
  double expansion = 1.0;

  /// Internal storage - the nominal transforms of the selected elements
  std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>
      nominalTransforms;

  /// The call operator to generate the transform store
  std::shared_ptr<Acts::ITransformStore> operator()(
      std::function<double()>& rng);
};

}  // namespace ActsExamples
