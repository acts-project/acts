// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <vector>

namespace ActsExamples {

/// A simultated vertex e.g. from a physics process.
struct SimVertex {
  using Scalar = Acts::ActsScalar;
  using Vector4 = Acts::ActsVector<4>;

  /// The vertex four-position.
  Vector4 position4 = Vector4::Zero();
  /// The vertex process type.
  ActsFatras::ProcessType process = ActsFatras::ProcessType::eUndefined;
  /// The incoming particles into the vertex
  std::vector<ActsFatras::Particle> incoming = {};
  /// The outgoing particles from the vertex
  std::vector<ActsFatras::Particle> outgoing = {};

  /// Construct the vertex from a position and optional process type.
  ///
  /// @param position4_ the vertex four-position
  /// @param process_ the process type that generated this vertex
  ///
  /// Associated particles are left empty by default and must be filled by the
  /// user after construction.
  SimVertex(const Vector4& position4_, ActsFatras::ProcessType process_ =
                                           ActsFatras::ProcessType::eUndefined)
      : position4(position4_), process(process_) {}
  // explicitly default rule-of-five.
  SimVertex() = default;
  SimVertex(const SimVertex&) = default;
  SimVertex(SimVertex&&) = default;
  SimVertex& operator=(const SimVertex&) = default;
  SimVertex& operator=(SimVertex&&) = default;

  /// The vertex three-position.
  auto position() const { return position4.head<3>(); }
  /// The vertex time.
  Scalar time() const { return position4[3]; }
};

}  // namespace ActsExamples
