// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"

namespace ActsExamples {

/// This is an abstract base class for EDM4hep output converters
/// It's main job is to enforce the presence of  getter for the collection names
/// that converter will write to. This is used by the @c Acts::PodioWriter to
/// pick up the collections to write, as coordinated in python.
class EDM4hepOutputConverter : public IAlgorithm {
 public:
  using IAlgorithm::IAlgorithm;

  /// Get the collection names that this converter writes to.
  ///
  /// @return The collection names
  virtual std::vector<std::string> collections() const = 0;
};

}  // namespace ActsExamples
