// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsPlugins/Arrow/Export.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

/// Abstract base class for Arrow output converters.
///
/// Its only job is to advertise, via @c collections(), the whiteboard keys
/// this converter will emit tables under. The @c ParquetWriter reads those
/// keys from its Python-side configuration; this mirrors how
/// @c PodioOutputConverter declares podio collection names.
class ACTS_ARROW_EXPORT ArrowOutputConverter : public IAlgorithm {
 public:
  using IAlgorithm::IAlgorithm;

  /// Get the whiteboard keys that this converter will write tables to.
  ///
  /// @return The collection names.
  virtual std::vector<std::string> collections() const = 0;
};

}  // namespace ActsExamples
