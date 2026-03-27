// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DetectorCommons/Aligned.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetectorElement.hpp"

namespace ActsExamples {
/// Define the aligned telescope detector element and factory type
using AlignedTelescopeDetectorElement = Aligned<TelescopeDetectorElement>;

class AlignedTelescopeDetector : public TelescopeDetector {
 public:
  using Config = TelescopeDetector::Config;
  explicit AlignedTelescopeDetector(const Config& cfg);
};

}  // namespace ActsExamples
