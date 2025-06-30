// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DetectorCommons/Aligned.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorElement.hpp"

namespace ActsExamples {
/// Define the aligned DD4hep detector element and factory type
using AlignedGenericDetectorElement = Aligned<GenericDetectorElement>;

class AlignedGenericDetector : public GenericDetector {
 public:
  using Config = GenericDetector::Config;
  explicit AlignedGenericDetector(const Config& cfg);
};

}  // namespace ActsExamples
