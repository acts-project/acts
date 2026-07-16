// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/IDesign.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"

#include <string_view>

namespace ActsExamples {

class DigitizationDesign final : public Acts::IDesign {
 public:
  explicit DigitizationDesign(const DigitizationAlgorithm::Digitizer* digitizer)
      : m_digitizer(digitizer) {}

  const DigitizationAlgorithm::Digitizer& digitizer() const {
    return *m_digitizer;
  }

  std::string_view name() const override { return "DigitizationDesign"; }

 private:
  const DigitizationAlgorithm::Digitizer* m_digitizer;  // non-owning
};

}  // namespace ActsExamples
