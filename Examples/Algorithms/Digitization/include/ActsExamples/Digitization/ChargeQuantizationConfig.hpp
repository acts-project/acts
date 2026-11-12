// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace ActsExamples {

/// Linear readout of accepted, aggregated channel activations.
/// Zero is reserved for absent cells; retained cells have codes 1..2^bits-1.
struct ChargeQuantizationConfig {
  /// At most 24 bits so integer codes remain exact in float32 cell writers.
  unsigned int bits = 4;
  double offset = 0.;
  double step = 1.;

  bool operator==(const ChargeQuantizationConfig &) const = default;

  /// Reject unsupported bit depths and invalid activation binning.
  /// Called once per configuration entry by DigitizationAlgorithm at startup.
  void validate() const;

  /// Encode once after same-channel summation. Boundary ties go upwards.
  double quantize(double activation) const;
};

}  // namespace ActsExamples
