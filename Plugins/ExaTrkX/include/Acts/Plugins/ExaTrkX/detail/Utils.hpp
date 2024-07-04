// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <ostream>

#include <torch/torch.h>

namespace Acts::detail {

struct TensorDetails {
  const torch::Tensor &tensor;
  TensorDetails(const torch::Tensor &t) : tensor(t) {}
};

inline std::ostream &operator<<(std::ostream &os, const TensorDetails &t) {
  os << t.tensor.dtype() << ", " << t.tensor.sizes();
  if (at::isnan(t.tensor).any().item<bool>()) {
    os << ", contains NaNs";
  } else {
    os << ", no NaNs";
  }
  return os;
}

}  // namespace Acts::detail
