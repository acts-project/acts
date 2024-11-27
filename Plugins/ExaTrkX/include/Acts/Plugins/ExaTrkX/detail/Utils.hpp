// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

template <typename It>
struct RangePrinter {
  It begin;
  It end;

  RangePrinter(It a, It b) : begin(a), end(b) {}
};

template <class It>
RangePrinter(It b, It e) -> RangePrinter<It>;

template <typename It>
inline std::ostream &operator<<(std::ostream &os, const RangePrinter<It> &r) {
  for (auto it = r.begin; it != r.end; ++it) {
    os << *it << " ";
  }
  return os;
}

}  // namespace Acts::detail
