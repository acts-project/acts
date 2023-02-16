// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/Smearers.hpp"

namespace ActsExamples::Digitization {

bool compareSmearers(
    const ActsFatras::SingleParameterSmearFunction<RandomEngine> &a,
    const ActsFatras::SingleParameterSmearFunction<RandomEngine> &b) {
  // First check if both have a target
  if (static_cast<bool>(a) != static_cast<bool>(b)) {
    return false;
  }

  // If both have no target, also equal
  if (not a && not b) {
    return true;
  }

  // I think this is the only way to really compare the smearers. Not
  // sure if this might be so slow to just avoid this and always?
  if (a.target<Exact>() == b.target<Exact>()) {
    return true;
  }
  if (a.target<Gauss>() == b.target<Gauss>()) {
    return true;
  }
  if (a.target<GaussTrunc>() == b.target<GaussTrunc>()) {
    return true;
  }
  if (a.target<GaussClipped>() == b.target<GaussClipped>()) {
    return true;
  }
  if (a.target<Uniform>() == b.target<Uniform>()) {
    return true;
  }
  if (a.target<Digital>() == b.target<Digital>()) {
    return true;
  }

  return false;
}

}  // namespace ActsExamples::Digitization
