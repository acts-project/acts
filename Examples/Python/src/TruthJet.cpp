// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Jets/TruthJetAlgorithm.hpp"

#include <cstddef>
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addTruthJet(Context& ctx) {
  auto mex = ctx.get("examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TruthJetAlgorithm, mex,
                                "TruthJetAlgorithm", inputTruthParticles,
                                outputJets, jetPtMin);
}  // addTruthJet
}  // namespace Acts::Python
