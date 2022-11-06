// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/TrackFittingChi2Options.hpp"

#include "ActsExamples/Utilities/Options.hpp"

using namespace boost::program_options;

void ActsExamples::Options::addFittingChi2Options(
    boost::program_options::options_description& opt) {
  opt.add_options()("chi2-updates", value<unsigned int>()->default_value(1),
                    "number of update steps for the Chi2Fitter");
}
