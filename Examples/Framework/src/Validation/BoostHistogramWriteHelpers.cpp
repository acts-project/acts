// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/BoostHistogramWriteHelpers.hpp"

#include "ActsExamples/Validation/BoostHistogramToRootConverter.hpp"
#include "ActsExamples/Validation/BoostHistogramWrappers.hpp"

#include <TEfficiency.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

namespace ActsExamples::BoostHistogramWriteHelpers {

void write(const BoostHistogram1D& hist) {
  TH1F* rootHist = BoostHistogramToRoot::toRoot(hist);
  rootHist->Write();
  delete rootHist;
}

void write(const BoostHistogram2D& hist) {
  TH2F* rootHist = BoostHistogramToRoot::toRoot(hist);
  rootHist->Write();
  delete rootHist;
}

void write(const BoostProfileHistogram& hist) {
  TProfile* rootProfile = BoostHistogramToRoot::toRoot(hist);
  rootProfile->Write();
  delete rootProfile;
}

void write(const BoostEfficiency1D& hist) {
  TEfficiency* rootEff = BoostHistogramToRoot::toRoot(hist);
  rootEff->Write();
  delete rootEff;
}

void write(const BoostEfficiency2D& hist) {
  TEfficiency* rootEff = BoostHistogramToRoot::toRoot(hist);
  rootEff->Write();
  delete rootEff;
}

}  // namespace ActsExamples::BoostHistogramWriteHelpers
