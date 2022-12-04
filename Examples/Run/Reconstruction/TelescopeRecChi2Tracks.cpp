// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/TelescopeDetectorWithOptions.hpp"

#include "RecChi2Tracks.hpp"

int main(int argc, char* argv[]) {
  return runRecChi2Tracks(
      argc, argv,
      std::make_shared<ActsExamples::TelescopeDetectorWithOptions>());
}
