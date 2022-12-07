// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "HoughExample.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Detector/GenericDetectorWithOptions.hpp"

int main(int argc, char* argv[]) {
  return runHoughExample(
      argc, argv, std::make_shared<ActsExamples::GenericDetectorWithOptions>());
}
