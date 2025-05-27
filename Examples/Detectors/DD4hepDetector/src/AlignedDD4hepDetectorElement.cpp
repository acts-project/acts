// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/AlignedDD4hepDetectorElement.hpp"

std::shared_ptr<ActsExamples::AlignedDD4hepDetectorElement>
ActsExamples::alignedDD4hepDetectorElementFactory(
    const dd4hep::DetElement detElement, const std::string& axes, double scalor,
    bool isDisc, std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  return std::make_shared<AlignedDD4hepDetectorElement>(
      detElement, axes, scalor, isDisc, std::move(material));
};
