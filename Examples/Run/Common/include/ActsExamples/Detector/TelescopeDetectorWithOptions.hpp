// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
#include "ActsExamples/Utilities/Options.hpp"

namespace ActsExamples {

class TelescopeDetectorWithOptions : public IBaseDetector {
  Telescope::TelescopeDetector m_detector;

 public:
  void addOptions(ActsExamples::Options::Description& desc) const override;

  auto finalize(const boost::program_options::variables_map& vm,
                std::shared_ptr<const Acts::IMaterialDecorator> /*mdecorator*/)
      -> std::pair<TrackingGeometryPtr, ContextDecorators> override;
};
}  // namespace ActsExamples
