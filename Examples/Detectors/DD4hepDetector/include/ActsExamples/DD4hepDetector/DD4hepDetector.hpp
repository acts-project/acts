// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <memory>
#include <vector>

namespace ActsExamples {
namespace DD4hep {

struct DD4hepDetector : public IBaseDetector {
  std::shared_ptr<DD4hepGeometryService> geometryService;

  void addOptions(
      boost::program_options::options_description& opt) const override;

  std::pair<IBaseDetector::TrackingGeometryPtr, ContextDecorators> finalize(
      const boost::program_options::variables_map& vm,
      std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) override;

  std::pair<IBaseDetector::TrackingGeometryPtr, ContextDecorators> finalize(
      DD4hepGeometryService::Config cfg,
      std::shared_ptr<const Acts::IMaterialDecorator> mdecorator);
};

}  // namespace DD4hep
}  // namespace ActsExamples
