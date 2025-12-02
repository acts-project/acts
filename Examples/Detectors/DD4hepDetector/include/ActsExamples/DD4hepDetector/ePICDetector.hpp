// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"

namespace Acts {
class GeometryContext;
}

namespace ActsPlugins {
class DD4hepDetectorElement;
}

namespace ActsExamples {

class ePICDetector final : public DD4hepDetectorBase {
 public:
  struct Config : public DD4hepDetectorBase::Config {
    using ElementFactory =
        std::function<std::shared_ptr<ActsPlugins::DD4hepDetectorElement>(
            const dd4hep::DetElement& element, const std::string& axes,
            double scale)>;

    ElementFactory detectorElementFactory = defaultDetectorElementFactory;
  };

  static std::shared_ptr<ActsPlugins::DD4hepDetectorElement>
  defaultDetectorElementFactory(const dd4hep::DetElement& element,
                                const std::string& axes, double scale);

  explicit ePICDetector(const Config& cfg,
                        const Acts::GeometryContext& gctx);

  const Config& config() const override;

 private:
  void construct(const Acts::GeometryContext& gctx);

  Config m_cfg;
};

}  // namespace ActsExamples
