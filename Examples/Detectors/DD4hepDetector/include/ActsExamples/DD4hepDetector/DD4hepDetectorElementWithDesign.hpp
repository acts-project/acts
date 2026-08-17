// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/ISensorDesign.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"

#include <memory>

namespace ActsExamples {

/// DD4hepDetectorElement extended with a sensor design pointer.
/// Created by the custom factory below instead of the default one.
class DD4hepDetectorElementWithDesign
    : public ActsPlugins::DD4hepDetectorElement {
 public:
  // Inherit constructors
  using ActsPlugins::DD4hepDetectorElement::DD4hepDetectorElement;

  /// Access the sensor design attached to this element, can be nullptr
  const Acts::ISensorDesign* sensorDesign() const { return m_design.get(); }

  /// Attach a design — const because m_design is mutable (setup-time
  /// annotation)
  void assignDesign(std::shared_ptr<const Acts::ISensorDesign> design) const {
    m_design = std::move(design);
  }

  /// Factory function — same signature as DD4hepLayerBuilder::ElementFactory
  static std::shared_ptr<ActsPlugins::DD4hepDetectorElement> factory(
      const dd4hep::DetElement& detElement, TGeoAxes axes, double scale,
      std::shared_ptr<const Acts::ISurfaceMaterial> material) {
    return std::make_shared<DD4hepDetectorElementWithDesign>(
        detElement, axes, scale, std::move(material));
  }

 private:
  mutable std::shared_ptr<const Acts::ISensorDesign> m_design{nullptr};
};

/// Drop-in replacement for DD4hepLayerBuilder::defaultDetectorElementFactory.
/// Creates DD4hepDetectorElementWithDesign so surfaces carry a design pointer.
inline std::shared_ptr<ActsPlugins::DD4hepDetectorElement>
defaultDetectorElementFactoryWithDesign(
    const dd4hep::DetElement& detElement, ActsPlugins::TGeoAxes axes,
    double scale, std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  return std::make_shared<DD4hepDetectorElementWithDesign>(
      detElement, axes, scale, std::move(material));
}

}  // namespace ActsExamples
