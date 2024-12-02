// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/MultiIndex.hpp"

namespace ActsExamples {

class ITkIdentifier {
  Acts::MultiIndex<std::size_t, 1, 2, 20, 1, 19, 20, 1> m_identifier{};

 public:
  ITkIdentifier(int hardware, int barrelEndcap, int layerWheel, int etaModule,
                int phiModule, int side);

  /// Access the hardware specifier (pixel=0, strip=1)
  int hardware() const;

  /// Access the barrel-endcap specifier (-2,0,2)
  int barrelEndcap() const;

  /// Access the layer specifier
  int layerWheel() const;

  /// Access the phi module specifier
  int phiModule() const;

  /// Access the eta module specifier
  int etaModule() const;

  /// Access the side (for double sided strip modules)
  int side() const;

  /// A unique identifier that represents the combination of specifiers
  std::size_t value() const;
};

std::ostream& operator<<(std::ostream& os, const ITkIdentifier& id);

/// Wrapper detector element that stores additional information for the ITk
/// Probably not optimal performancewise due to pointer chasing...
class ITkDetectorElement : public Acts::DetectorElementBase {
 public:
  ITkDetectorElement(std::shared_ptr<DetectorElementBase> element, int hardware,
                     int barrelEndcap, int layerWheel, int etaModule,
                     int phiModule, int side)
      : m_element(element),
        m_identifier(hardware, barrelEndcap, layerWheel, etaModule, phiModule,
                     side) {
    m_element->surface().assignDetectorElement(*this);
  }

  ITkIdentifier identifier() const { return m_identifier; }

  double thickness() const override { return m_element->thickness(); }

  Acts::Surface& surface() override { return m_element->surface(); }

  const Acts::Surface& surface() const override { return m_element->surface(); }

  const Acts::Transform3& transform(
      const Acts::GeometryContext& gctx) const override {
    return m_element->transform(gctx);
  }

 private:
  std::shared_ptr<DetectorElementBase> m_element;
  ITkIdentifier m_identifier;
};

}  // namespace ActsExamples
