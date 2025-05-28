// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"

#include <utility>

#include <DD4hep/Alignments.h>
#include <DD4hep/DetElement.h>
#include <DD4hep/Volumes.h>

Acts::DD4hepDetectorElement::DD4hepDetectorElement(
    const dd4hep::DetElement detElement, const std::string& axes, double scalor,
    bool /*isDisc*/, std::shared_ptr<const ISurfaceMaterial> material)
    : TGeoDetectorElement(
          static_cast<TGeoDetectorElement::Identifier>(detElement.volumeID()),
          *(detElement.placement().ptr()),
          detElement.nominal().worldTransformation(), axes, scalor,
          std::move(material)),
      m_detElement(detElement) {}

const Acts::Transform3& Acts::DD4hepDetectorElement::transform(
    const Acts::GeometryContext& gctx) const {
  // Treating non-empty context => contextual alignment present
  if (gctx.hasValue()) {
    const ContextType* dd4hepCtx =
        gctx.maybeGet<DD4hepDetectorElement::ContextType>();
    // Check if a contextual transform is available for this detector element
    auto trfPtr = (dd4hepCtx != nullptr) ? dd4hepCtx->contextualTransform(*this)
                                         : nullptr;
    if (trfPtr != nullptr) {
      return *trfPtr;
    }
  }
  // Empty context, return nominal transform
  return TGeoDetectorElement::transform(gctx);
}
