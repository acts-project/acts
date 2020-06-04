// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <iostream>
#include <utility>

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelpers.hpp"
#include "Acts/Plugins/TGeo/TGeoSurfaceConverter.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "TGeoArb8.h"
#include "TGeoBBox.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoTrd2.h"
#include "TGeoTube.h"

using Line2D = Eigen::Hyperplane<double, 2>;

Acts::TGeoDetectorElement::TGeoDetectorElement(
    const Identifier& identifier, TGeoNode* tGeoDetElement,
    const TGeoMatrix* mGlobal, const std::string& axes, double scalor,
    bool isDisc, std::shared_ptr<const Acts::ISurfaceMaterial> material,
    std::shared_ptr<const Acts::DigitizationModule> digitizationModule)
    : Acts::IdentifiedDetectorElement(),
      m_detElement(tGeoDetElement),
      m_identifier(identifier),
      m_digitizationModule(std::move(digitizationModule)) {
  using namespace TGeoPrimitivesHelpers;

  // Create temporary local non const surface (to allow setting the material)
  std::shared_ptr<Surface> surface = nullptr;
  // Get the placement and orientation in respect to its mother
  const TGeoMatrix* nodeTransform = (m_detElement->GetMatrix());

  const Double_t* rotation = nullptr;
  const Double_t* translation = nullptr;

  if (mGlobal != nullptr) {
    // the new combined translation
    TGeoHMatrix nTransform =
        TGeoCombiTrans(*mGlobal) * TGeoCombiTrans(*nodeTransform);
    std::string nName = tGeoDetElement->GetName();
    std::string suffix = "_transform";
    nTransform.SetName((nName + suffix).c_str());
    translation = nTransform.GetTranslation();
    rotation = nTransform.GetRotationMatrix();
  } else {
    translation = (nodeTransform->GetTranslation());
    rotation = (nodeTransform->GetRotationMatrix());
  }
  // Simply call the construct method
  construct(rotation, translation, axes, scalor, isDisc, std::move(material));
}

Acts::TGeoDetectorElement::TGeoDetectorElement(
    const Identifier& identifier, const TGeoMatrix& transform,
    TGeoNode* tGeoDetElement, const std::string& axes, double scalor,
    bool isDisc, std::shared_ptr<const Acts::ISurfaceMaterial> material,
    std::shared_ptr<const Acts::DigitizationModule> digitizationModule)
    : Acts::IdentifiedDetectorElement(),
      m_detElement(tGeoDetElement),
      m_identifier(identifier),
      m_digitizationModule(std::move(digitizationModule)) {
  // get the placement and orientation in respect to its mother
  const Double_t* rotation = transform.GetRotationMatrix();
  const Double_t* translation = transform.GetTranslation();
  // Simply call the construct method
  construct(rotation, translation, axes, scalor, isDisc, std::move(material));
}

void Acts::TGeoDetectorElement::construct(
    const Double_t* rotation, const Double_t* translation,
    const std::string& axes, double scalor, bool isDisc,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  using namespace TGeoPrimitivesHelpers;

  // create temporary local non const surface (to allow setting the material)
  std::shared_ptr<Surface> surface = nullptr;

  auto sensor = m_detElement->GetVolume();
  auto tgShape = sensor->GetShape();

  if (not isDisc) {
    auto cylinderComps = TGeoSurfaceConverter::cylinderComponents(
        *tgShape, rotation, translation, axes, scalor);
    auto cylinderBounds = std::get<0>(cylinderComps);
    if (cylinderBounds != nullptr) {
      m_transform = std::get<1>(cylinderComps);
      m_bounds = cylinderBounds;
      m_thickness = std::get<2>(cylinderComps);
      m_surface = Surface::makeShared<CylinderSurface>(cylinderBounds, *this);
    }
  }

  // Check next if you do not have a surface
  if (surface == nullptr) {
    auto discComps = TGeoSurfaceConverter::discComponents(
        *tgShape, rotation, translation, axes, scalor);
    auto discBounds = std::get<0>(discComps);
    if (discBounds != nullptr) {
      m_bounds = discBounds;
      m_transform = std::get<1>(discComps);
      m_thickness = std::get<2>(discComps);
      m_surface = Surface::makeShared<DiscSurface>(discBounds, *this);
    }
  }

  // Check next if you do not have a surface
  if (surface == nullptr) {
    auto planeComps = TGeoSurfaceConverter::planeComponents(
        *tgShape, rotation, translation, axes, scalor);
    auto planeBounds = std::get<0>(planeComps);
    if (planeBounds != nullptr) {
      m_bounds = planeBounds;
      m_transform = std::get<1>(planeComps);
      m_thickness = std::get<2>(planeComps);
      m_surface = Surface::makeShared<PlaneSurface>(planeBounds, *this);
    }
  }

  // set the asscoiated material (non const method)
  if (surface != nullptr) {
    m_surface->assignSurfaceMaterial(std::move(material));
  }
}

Acts::TGeoDetectorElement::~TGeoDetectorElement() = default;
