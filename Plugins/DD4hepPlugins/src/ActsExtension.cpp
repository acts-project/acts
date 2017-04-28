// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/ActsExtension.hpp"
#include <boost/algorithm/string.hpp>
#include "ACTS/Digitization/CartesianSegmentation.hpp"
#include "ACTS/Digitization/DigitizationModule.hpp"
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/TrapezoidBounds.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "DD4hep/CartesianGridXY.h"

Acts::ActsExtension::ActsExtension(const Config& cfg) : Acts::IActsExtension()
{
  setConfiguration(cfg);
}

Acts::ActsExtension::ActsExtension(DD4hep::Geometry::Segmentation segmentation,
                                   DD4hep::Geometry::Volume       volume,
                                   std::string                    axes)
  : Acts::IActsExtension()
{
  setSegmentation(segmentation, volume, axes);
}

Acts::ActsExtension::ActsExtension(const ActsExtension& det,
                                   const DD4hep::Geometry::DetElement&)
  : Acts::IActsExtension(), m_cfg(det.m_cfg)
{
}

void
Acts::ActsExtension::setConfiguration(const Acts::ActsExtension::Config& config)
{
  // @todo check consistency
  // copy the configuration
  m_cfg = config;
}

void
Acts::ActsExtension::setSegmentation(
    DD4hep::Geometry::Segmentation segmentation,
    DD4hep::Geometry::Volume       volume,
    std::string                    axes)
{
  m_cfg.axes = axes;
  // find out the shape and create the acts bounds
  std::shared_ptr<const PlanarBounds> bounds    = nullptr;
  double                              thickness = 0.;
  double                              scalor    = units::_cm;
  double                              l0;
  double                              l1;

  TGeoBBox* box       = dynamic_cast<TGeoBBox*>(volume->GetShape());
  TGeoTrd2* trapezoid = dynamic_cast<TGeoTrd2*>(volume->GetShape());

  //
  if (boost::iequals(m_cfg.axes, "XYZ")) {
    if (trapezoid) {
      // bounds with x/y
      double x1 = scalor * trapezoid->GetDx1();
      double x2 = scalor * trapezoid->GetDx2();
      double y  = scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2());
      l0        = x2;
      l1        = y;
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(x1, x2, y);
      // thickness
      thickness = scalor * trapezoid->GetDz();
      // assign them
      bounds = trapezoidBounds;
    } else {
      l0 = scalor * box->GetDX();
      l1 = scalor * box->GetDY();
      // bounds with x/y
      auto rectangleBounds = std::make_shared<const RectangleBounds>(l0, l1);
      // thickness
      thickness = scalor * box->GetDZ();
      // assign them
      bounds = rectangleBounds;
    }
  } else if (boost::iequals(m_cfg.axes, "XZY")) {
    if (trapezoid) {
      double x1 = scalor * trapezoid->GetDx1();
      double x2 = scalor * trapezoid->GetDx2();
      double y  = scalor * trapezoid->GetDz();
      l0        = x2;
      l1        = y;

      // bounds with x/z
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(x1, x2, y);
      // thickness
      thickness = scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2());
      // assign them
      bounds = trapezoidBounds;
    } else {
      l0 = scalor * box->GetDX();
      l1 = scalor * box->GetDZ();
      // bounds with x/z
      auto rectangleBounds = std::make_shared<const RectangleBounds>(l0, l1);
      // thickness
      thickness = scalor * box->GetDY();
      // assign them
      bounds = rectangleBounds;
    }

  } else if (boost::iequals(m_cfg.axes, "YZX")) {
    if (trapezoid) {
      double x1 = scalor * trapezoid->GetDy1();
      double x2 = scalor * trapezoid->GetDy2();
      double y  = scalor * trapezoid->GetDz();
      l0        = x2;
      l1        = y;
      // bounds with y/z
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(x1, x2, y);
      // thickness
      thickness = scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2());
      // assign them
      bounds = trapezoidBounds;
    } else {
      l0 = scalor * box->GetDY();
      l1 = scalor * box->GetDZ();
      // bounds with y/z
      auto rectangleBounds = std::make_shared<const RectangleBounds>(l0, l1);
      // thickness
      thickness = scalor * box->GetDX();
      // assign them
      bounds = rectangleBounds;
    }
  } else if (boost::iequals(m_cfg.axes, "YXZ")) {
    if (trapezoid) {
      double x1 = scalor * trapezoid->GetDy1();
      double x2 = scalor * trapezoid->GetDy2();
      double y  = scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2());
      l0        = x2;
      l1        = y;
      // bounds with y/x
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(x1, x2, y);
      // thickness
      thickness = scalor * trapezoid->GetDz();
      // assign them
      bounds = trapezoidBounds;
    } else {
      l1 = scalor * box->GetDY();
      l1 = scalor * box->GetDX();
      // bounds with y/x
      auto rectangleBounds = std::make_shared<const RectangleBounds>(l0, l1);
      // thickness
      thickness = scalor * box->GetDZ();
      // assign them
      bounds = rectangleBounds;
    }
  } else if (boost::iequals(m_cfg.axes, "ZYX")) {
    if (trapezoid) {
      double x1 = scalor * trapezoid->GetDz();
      double x2 = scalor * trapezoid->GetDz();
      double y  = scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2());
      l0        = x2;
      l1        = y;
      // bounds with z/y
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(x1, x2, y);
      // thickness
      thickness = scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2());
      // assign them
      bounds = trapezoidBounds;
    } else {
      l0 = scalor * box->GetDZ();
      l1 = scalor * box->GetDY();
      // bounds with z/y
      auto rectangleBounds = std::make_shared<const RectangleBounds>(l0, l1);
      // thickness
      thickness = scalor * box->GetDX();
      // assign them
      bounds = rectangleBounds;
    }
  } else {
    if (trapezoid) {
      double x1 = scalor * trapezoid->GetDz();
      double x2 = scalor * trapezoid->GetDz();
      double y  = scalor * 0.5 * (trapezoid->GetDx1() + trapezoid->GetDx2());

      // bounds with z/x
      auto trapezoidBounds = std::make_shared<const TrapezoidBounds>(x1, x2, y);
      // thickness
      thickness = scalor * 0.5 * (trapezoid->GetDy1() + trapezoid->GetDy2());
      // assign them
      bounds = trapezoidBounds;
    } else {
      l0 = scalor * box->GetDZ();
      l1 = scalor * box->GetDX();
      // bounds with z/x
      auto rectangleBounds = std::make_shared<const RectangleBounds>(l0, l1);
      // thickness
      thickness = scalor * box->GetDY();
      // assign them
      bounds = rectangleBounds;
    }
  }

  DD4hep::Geometry::CartesianGridXY cartesianGrid = segmentation;
  if (cartesianGrid.isValid()) {
    // the Acts segmentation of the DigitizationModule

    size_t bins0
        = (cartesianGrid.gridSizeX() != 0) ? l0 / cartesianGrid.gridSizeX() : 0;
    size_t bins1
        = (cartesianGrid.gridSizeY() != 0) ? l1 / cartesianGrid.gridSizeY() : 0;

    std::shared_ptr<const CartesianSegmentation> actsSegmentation
        = std::make_shared<const CartesianSegmentation>(bounds, bins0, bins1);
    // finally create the digitization module
    // @todo set lorentz angle
    m_digiModule = std::make_shared<const DigitizationModule>(
        actsSegmentation, thickness, 1, 0);
  }
}
