// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "DD4hepTestsHelper.hpp"

#include "Acts/Plugins/DD4hep/DD4hepConversionHelpers.hpp"
#include "Acts/Surfaces/Surface.hpp"

dd4hep::Transform3D DD4hepTestsHelper::createTransform(
    const xml_comp_t& x_det_comp) {
  // Build the transform - center def
  double cx = Acts::getAttrValueOr<double>(x_det_comp, "cx", 0.);
  double cy = Acts::getAttrValueOr<double>(x_det_comp, "cy", 0.);
  double cz = Acts::getAttrValueOr<double>(x_det_comp, "cz", 0.);

  double xx = Acts::getAttrValueOr<double>(x_det_comp, "xx", 1.);
  double xy = Acts::getAttrValueOr<double>(x_det_comp, "xy", 0.);
  double xz = Acts::getAttrValueOr<double>(x_det_comp, "xz", 0.);

  double yx = Acts::getAttrValueOr<double>(x_det_comp, "yx", 0.);
  double yy = Acts::getAttrValueOr<double>(x_det_comp, "yy", 1.);
  double yz = Acts::getAttrValueOr<double>(x_det_comp, "yz", 0.);

  Position xAxis(xx, xy, xz);
  Position yAxis(yx, yy, yz);
  Position zAxis = xAxis.Cross(yAxis);
  double zx = zAxis.X();
  double zy = zAxis.Y();
  double zz = zAxis.Z();

  // Create the transform
  return Transform3D(xx, yx, zx, cx, xy, yy, zy, cy, xz, yz, zz, cz);
}

std::string DD4hepTestsHelper::transformToXML(const Acts::Transform3& tf,
                                              const std::array<int, 2u>& axes) {
  auto tr = tf.translation();
  auto rot = tf.rotation();

  std::stringstream sxml;
  sxml << "cx=\"" << tr[0u] << "*mm\" ";
  sxml << "cy=\"" << tr[1u] << "*mm\" ";
  sxml << "cz=\"" << tr[2u] << "*mm\" ";

  sxml << "xx=\"" << rot.col(axes[0u])[0u] << "\" ";
  sxml << "xy=\"" << rot.col(axes[0u])[1u] << "\" ";
  sxml << "xz=\"" << rot.col(axes[0u])[2u] << "\" ";
  sxml << "yx=\"" << rot.col(axes[1u])[0u] << "\" ";
  sxml << "yy=\"" << rot.col(axes[1u])[1u] << "\" ";
  sxml << "yz=\"" << rot.col(axes[1u])[2u] << "\" ";

  return sxml.str();
}

std::string DD4hepTestsHelper::surfaceToXML(const Acts::GeometryContext& gctx,
                                            const Acts::Surface& surface,
                                            const Acts::Transform3& ref) {
  // The xml to be translated
  std::stringstream sxml;
  auto boundValues = surface.bounds().values();

  std::array<int, 2u> axes = {0, 1};
  // Change/adapt the behavior
  switch (surface.bounds().type()) {
    case Acts::SurfaceBounds::eRectangle: {
      sxml << "<box ";
      double dx = (boundValues[2u] - boundValues[0u]);
      double dy = (boundValues[3u] - boundValues[1u]);
      double dz = 0.125;
      sxml << "dx=\"" << dx << "*mm\" ";
      sxml << "dy=\"" << dy << "*mm\" ";
      sxml << "dz=\"" << dz << "*mm\" ";
    }; break;
    case Acts::SurfaceBounds::eTrapezoid: {
      axes = {2, 0};

      sxml << "<trap ";
      double hxmin = boundValues[0u];
      double hxmax = boundValues[1u];
      double dy = 2 * boundValues[2u];
      double dz = 0.125;
      sxml << "x1=\"" << hxmin << "*mm\" ";
      sxml << "x2=\"" << hxmax << "*mm\" ";
      sxml << "dy=\"" << dy << "*mm\" ";
      sxml << "dz=\"" << dz << "*mm\" ";
    }; break;
    default:
      break;
  }

  // Unwind the placement you have already
  auto relTransform = ref * surface.transform(gctx);
  sxml << transformToXML(relTransform, axes);
  sxml << " material=\"Air\"";
  sxml << " sensitive=\"true\"/>";
  return sxml.str();
}
