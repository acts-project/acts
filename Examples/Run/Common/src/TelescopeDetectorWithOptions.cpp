// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/TelescopeDetectorWithOptions.hpp"

#include "ActsExamples/Utilities/Options.hpp"

#include <boost/program_options.hpp>

namespace ActsExamples {

void TelescopeDetectorWithOptions::addOptions(
    ActsExamples::Options::Description& desc) const {
  using boost::program_options::value;
  using namespace Options;

  auto opt = desc.add_options();
  opt("geo-tele-positions",
      value<VariableReals>()->default_value({{0, 30, 60, 120, 150, 180}}),
      "Telescope detector Input: the layers positions in the longidutinal "
      "direction in mm");
  opt("geo-tele-offsets", value<Reals<2>>()->default_value({{0, 0}}),
      "Telescope detector Input: the layers offsets in the transverse plane "
      "in "
      "mm. Same values for "
      "all layers");
  opt("geo-tele-bounds", value<Reals<2>>()->default_value({{25, 100}}),
      "Telescope detector Input: the values for surface bounds in mm: "
      "(halfX, halfY) - plane surface, (minR, maxR) - disc surface");
  opt("geo-tele-thickness", value<double>()->default_value(80),
      "Telescope detector Input: the silicon material thickness of "
      "each layer in um. Same value for all layers");
  opt("geo-tele-surface", value<int>()->default_value(0),
      "Telescope detector Input: the detector surface type: 0 - plane "
      "surface, "
      "1 - disc surface");
  opt("geo-tele-alignaxis", value<int>()->default_value(2),
      "Telescope detector Input: the detector is placed along which "
      "axis: 0 - x axis, 1 - y axis, 2 - z aixs");
}

auto TelescopeDetectorWithOptions::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> /*mdecorator*/)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  Telescope::TelescopeDetector::Config cfg;

  cfg.positions = vm["geo-tele-positions"]
                      .template as<ActsExamples::Options::VariableReals>()
                      .values;
  cfg.offsets =
      vm["geo-tele-offsets"].template as<ActsExamples::Options::Reals<2>>();
  // The bounds values are taken as (halfX, halfY) for plane surface and
  // (minR, maxR) for disc surface
  cfg.bounds =
      vm["geo-tele-bounds"].template as<ActsExamples::Options::Reals<2>>();
  // Translate the thickness in unit of mm
  cfg.thickness = vm["geo-tele-thickness"].template as<double>() * 0.001;
  cfg.surfaceType = vm["geo-tele-surface"].template as<int>();
  cfg.binValue = vm["geo-tele-alignaxis"].template as<int>();

  return m_detector.finalize(cfg, {});
}

}  // namespace ActsExamples
