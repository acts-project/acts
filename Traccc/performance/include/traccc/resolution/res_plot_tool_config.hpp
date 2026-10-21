// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/utils/helpers.hpp"

// System include(s).
#include <map>
#include <string>
#include <vector>

namespace traccc {

/// @brief Configuration structure for @c traccc::res_plot_tool
struct res_plot_tool_config {
  /// parameter sets to do plots
  std::vector<std::string> param_names = {"d0",  "z0", "phi",  "theta",
                                          "qop", "t",  "qopT", "qopz"};

  /// Binning setups
  std::map<std::string, plot_helpers::binning> var_binning = {
      {"pull", plot_helpers::binning("pull", 100, -5, 5)},
      {"residual_d0", plot_helpers::binning("r_{d0} [mm]", 100, -0.5f, 0.5f)},
      {"residual_z0", plot_helpers::binning("r_{z0} [mm]", 100, -0.5f, 0.5f)},
      {"residual_phi",
       plot_helpers::binning("r_{#phi} [rad]", 100, -0.01f, 0.01f)},
      {"residual_theta",
       plot_helpers::binning("r_{#theta} [rad]", 100, -0.01f, 0.01f)},
      {"residual_qop",
       plot_helpers::binning("r_{q/p} [c/GeV]", 100, -1.f, 1.f)},
      {"residual_t", plot_helpers::binning("r_{t} [s]", 100, -1000.f, 1000.f)},
      {"residual_qopT",
       plot_helpers::binning("r_{q/p_{T}} [c/GeV]", 100, -1.f, 1.f)},
      {"residual_qopz",
       plot_helpers::binning("r_{q/p_{z}} [c/GeV]", 100, -1.f, 1.f)},
      {"resolution_d0",
       plot_helpers::binning("#sigma_{d0} [mm]", 100, 0.f, 0.5f)},
      {"resolution_z0",
       plot_helpers::binning("#sigma_{z0} [mm]", 100, 0.f, 0.5f)},
      {"resolution_phi",
       plot_helpers::binning("#sigma_{#phi} [rad]", 100, 0.f, 0.01f)},
      {"resolution_theta",
       plot_helpers::binning("#sigma_{#theta} [rad]", 100, 0.f, 0.01f)},
      {"resolution_qop",
       plot_helpers::binning("#sigma_{q/p} [c/GeV]", 100, 0.f, 0.1f)},
      {"resolution_t",
       plot_helpers::binning("#sigma_{t} [s]", 100, -1000.f, 1000.f)},
      {"resolution_qopT",
       plot_helpers::binning("#sigma_{q/p_{T}} [c/GeV]", 100, 0.f, 0.1f)},
      {"resolution_qopz",
       plot_helpers::binning("#sigma_{q/p_{z}} [c/GeV]", 100, 0.f, 0.1f)},
      {"Eta", plot_helpers::binning("#eta", 40, -4.f, 4.f)},
      {"Pt", plot_helpers::binning("p_{T} [GeV/c]", 40, 0.f, 100.f)},
  };
};

}  // namespace traccc
