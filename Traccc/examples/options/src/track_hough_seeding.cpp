// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/options/track_hough_seeding.hpp"

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace traccc::opts {

/// Convenience namespace int16_thand
namespace po = boost::program_options;

static std::string vector_as_text(const auto& values) {
  std::ostringstream out;
  const char* separator = "";

  for (const auto& value : values) {
    out << separator << value;
    separator = " ";
  }

  return out.str();
}

track_hough_seeding::track_hough_seeding() : interface("HOUGH Options") {
  m_desc.add_options()("useHough",
                       boost::program_options::bool_switch(&useHough),
                       "use hough algorithm");

  m_desc.add_options()(
      "n_wedges",
      po::value(&hough_config.n_wedges)->default_value(hough_config.n_wedges),
      "Number of wedges");

  m_desc.add_options()(
      "n_cols",
      po::value(&hough_config.n_cols)->default_value(hough_config.n_cols),
      "Number of columns");

  m_desc.add_options()(
      "n_rows",
      po::value(&hough_config.n_rows)->default_value(hough_config.n_rows),
      "Number of rows");

  m_desc.add_options()("min_range_x",
                       po::value(&hough_config.min_range_x)
                           ->default_value(hough_config.min_range_x),
                       "Minimum x range");

  m_desc.add_options()("max_range_x",
                       po::value(&hough_config.max_range_x)
                           ->default_value(hough_config.max_range_x),
                       "Maximum x range");

  m_desc.add_options()("min_range_y",
                       po::value(&hough_config.min_range_y)
                           ->default_value(hough_config.min_range_y),
                       "Minimum y range");

  m_desc.add_options()("max_range_y",
                       po::value(&hough_config.max_range_y)
                           ->default_value(hough_config.max_range_y),
                       "Maximum y range");

  m_desc.add_options()("interaction_center",
                       po::value(&hough_config.interaction_center)
                           ->default_value(hough_config.interaction_center),
                       "Interaction center");

  m_desc.add_options()("interaction_window",
                       po::value(&hough_config.interaction_window)
                           ->default_value(hough_config.interaction_window),
                       "Interaction window");

  m_desc.add_options()(
      "eta_width",
      po::value(&hough_config.eta_width)->default_value(hough_config.eta_width),
      "Eta width");

  m_desc.add_options()(
      "z_window",
      po::value(&hough_config.z_window)->default_value(hough_config.z_window),
      "Z window");

  m_desc.add_options()("wedges",
                       po::value(&hough_config.wedges)
                           ->multitoken()
                           ->default_value(hough_config.wedges,
                                           vector_as_text(hough_config.wedges)),
                       "Wedge values (space-separated)");

  m_desc.add_options()(
      "min_contributions",
      po::value(&hough_config.min_contributions)
          ->multitoken()
          ->default_value(hough_config.min_contributions,
                          vector_as_text(hough_config.min_contributions)),
      "Minimum contributions per wedge (space-separated)");

  m_desc.add_options()("module_veto",
                       po::value(&hough_config.module_veto)
                           ->default_value(hough_config.module_veto),
                       "Module veto bitmask (decimal)");

  m_desc.add_options()(
      "shift_map",
      po::value(&hough_config.shift_map)
          ->multitoken()
          ->default_value(hough_config.shift_map,
                          vector_as_text(hough_config.shift_map)),
      "Bit position per detector volume");

  m_desc.add_options()("nn_model_path",
                       po::value(&hough_config.nn_model_path)
                           ->default_value(hough_config.nn_model_path),
                       "ONNX model of the peak filter (empty: keep all peaks)");
}

track_hough_seeding::operator hough_seeding_config() const {
  return hough_config;
}

void track_hough_seeding::read(const boost::program_options::variables_map&) {
  if (!useHough) {
    return;
  }

  if (hough_config.n_wedges == 0u || hough_config.n_cols == 0u ||
      hough_config.n_rows == 0u) {
    throw std::runtime_error("n_wedges, n_cols, and n_rows must be positive");
  }

  if (hough_config.wedges.size() != hough_config.n_wedges) {
    throw std::runtime_error("wedges must contain n_wedges values");
  }

  if (hough_config.min_contributions.size() != hough_config.n_wedges) {
    throw std::runtime_error("min_contributions must contain n_wedges values");
  }

  if (!(hough_config.min_range_x < hough_config.max_range_x)) {
    throw std::runtime_error("min_range_x must be less than max_range_x");
  }

  if (!(hough_config.min_range_y < hough_config.max_range_y)) {
    throw std::runtime_error("min_range_y must be less than max_range_y");
  }

  for (const unsigned int shift : hough_config.shift_map) {
    if (shift != hough_config.invalid && shift > 64) {
      throw std::runtime_error("shift_map entries must be between 0 and 64");
    }
  }
}

std::unique_ptr<configuration_printable> track_hough_seeding::as_printable()
    const {
  auto cat = std::make_unique<configuration_category>(m_description);

  cat->add_child(std::make_unique<configuration_kv_pair>(
      "using hough algorithm ", std::to_string(useHough)));

  cat->add_child(std::make_unique<configuration_kv_pair>(
      "n_wedges", std::to_string(hough_config.n_wedges)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "n_cols", std::to_string(hough_config.n_cols)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "n_rows", std::to_string(hough_config.n_rows)));

  cat->add_child(std::make_unique<configuration_kv_pair>(
      "min_range_x", std::format("{:.5f}", hough_config.min_range_x)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "max_range_x", std::format("{:.5f}", hough_config.max_range_x)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "min_range_y", std::format("{:.5f}", hough_config.min_range_y)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "max_range_y", std::format("{:.5f}", hough_config.max_range_y)));

  cat->add_child(std::make_unique<configuration_kv_pair>(
      "interaction_center",
      std::format("{:.5f}", hough_config.interaction_center)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "interaction_window",
      std::format("{:.5f}", hough_config.interaction_window)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "eta_width", std::format("{:.5f}", hough_config.eta_width)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "z_window", std::format("{:.5f}", hough_config.z_window)));

  cat->add_child(std::make_unique<configuration_kv_pair>(
      "wedges", vector_as_text(hough_config.wedges)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "min_contributions", vector_as_text(hough_config.min_contributions)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "module_veto", std::to_string(hough_config.module_veto)));

  cat->add_child(std::make_unique<configuration_kv_pair>(
      "shift_map", vector_as_text(hough_config.shift_map)));

  cat->add_child(std::make_unique<configuration_kv_pair>(
      "nn_model_path", hough_config.nn_model_path));

  return cat;
}

}  // namespace traccc::opts
