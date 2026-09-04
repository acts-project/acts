// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Detray test include(s)
#include "detray/test/framework/test_configuration.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <limits>
#include <memory>
#include <string>

namespace detray::test {

/// @brief Configuration for a detector scan test.
template <concepts::algebra algebra_t>
struct material_validation_config
    : public detray::test::configuration<dscalar<algebra_t>> {
  using scalar_type = dscalar<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using base_type = detray::test::configuration<scalar_type>;

  /// Name of the test
  std::string m_name{"material_validation"};
  /// Vecmem memory resource for the device allocations
  vecmem::memory_resource *m_dev_mr{nullptr};
  /// Name of the output file, containing the complete ray material traces
  std::string m_material_file{"navigation_material_trace"};
  /// The maximal number of test tracks to run
  std::size_t m_n_tracks{detray::detail::invalid_value<std::size_t>()};
  /// Allowed relative discrepancy between truth and navigation material
  scalar_type m_rel_error{0.001f};

  /// Getters
  /// @{
  const std::string &name() const { return m_name; }
  vecmem::memory_resource *device_mr() const { return m_dev_mr; }
  const std::string &material_file() const { return m_material_file; }
  std::size_t n_tracks() const { return m_n_tracks; }
  scalar_type relative_error() const { return m_rel_error; }
  /// @}

  /// Setters
  /// @{
  material_validation_config &name(const std::string &n) {
    m_name = n;
    return *this;
  }
  material_validation_config &device_mr(vecmem::memory_resource *mr) {
    m_dev_mr = mr;
    return *this;
  }
  material_validation_config &material_file(const std::string &f) {
    m_material_file = f;
    return *this;
  }
  material_validation_config &n_tracks(std::size_t n) {
    m_n_tracks = n;
    return *this;
  }
  material_validation_config &relative_error(scalar_type re) {
    assert(re > 0.f);
    m_rel_error = re;
    return *this;
  }
  /// @}
};

}  // namespace detray::test
