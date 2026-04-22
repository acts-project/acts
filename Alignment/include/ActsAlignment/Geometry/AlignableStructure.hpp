// This file is part of the ACTS project.
//
// Copyright (C) 2026 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// #include "Acts/Definitions/Alignment.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <map>
#include <memory>
#include <vector>

namespace ActsAlignment {

/**
 * @brief A class that groups multiple surfaces to be aligned together.
 *
 * AlignableStructure acts as a logical parent for a set of modules (surfaces).
 * It manages the 6 standard degrees of freedom (Center/Rotation) and
 * constraints.
 */
class AlignableStructure {
 public:
  /**
   * @brief Constructor
   * @param id The geometric identifier for this structure
   */
  explicit AlignableStructure(Acts::GeometryIdentifier id)
      : m_id(id), m_surfaces(), m_constraints(), m_mask(AlignmentMask::None) {}

  /// Disallow copy construction to ensure identity uniqueness
  AlignableStructure(const AlignableStructure&) = delete;
  AlignableStructure& operator=(const AlignableStructure&) = delete;

  /**
   * @brief Add a surface to this alignable structure
   * @param surface The surface to add
   */
  void addSurface(std::shared_ptr<Acts::Surface> surface) {
    m_surfaces.push_back(std::move(surface));
  }

  /**
   * @brief Get the geometric identifier
   * @return The ID of the structure
   */
  Acts::GeometryIdentifier geometryId() const { return m_id; }

  /**
   * @brief Access the list of associated surfaces
   * @return A vector of shared pointers to surfaces
   */
  const std::vector<std::shared_ptr<Acts::Surface>>& surfaces() const {
    return m_surfaces;
  }

  /**
   * @brief Access the constraints map (Rigidity)
   *
   * The map key is the standard alignment index.
   * The value is the variance constraint (W matrix entry).
   * @return Reference to the constraints map
   */
  std::map<Acts::AlignmentIndices, double>& constraints() {
    return m_constraints;
  }
  const std::map<Acts::AlignmentIndices, double>& constraints() const {
    return m_constraints;
  }

  /**
   * @brief Access the alignment mask (floating DoFs)
   * @return Reference to the mask
   */
  AlignmentMask& alignmentMask() { return m_mask; }
  AlignmentMask alignmentMask() const { return m_mask; }

 private:
  /// The unique identifier for this structure
  Acts::GeometryIdentifier m_id;

  /// The collection of sensors governed by this structure
  std::vector<std::shared_ptr<Acts::Surface>> m_surfaces;

  /// Map of parameter index -> constraint value (rigidity)
  std::map<Acts::AlignmentIndices, double> m_constraints;

  /// Bitmask for determining which of the 6 DoFs are floating
  AlignmentMask m_mask{AlignmentMask::None};
};

}  // namespace ActsAlignment
