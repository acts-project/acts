// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"

namespace Acts::Experimental {

/// Collection of pointers to the space point container and its
/// additional columns. This is used as a basket to pass around
/// the input data for the triplet seed finder.
class SpacePointContainerPointers {
 public:
  /// Minimal input: space points and r column.
  SpacePointContainerPointers(
      const SpacePointContainer2& spacePoints,
      const SpacePointContainer2::DenseColumn<float>& rColumn)
      : m_spacePoints(&spacePoints), m_rColumn(&rColumn) {}

  /// Space points, r column, and variance columns.
  SpacePointContainerPointers(
      const SpacePointContainer2& spacePoints,
      const SpacePointContainer2::DenseColumn<float>& rColumn,
      const SpacePointContainer2::DenseColumn<float>& varianceRColumn,
      const SpacePointContainer2::DenseColumn<float>& varianceZColumn)
      : m_spacePoints(&spacePoints),
        m_rColumn(&rColumn),
        m_varianceRColumn(&varianceRColumn),
        m_varianceZColumn(&varianceZColumn) {}

  /// Space points, r column, variance columns, and strip columns.
  SpacePointContainerPointers(
      const SpacePointContainer2& spacePoints,
      const SpacePointContainer2::DenseColumn<float>& rColumn,
      const SpacePointContainer2::DenseColumn<float>& varianceRColumn,
      const SpacePointContainer2::DenseColumn<float>& varianceZColumn,
      const SpacePointContainer2::DenseColumn<Vector3>& topStripVectorColumn,
      const SpacePointContainer2::DenseColumn<Vector3>& bottomStripVectorColumn,
      const SpacePointContainer2::DenseColumn<Vector3>&
          stripCenterDistanceColumn,
      const SpacePointContainer2::DenseColumn<Vector3>&
          topStripCenterPositionColumn)
      : m_spacePoints(&spacePoints),
        m_rColumn(&rColumn),
        m_varianceRColumn(&varianceRColumn),
        m_varianceZColumn(&varianceZColumn),
        m_topStripVectorColumn(&topStripVectorColumn),
        m_bottomStripVectorColumn(&bottomStripVectorColumn),
        m_stripCenterDistanceColumn(&stripCenterDistanceColumn),
        m_topStripCenterPositionColumn(&topStripCenterPositionColumn) {}

  /// Pointer to the copied-from index column, if available.
  const SpacePointContainer2::DenseColumn<SpacePointIndex2>*
      copiedFromIndexColumn = nullptr;

  [[nodiscard]] const SpacePointContainer2& spacePoints() const {
    return *m_spacePoints;
  }
  [[nodiscard]] const SpacePointContainer2::DenseColumn<float>& rColumn()
      const {
    return *m_rColumn;
  }
  [[nodiscard]] const SpacePointContainer2::DenseColumn<float>&
  varianceRColumn() const {
    return *m_varianceRColumn;
  }
  [[nodiscard]] const SpacePointContainer2::DenseColumn<float>&
  varianceZColumn() const {
    return *m_varianceZColumn;
  }
  [[nodiscard]] const SpacePointContainer2::DenseColumn<Vector3>&
  topStripVectorColumn() const {
    return *m_topStripVectorColumn;
  }
  [[nodiscard]] const SpacePointContainer2::DenseColumn<Vector3>&
  bottomStripVectorColumn() const {
    return *m_bottomStripVectorColumn;
  }
  [[nodiscard]] const SpacePointContainer2::DenseColumn<Vector3>&
  stripCenterDistanceColumn() const {
    return *m_stripCenterDistanceColumn;
  }
  [[nodiscard]] const SpacePointContainer2::DenseColumn<Vector3>&
  topStripCenterPositionColumn() const {
    return *m_topStripCenterPositionColumn;
  }

  [[nodiscard]] bool hasVarianceColumns() const {
    return m_varianceRColumn != nullptr && m_varianceZColumn != nullptr;
  }
  [[nodiscard]] bool hasStripColumns() const {
    return m_topStripVectorColumn != nullptr &&
           m_bottomStripVectorColumn != nullptr &&
           m_stripCenterDistanceColumn != nullptr &&
           m_topStripCenterPositionColumn != nullptr;
  }
  [[nodiscard]] bool hasCopiedFromIndexColumn() const {
    return copiedFromIndexColumn != nullptr;
  }

 private:
  const SpacePointContainer2* m_spacePoints = nullptr;
  const SpacePointContainer2::DenseColumn<float>* m_rColumn = nullptr;

  const SpacePointContainer2::DenseColumn<float>* m_varianceRColumn = nullptr;
  const SpacePointContainer2::DenseColumn<float>* m_varianceZColumn = nullptr;

  const SpacePointContainer2::DenseColumn<Vector3>* m_topStripVectorColumn =
      nullptr;
  const SpacePointContainer2::DenseColumn<Vector3>* m_bottomStripVectorColumn =
      nullptr;
  const SpacePointContainer2::DenseColumn<Vector3>*
      m_stripCenterDistanceColumn = nullptr;
  const SpacePointContainer2::DenseColumn<Vector3>*
      m_topStripCenterPositionColumn = nullptr;
};

}  // namespace Acts::Experimental
