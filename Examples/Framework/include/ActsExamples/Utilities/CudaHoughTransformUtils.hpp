// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#ifdef ACTS_ENABLE_CUDA


#include "Acts/Seeding/HoughTransformUtils.hpp"
#include "ActsExamples/EventData/CudaMuonSpacePoint.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace ActsExamples::CudaHoughTransformUtils {

using YieldType = Acts::HoughTransformUtils::YieldType;
using CoordType = Acts::HoughTransformUtils::CoordType;
using HoughPlaneConfig = Acts::HoughTransformUtils::HoughPlaneConfig;
using HoughAxisRanges = Acts::HoughTransformUtils::HoughAxisRanges;

/// Bit mask encoding which logical detector layers contributed to one Hough
/// cell. 
///
/// One bit corresponds to one zero-based logical layer:
///
///   layer 0 -> bit 0 -> 00000001
using LayerMask = unsigned long long;

/// Non-owning device-side event-level batch of Hough planes.
struct CudaHoughPlaneBatchArrays {
  /// Weighted hit contribution per bucket/cell.
  YieldType* nHits = nullptr;

  /// Weighted unique-layer contribution per bucket/cell.
  YieldType* nLayers = nullptr;

  /// Bit mask of logical layers seen by each bucket/cell.
  LayerMask* layerMask = nullptr;

  std::uint32_t nBuckets = 0;
  std::uint32_t nBinsX = 0;
  std::uint32_t nBinsY = 0;
};

/// Host-side snapshot of one Hough cell.
class CudaHoughCell {
 public:
  CudaHoughCell() = default;

  CudaHoughCell(YieldType nHits, YieldType nLayers,
                LayerMask layerMask) noexcept
      : m_nHits{nHits}, m_nLayers{nLayers}, m_layerMask{layerMask} {}

  YieldType nHits() const noexcept { return m_nHits; }
  YieldType nLayers() const noexcept { return m_nLayers; }
  LayerMask layerMask() const noexcept { return m_layerMask; }

  bool hasLayer(unsigned layer) const noexcept;

 private:
  YieldType m_nHits = 0.0f;
  YieldType m_nLayers = 0.0f;
  LayerMask m_layerMask = 0ull;
};

/// Event-level CUDA Hough accumulator batch.
///
/// This class intentionally does not store per-cell hit identifiers. It stores:
///
///   nHits[bucket, cell]
///   nLayers[bucket, cell]
///   layerMask[bucket, cell]
///
/// Hit association is done later after peak finding, since otherwise there would be need for large prealocation.
class CudaHoughPlaneBatch {
 public:
  using size_type = std::size_t;
  using Index = std::array<size_type, 2>;

  CudaHoughPlaneBatch(const HoughPlaneConfig& cfg, size_type nBuckets);

  CudaHoughPlaneBatch(const CudaHoughPlaneBatch&) = delete;
  CudaHoughPlaneBatch& operator=(const CudaHoughPlaneBatch&) = delete;

  CudaHoughPlaneBatch(CudaHoughPlaneBatch&& other) noexcept;
  CudaHoughPlaneBatch& operator=(CudaHoughPlaneBatch&& other) noexcept;

  ~CudaHoughPlaneBatch() noexcept;

  size_type nBuckets() const noexcept { return m_nBuckets; }
  size_type nBinsX() const noexcept { return m_cfg.nBinsX; }
  size_type nBinsY() const noexcept { return m_cfg.nBinsY; }
  size_type nCellsPerBucket() const noexcept { return nBinsX() * nBinsY(); }
  size_type totalCells() const noexcept {
    return nBuckets() * nCellsPerBucket();
  }

  bool empty() const noexcept { return totalCells() == 0; }

  /// Row-major flat index inside the whole batch:
  ///
  ///   globalBin = bucket * nCellsPerBucket + yBin * nBinsX + xBin
  size_type globalBin(size_type bucket, size_type xBin, size_type yBin) const;

  /// Reverse mapping from global batch bin to {xBin, yBin}.
  Index axisBins(size_type globalBin) const;

  /// CPU-side direct bin fill. Useful for vdalidation.
  void fillBin(size_type bucket, size_type xBin, size_type yBin,
               unsigned layer, YieldType weight = 1.0f);

  /// CPU reference fill for MDT eta Hough, all buckets in the event.
  void fillEtaDriftCirclesHost(const CudaMuonSpacePointContainer& spacePoints,
                               const HoughAxisRanges& axisRanges,
                               double widthScale = 3.0,
                               double maxWidth = 1.0,
                               YieldType weight = 1.0f);

void fillEtaDriftCirclesOnDevice(CudaMuonSpacePointContainer& spacePoints,
                                 const HoughAxisRanges& axisRanges,
                                 double widthScale = 3.0,
                                 double maxWidth = 1.0,
                                 YieldType weight = 1.0f,
                                 std::uint32_t threadsPerBlock = 128,
                                 std::uint32_t num_blocks = 0); // 0 is auto use number of SMs

  /// Reset host data and device data if allocated.
  void reset();

  YieldType nHits(size_type bucket, size_type xBin, size_type yBin) const;
  YieldType nLayers(size_type bucket, size_type xBin, size_type yBin) const;
  LayerMask layerMask(size_type bucket, size_type xBin, size_type yBin) const;

  bool hasLayer(size_type bucket, size_type xBin, size_type yBin,
                unsigned layer) const;

  std::vector<unsigned> layers(size_type bucket, size_type xBin,
                               size_type yBin) const;

  CudaHoughCell cell(size_type bucket, size_type xBin, size_type yBin) const;

  YieldType maxHits(size_type bucket) const;
  YieldType maxLayers(size_type bucket) const;

  Index locMaxHits(size_type bucket) const;
  Index locMaxLayers(size_type bucket) const;

  std::vector<size_type> nonEmptyBins(size_type bucket) const;

  void moveToDevice();
  void moveToHost();
  void clearDevice() noexcept;

  bool isOnDevice() const noexcept { return m_onDevice; }

  CudaHoughPlaneBatchArrays deviceArrays() const noexcept { return m_device; }

 private:
  HoughPlaneConfig m_cfg{};
  size_type m_nBuckets = 0;

  std::vector<YieldType> m_hostHits{};
  std::vector<YieldType> m_hostLayers{};
  std::vector<LayerMask> m_hostLayerMask{};

  CudaHoughPlaneBatchArrays m_device{};
  bool m_onDevice = false;

  size_type uncheckedGlobalBin(size_type bucket, size_type xBin,
                               size_type yBin) const noexcept {
    return bucket * nCellsPerBucket() + yBin * nBinsX() + xBin;
  }

  void checkBucket(size_type bucket) const;
  void checkIndices(size_type xBin, size_type yBin) const;
  void checkGlobalBin(size_type globalBin) const;
  void checkSpacePointBuckets(
      const CudaMuonSpacePointContainer& spacePoints) const;
};

}  // namespace ActsExamples::CudaHoughTransformUtils

#endif
