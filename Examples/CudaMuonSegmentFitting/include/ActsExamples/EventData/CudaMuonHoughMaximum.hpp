// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/HoughTransformUtils.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include <cuda_runtime.h>

namespace ActsExamples {

using YieldType = Acts::HoughTransformUtils::YieldType;
using CoordType = Acts::HoughTransformUtils::CoordType;
using LayerMask = unsigned long long;

/// @brief Non-owning device-side view of a batch of Hough maxima.
///
/// Maximum slots are stored bucket-by-bucket:
///
///   bucket 0: maximum 0, maximum 1, ..., maximum N - 1
///   bucket 1: maximum 0, maximum 1, ..., maximum N - 1
struct CudaHoughMaximumBatchArrays {
  CoordType* tanBeta = nullptr;
  CoordType* interceptY = nullptr;

  YieldType* nHits = nullptr;
  YieldType* nLayers = nullptr;
  LayerMask* layerMask = nullptr;

  std::uint32_t* xBin = nullptr;
  std::uint32_t* yBin = nullptr;

  /// Number of occupied maximum slots in each bucket.
  std::uint32_t* nMaxima = nullptr;

  /// Number of input space points associated with each maximum.
  std::uint32_t* nAssociatedHits = nullptr;

  /// CSR offsets into associatedHitIndices.
  /// Maximum slot i owns [associatedHitOffsets[i], associatedHitOffsets[i+1])
  /// Array as max capacity + 1
  std::uint32_t* associatedHitOffsets = nullptr;

  /// Flat list of indices into CudaMuonSpacePointContainer.
  std::uint32_t* associatedHitIndices = nullptr;

  std::uint32_t nBuckets = 0;
  std::uint32_t capacityPerBucket = 0;
  /// Full size of the associated-hit list.
  std::uint32_t totalAssociatedHits = 0;

  /// Return the flat array index for a bucket and maximum slot.
  __host__ __device__ std::uint32_t index(
      std::uint32_t bucket, std::uint32_t maximum) const noexcept {
    return bucket * capacityPerBucket + maximum;
  }
};

/// @brief Owning host/device storage for Hough-maximum slots per bucket.
///
/// The capacity per bucket specifies the maximum number of maxima that may be
/// stored for each bucket. The actual number stored is tracked independently
/// by the nMaxima counter for each bucket.
class CudaHoughMaximumBatch {
 public:
  using size_type = std::size_t;

  CudaHoughMaximumBatch(size_type nBuckets, size_type capacityPerBucket);

  /// Copy creates problem with cuda memory ownership
  CudaHoughMaximumBatch(const CudaHoughMaximumBatch&) = delete;
  CudaHoughMaximumBatch& operator=(const CudaHoughMaximumBatch&) = delete;

  CudaHoughMaximumBatch(CudaHoughMaximumBatch&& other) noexcept;
  CudaHoughMaximumBatch& operator=(CudaHoughMaximumBatch&& other) noexcept;

  ~CudaHoughMaximumBatch() noexcept;

  size_type capacityPerBucket() const noexcept { return m_capacityPerBucket; }

  size_type nBuckets() const noexcept { return m_nBuckets; }

  size_type totalCapacity() const noexcept {
    return nBuckets() * capacityPerBucket();
  }

  bool empty() const noexcept { return nBuckets() == 0; }

  /// Return the number of stored maxima in one bucket.
  size_type nMaxima(size_type bucket) const;

  CoordType tanBeta(size_type bucket, size_type maximum) const;
  CoordType interceptY(size_type bucket, size_type maximum) const;

  YieldType nHits(size_type bucket, size_type maximum) const;
  YieldType nLayers(size_type bucket, size_type maximum) const;

  LayerMask layerMask(size_type bucket, size_type maximum) const;

  size_type xBin(size_type bucket, size_type maximum) const;
  size_type yBin(size_type bucket, size_type maximum) const;

  /// Copy only the metadata required to allocate hit-association storage.
  ///
  /// Copies nMaxima and nAssociatedHits from the device.
  void copyAssociationMetadataToHost(cudaStream_t stream);

  /// Calculate CSR offsets from the copied association counts and allocate
  /// exact associated-hit storage on the device.
  void allocateAssociationStorage(cudaStream_t stream);

  /// Copy the associated space-point indices from the device to the host.
  ///
  /// This is called after the association-fill kernel has completed.
  void copyAssociatedHitIndicesToHost(cudaStream_t stream);

  /// Return the number of associated input space points for one maximum.
  size_type nAssociatedHits(size_type bucket, size_type maximum) const;

  /// Return all original space-point indices associated with one maximum.
  std::span<const std::uint32_t> associatedHitIndices(size_type bucket,
                                                      size_type maximum) const;

  /// Return the total exact size of associatedHitIndices.
  size_type totalAssociatedHits() const noexcept {
    return m_hostAssociatedHitIndices.size();
  }

  bool hasAssociationStorage() const noexcept {
    return m_associationStorageAllocated;
  }

  /// Allocate device output storage for use on the supplied stream.
  void moveToDevice(cudaStream_t stream);

  /// Copy device data back into host storage.
  void moveToHost(cudaStream_t stream);

  /// Release all owned device storage.
  void clearDevice() noexcept;

  bool isOnDevice() const noexcept { return m_onDevice; }

  CudaHoughMaximumBatchArrays deviceArrays() const noexcept { return m_device; }

 private:
  size_type m_nBuckets = 0;
  size_type m_capacityPerBucket = 0;

  std::vector<CoordType> m_hostTanBeta{};
  std::vector<CoordType> m_hostInterceptY{};

  std::vector<YieldType> m_hostHits{};
  std::vector<YieldType> m_hostLayers{};
  std::vector<LayerMask> m_hostLayerMask{};

  std::vector<std::uint32_t> m_hostXBin{};
  std::vector<std::uint32_t> m_hostYBin{};

  /// One counter per bucket, not one counter per maximum slot.
  std::vector<std::uint32_t> m_hostNMaxima{};

  /// One exact association count per maximum slot.
  std::vector<std::uint32_t> m_hostNAssociatedHits{};

  /// CSR offsets. Size is totalCapacity() + 1 after allocation.
  std::vector<std::uint32_t> m_hostAssociatedHitOffsets{};

  /// Exact flat associated-hit storage.
  std::vector<std::uint32_t> m_hostAssociatedHitIndices{};

  /// Whether nMaxima and nAssociatedHits have been copied from the device.
  bool m_associationMetadataOnHost = false;

  /// Whether exact device-side offsets and index storage have been allocated.
  bool m_associationStorageAllocated = false;

  /// Whether the final associated indices have been copied back to the host.
  bool m_associatedHitIndicesOnHost = false;

  CudaHoughMaximumBatchArrays m_device{};
  bool m_onDevice = false;

  size_type slotIndex(size_type bucket, size_type maximum) const noexcept {
    return bucket * capacityPerBucket() + maximum;
  }

  void checkBucket(size_type bucket) const;
  void checkMaximum(size_type bucket, size_type maximum) const;

  /// Build host-side CSR offsets and resize the exact hit-index storage.
  void prepareAssociationStorageHost();

  /// Release only the variable-size association offset/index storage.
  void clearAssociationStorage() noexcept;
};

}  // namespace ActsExamples
