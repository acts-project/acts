// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Utilities/CudaUtilities.hpp"

#include <limits>
#include <stdexcept>
#include <utility>

#include <cuda_runtime_api.h>

namespace ActsExamples {

template <std::size_t MaximaPerBucket>
CudaHoughMaximumBatch<MaximaPerBucket>::CudaHoughMaximumBatch(
    size_type nBuckets)
    : m_nBuckets{nBuckets} {
  if (m_nBuckets == 0u) {
    throw std::invalid_argument(
        "CudaHoughMaximumBatch requires non-zero nBuckets");
  }

  constexpr size_type maxUint32 =
      static_cast<size_type>(std::numeric_limits<std::uint32_t>::max());

  if (m_nBuckets > maxUint32) {
    throw std::overflow_error(
        "CudaHoughMaximumBatch nBuckets must fit into std::uint32_t and be "
        "allocatable");
  }

  if (m_nBuckets > std::numeric_limits<size_type>::max() / MaximaPerBucket) {
    throw std::overflow_error(
        "CudaHoughMaximumBatch total capacity overflows std::size_t and be "
        "allocatable");
  }

  const size_type capacity = totalCapacity();

  if (capacity > maxUint32) {
    throw std::overflow_error(
        "CudaHoughMaximumBatch total capacity must fit into std::uint32_t and "
        "be allocatable");
  }

  m_hostTanBeta.resize(capacity, CoordType{0.0});
  m_hostInterceptY.resize(capacity, CoordType{0.0});

  m_hostHits.resize(capacity, YieldType{0.0});
  m_hostLayers.resize(capacity, YieldType{0.0});
  m_hostLayerMask.resize(capacity, LayerMask{0ull});

  m_hostXBin.resize(capacity, 0u);
  m_hostYBin.resize(capacity, 0u);

  m_hostNMaxima.resize(m_nBuckets, 0u);
  m_hostNAssociatedHits.resize(capacity, 0u);
}

template <std::size_t MaximaPerBucket>
CudaHoughMaximumBatch<MaximaPerBucket>::CudaHoughMaximumBatch(
    CudaHoughMaximumBatch&& other) noexcept
    : m_nBuckets{std::exchange(other.m_nBuckets, 0u)},
      m_hostTanBeta{std::move(other.m_hostTanBeta)},
      m_hostInterceptY{std::move(other.m_hostInterceptY)},
      m_hostHits{std::move(other.m_hostHits)},
      m_hostLayers{std::move(other.m_hostLayers)},
      m_hostLayerMask{std::move(other.m_hostLayerMask)},
      m_hostXBin{std::move(other.m_hostXBin)},
      m_hostYBin{std::move(other.m_hostYBin)},
      m_hostNMaxima{std::move(other.m_hostNMaxima)},
      m_hostNAssociatedHits{std::move(other.m_hostNAssociatedHits)},
      m_hostAssociatedHitOffsets{std::move(other.m_hostAssociatedHitOffsets)},
      m_hostAssociatedHitIndices{std::move(other.m_hostAssociatedHitIndices)},
      m_associationMetadataOnHost{
          std::exchange(other.m_associationMetadataOnHost, false)},
      m_associationStorageAllocated{
          std::exchange(other.m_associationStorageAllocated, false)},
      m_associatedHitIndicesOnHost{
          std::exchange(other.m_associatedHitIndicesOnHost, false)},
      m_device{std::exchange(other.m_device, CudaHoughMaximumBatchArrays{})},
      m_onDevice{std::exchange(other.m_onDevice, false)} {}

template <std::size_t MaximaPerBucket>
CudaHoughMaximumBatch<MaximaPerBucket>&
CudaHoughMaximumBatch<MaximaPerBucket>::operator=(
    CudaHoughMaximumBatch&& other) noexcept {
  if (this == &other) {
    return *this;
  }

  clearDevice();

  m_nBuckets = std::exchange(other.m_nBuckets, 0u);

  m_hostTanBeta = std::move(other.m_hostTanBeta);
  m_hostInterceptY = std::move(other.m_hostInterceptY);

  m_hostHits = std::move(other.m_hostHits);
  m_hostLayers = std::move(other.m_hostLayers);
  m_hostLayerMask = std::move(other.m_hostLayerMask);

  m_hostXBin = std::move(other.m_hostXBin);
  m_hostYBin = std::move(other.m_hostYBin);

  m_hostNMaxima = std::move(other.m_hostNMaxima);
  m_hostNAssociatedHits = std::move(other.m_hostNAssociatedHits);
  m_hostAssociatedHitOffsets = std::move(other.m_hostAssociatedHitOffsets);
  m_hostAssociatedHitIndices = std::move(other.m_hostAssociatedHitIndices);

  m_associationMetadataOnHost =
      std::exchange(other.m_associationMetadataOnHost, false);
  m_associationStorageAllocated =
      std::exchange(other.m_associationStorageAllocated, false);
  m_associatedHitIndicesOnHost =
      std::exchange(other.m_associatedHitIndicesOnHost, false);

  m_device = std::exchange(other.m_device, CudaHoughMaximumBatchArrays{});
  m_onDevice = std::exchange(other.m_onDevice, false);

  return *this;
}

template <std::size_t MaximaPerBucket>
CudaHoughMaximumBatch<MaximaPerBucket>::~CudaHoughMaximumBatch() noexcept {
  clearDevice();
}

template <std::size_t MaximaPerBucket>
typename CudaHoughMaximumBatch<MaximaPerBucket>::size_type
CudaHoughMaximumBatch<MaximaPerBucket>::nMaxima(size_type bucket) const {
  checkBucket(bucket);

  const size_type count = static_cast<size_type>(m_hostNMaxima[bucket]);

  if (count > capacityPerBucket()) {
    throw std::runtime_error(
        "CudaHoughMaximumBatch contains an invalid maximum count");
  }

  return count;
}

template <std::size_t MaximaPerBucket>
CoordType CudaHoughMaximumBatch<MaximaPerBucket>::tanBeta(
    size_type bucket, size_type maximum) const {
  checkMaximum(bucket, maximum);
  return m_hostTanBeta[slotIndex(bucket, maximum)];
}

template <std::size_t MaximaPerBucket>
CoordType CudaHoughMaximumBatch<MaximaPerBucket>::interceptY(
    size_type bucket, size_type maximum) const {
  checkMaximum(bucket, maximum);
  return m_hostInterceptY[slotIndex(bucket, maximum)];
}

template <std::size_t MaximaPerBucket>
YieldType CudaHoughMaximumBatch<MaximaPerBucket>::nHits(
    size_type bucket, size_type maximum) const {
  checkMaximum(bucket, maximum);
  return m_hostHits[slotIndex(bucket, maximum)];
}

template <std::size_t MaximaPerBucket>
YieldType CudaHoughMaximumBatch<MaximaPerBucket>::nLayers(
    size_type bucket, size_type maximum) const {
  checkMaximum(bucket, maximum);
  return m_hostLayers[slotIndex(bucket, maximum)];
}

template <std::size_t MaximaPerBucket>
LayerMask CudaHoughMaximumBatch<MaximaPerBucket>::layerMask(
    size_type bucket, size_type maximum) const {
  checkMaximum(bucket, maximum);
  return m_hostLayerMask[slotIndex(bucket, maximum)];
}

template <std::size_t MaximaPerBucket>
typename CudaHoughMaximumBatch<MaximaPerBucket>::size_type
CudaHoughMaximumBatch<MaximaPerBucket>::xBin(size_type bucket,
                                             size_type maximum) const {
  checkMaximum(bucket, maximum);
  return static_cast<size_type>(m_hostXBin[slotIndex(bucket, maximum)]);
}

template <std::size_t MaximaPerBucket>
typename CudaHoughMaximumBatch<MaximaPerBucket>::size_type
CudaHoughMaximumBatch<MaximaPerBucket>::yBin(size_type bucket,
                                             size_type maximum) const {
  checkMaximum(bucket, maximum);
  return static_cast<size_type>(m_hostYBin[slotIndex(bucket, maximum)]);
}

template <std::size_t MaximaPerBucket>
typename CudaHoughMaximumBatch<MaximaPerBucket>::size_type
CudaHoughMaximumBatch<MaximaPerBucket>::nAssociatedHits(
    size_type bucket, size_type maximum) const {
  if (!m_associationMetadataOnHost) {
    throw std::logic_error("Association metadata is not available on the host");
  }

  checkMaximum(bucket, maximum);
  return m_hostNAssociatedHits[slotIndex(bucket, maximum)];
}

template <std::size_t MaximaPerBucket>
std::span<const std::uint32_t>
CudaHoughMaximumBatch<MaximaPerBucket>::associatedHitIndices(
    size_type bucket, size_type maximum) const {
  if (!m_associatedHitIndicesOnHost) {
    throw std::logic_error(
        "Associated hit indices are not available on the host");
  }

  checkMaximum(bucket, maximum);

  const size_type slot = slotIndex(bucket, maximum);
  const size_type begin = m_hostAssociatedHitOffsets[slot];
  const size_type end = m_hostAssociatedHitOffsets[slot + 1u];

  return {m_hostAssociatedHitIndices.data() + begin, end - begin};
}

template <std::size_t MaximaPerBucket>
void CudaHoughMaximumBatch<MaximaPerBucket>::moveToDevice(
    cudaStream_t /*stream*/) {
  clearDevice();

  m_hostAssociatedHitOffsets.clear();
  m_hostAssociatedHitIndices.clear();

  m_associationMetadataOnHost = false;
  m_associationStorageAllocated = false;
  m_associatedHitIndicesOnHost = false;

  m_device.nBuckets = static_cast<std::uint32_t>(nBuckets());
  m_device.capacityPerBucket = static_cast<std::uint32_t>(capacityPerBucket());

  try {
    allocateDeviceColumn(m_device.tanBeta, totalCapacity());
    allocateDeviceColumn(m_device.interceptY, totalCapacity());

    allocateDeviceColumn(m_device.nHits, totalCapacity());
    allocateDeviceColumn(m_device.nLayers, totalCapacity());
    allocateDeviceColumn(m_device.layerMask, totalCapacity());

    allocateDeviceColumn(m_device.xBin, totalCapacity());
    allocateDeviceColumn(m_device.yBin, totalCapacity());

    allocateDeviceColumn(m_device.nMaxima, nBuckets());
    allocateDeviceColumn(m_device.nAssociatedHits, totalCapacity());
  } catch (...) {
    clearDevice();
    throw;
  }

  m_onDevice = true;
}

template <std::size_t MaximaPerBucket>
void CudaHoughMaximumBatch<MaximaPerBucket>::moveToHost(cudaStream_t stream) {
  if (!m_onDevice) {
    return;
  }

  copyColumnToHost(m_hostTanBeta, m_device.tanBeta, stream);
  copyColumnToHost(m_hostInterceptY, m_device.interceptY, stream);
  copyColumnToHost(m_hostHits, m_device.nHits, stream);
  copyColumnToHost(m_hostLayers, m_device.nLayers, stream);
  copyColumnToHost(m_hostLayerMask, m_device.layerMask, stream);
  copyColumnToHost(m_hostXBin, m_device.xBin, stream);
  copyColumnToHost(m_hostYBin, m_device.yBin, stream);
  copyColumnToHost(m_hostNMaxima, m_device.nMaxima, stream);
  copyColumnToHost(m_hostNAssociatedHits, m_device.nAssociatedHits, stream);
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  m_associationMetadataOnHost = true;
}

template <std::size_t MaximaPerBucket>
void CudaHoughMaximumBatch<MaximaPerBucket>::copyAssociationMetadataToHost(
    cudaStream_t stream) {
  if (!m_onDevice) {
    throw std::logic_error("CudaHoughMaximumBatch is not on the device");
  }

  copyColumnToHost(m_hostNMaxima, m_device.nMaxima, stream);
  copyColumnToHost(m_hostNAssociatedHits, m_device.nAssociatedHits, stream);
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  m_associationMetadataOnHost = true;
}

template <std::size_t MaximaPerBucket>
void CudaHoughMaximumBatch<MaximaPerBucket>::prepareAssociationStorageHost() {
  clearAssociationStorage();
  m_hostAssociatedHitOffsets.assign(totalCapacity() + 1u, 0u);
  m_hostAssociatedHitIndices.clear();

  std::uint64_t totalAssociatedHits = 0u;
  for (size_type bucket = 0u; bucket < nBuckets(); ++bucket) {
    const size_type maximaInBucket = nMaxima(bucket);
    for (size_type maximum = 0u; maximum < capacityPerBucket(); ++maximum) {
      const size_type slot = slotIndex(bucket, maximum);
      const std::uint32_t count = m_hostNAssociatedHits[slot];
      if (maximum >= maximaInBucket && count != 0u) {
        throw std::runtime_error(
            "Unoccupied maximum slot contains associated hits");
      }
      if (maximum < maximaInBucket) {
        totalAssociatedHits += count;
      }
      if (totalAssociatedHits > std::numeric_limits<std::uint32_t>::max()) {
        throw std::overflow_error(
            "Total associated hit count must fit into std::uint32_t");
      }
      m_hostAssociatedHitOffsets[slot + 1u] =
          static_cast<std::uint32_t>(totalAssociatedHits);
    }
  }
  m_hostAssociatedHitIndices.resize(
      static_cast<size_type>(totalAssociatedHits));
}

template <std::size_t MaximaPerBucket>
void CudaHoughMaximumBatch<MaximaPerBucket>::allocateAssociationStorage(
    cudaStream_t stream) {
  if (!m_onDevice) {
    throw std::logic_error("CudaHoughMaximumBatch is not on the device");
  }

  if (!m_associationMetadataOnHost) {
    throw std::logic_error(
        "Association metadata must be copied before allocation");
  }

  prepareAssociationStorageHost();

  try {
    allocateDeviceColumn(m_device.associatedHitOffsets,
                         m_hostAssociatedHitOffsets.size());
    copyColumnToDevice(m_device.associatedHitOffsets,
                       m_hostAssociatedHitOffsets, stream);
    ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

    allocateDeviceColumn(m_device.associatedHitIndices,
                         m_hostAssociatedHitIndices.size());

    m_device.totalAssociatedHits =
        static_cast<std::uint32_t>(m_hostAssociatedHitIndices.size());

    m_associationStorageAllocated = true;
    m_associatedHitIndicesOnHost = false;
  } catch (...) {
    clearAssociationStorage();
    m_hostAssociatedHitOffsets.clear();
    m_hostAssociatedHitIndices.clear();
    throw;
  }
}

template <std::size_t MaximaPerBucket>
void CudaHoughMaximumBatch<MaximaPerBucket>::copyAssociatedHitIndicesToHost(
    cudaStream_t stream) {
  if (!m_onDevice) {
    throw std::logic_error("CudaHoughMaximumBatch is not on the device");
  }

  if (!m_associationStorageAllocated) {
    throw std::logic_error("Association storage has not been allocated");
  }

  copyColumnToHost(m_hostAssociatedHitIndices, m_device.associatedHitIndices,
                   stream);
  ACTS_CUDA_CHECK(cudaStreamSynchronize(stream));

  m_associatedHitIndicesOnHost = true;
}

template <std::size_t MaximaPerBucket>
void CudaHoughMaximumBatch<
    MaximaPerBucket>::clearAssociationStorage() noexcept {
  freeDeviceColumn(m_device.associatedHitOffsets);
  freeDeviceColumn(m_device.associatedHitIndices);

  m_device.totalAssociatedHits = 0u;
  m_associationStorageAllocated = false;
}

template <std::size_t MaximaPerBucket>
void CudaHoughMaximumBatch<MaximaPerBucket>::clearDevice() noexcept {
  freeDeviceColumn(m_device.tanBeta);
  freeDeviceColumn(m_device.interceptY);

  freeDeviceColumn(m_device.nHits);
  freeDeviceColumn(m_device.nLayers);
  freeDeviceColumn(m_device.layerMask);

  freeDeviceColumn(m_device.xBin);
  freeDeviceColumn(m_device.yBin);

  freeDeviceColumn(m_device.nMaxima);
  freeDeviceColumn(m_device.nAssociatedHits);

  clearAssociationStorage();

  m_device = {};
  m_onDevice = false;
}

template <std::size_t MaximaPerBucket>
void CudaHoughMaximumBatch<MaximaPerBucket>::checkBucket(
    size_type bucket) const {
  if (bucket >= nBuckets()) {
    throw std::out_of_range("CudaHoughMaximumBatch bucket index out of range");
  }
}

template <std::size_t MaximaPerBucket>
void CudaHoughMaximumBatch<MaximaPerBucket>::checkMaximum(
    size_type bucket, size_type maximum) const {
  if (maximum >= nMaxima(bucket)) {
    throw std::out_of_range("CudaHoughMaximumBatch maximum index out of range");
  }
}

}  // namespace ActsExamples
