// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/CudaHoughTransformUtils.hpp"

#ifdef ACTS_ENABLE_CUDA

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cuda_runtime.h>

namespace {

using ActsExamples::CudaMuonSpacePointArrays;
using ActsExamples::CudaHoughTransformUtils::CudaHoughPlaneBatchArrays;
using ActsExamples::CudaHoughTransformUtils::HoughAxisRanges;
using ActsExamples::CudaHoughTransformUtils::LayerMask;
using ActsExamples::CudaHoughTransformUtils::YieldType;

void cudaCheck(cudaError_t status) {
  if (status != cudaSuccess) {
    throw std::runtime_error(cudaGetErrorString(status));
  }
}

template <typename T>
void allocateDeviceColumn(T*& deviceColumn, std::size_t size) {
  if (size == 0) {
    return;
  }

  cudaCheck(
      cudaMalloc(reinterpret_cast<void**>(&deviceColumn), size * sizeof(T)));
}

template <typename T>
void freeDeviceColumn(T*& deviceColumn) noexcept {
  if (deviceColumn != nullptr) {
    cudaFree(deviceColumn);
    deviceColumn = nullptr;
  }
}

template <typename T>
void copyColumnToDevice(T* deviceColumn, const std::vector<T>& hostColumn) {
  if (hostColumn.empty()) {
    return;
  }

  cudaCheck(cudaMemcpy(deviceColumn, hostColumn.data(),
                       hostColumn.size() * sizeof(T), cudaMemcpyHostToDevice));
}

template <typename T>
void copyColumnToHost(std::vector<T>& hostColumn, const T* deviceColumn) {
  if (hostColumn.empty() || deviceColumn == nullptr) {
    return;
  }

  cudaCheck(cudaMemcpy(hostColumn.data(), deviceColumn,
                       hostColumn.size() * sizeof(T), cudaMemcpyDeviceToHost));
}

/// Decode the zero-based layer index from raw MuonId.
__host__ __device__ inline unsigned layerIndexFromMuonId(std::uint32_t rawId) {
  static constexpr std::uint32_t fourBit = 0xFu;
  static constexpr std::uint32_t layerShift = 17u;

  return static_cast<unsigned>((rawId >> layerShift) & fourBit);
}

LayerMask layerBitHost(unsigned layer) {
#ifdef ACTS_ENABLE_CUDA_RUNTIME_CHECKS
  if (layer >= 8u * sizeof(LayerMask)) {
    throw std::out_of_range(
        "CudaHoughPlaneBatch supports at most 64 logical layers");
  }
#endif

  return LayerMask{1ull} << layer;
}

__device__ LayerMask layerBitDevice(unsigned layer) {

#ifdef ACTS_ENABLE_CUDA_RUNTIME_CHECKS
  if (layer >= 8u * sizeof(LayerMask)) {
    return LayerMask{0ull};
  }
#endif

  return LayerMask{1ull} << layer;
}

__device__ double binCenterDevice(double min, double max, unsigned nSteps,
                                  unsigned binIndex) {
  return min + (max - min) * 0.5 * (2.0 * binIndex + 1.0) / nSteps;
}

__device__ int binIndexDevice(double min, double max, unsigned nSteps,
                              double val) {
  return static_cast<int>((val - min) / (max - min) * nSteps);
}

__device__ void fillSharedBin(YieldType* sHits, YieldType* sLayers,
                              LayerMask* sMask, std::uint32_t nBinsX,
                              std::uint32_t xBin, std::uint32_t yBin,
                              unsigned layer, YieldType weight) {
  const std::uint32_t localBin = yBin * nBinsX + xBin;

  atomicAdd(&sHits[localBin], weight);

  // For example for 2 get 0..0010
  const LayerMask bit = layerBitDevice(layer);

#ifdef ACTS_ENABLE_CUDA_RUNTIME_CHECKS
  // Out of range check
  // Can be removed later
  if (bit == LayerMask{0ull}) {
    return;
  }
#endif

  // atomic returns old value
  const LayerMask oldMask = atomicOr(&sMask[localBin], bit);

  
  // if layer bit was 0 & 1 then we know this if first time and we add layer number
  if ((oldMask & bit) == LayerMask{0ull}) {
    atomicAdd(&sLayers[localBin], weight);
  }
}

__device__ void fillSharedYBand(YieldType* sHits, YieldType* sLayers,
                                LayerMask* sMask,
                                const CudaHoughPlaneBatchArrays batch,
                                const HoughAxisRanges ranges,
                                std::uint32_t xBin, double yCenter,
                                double yHalfWidth, unsigned layer,
                                YieldType weight) {
  int yBinDown = binIndexDevice(ranges.yMin, ranges.yMax, batch.nBinsY,
                                yCenter - yHalfWidth);
  int yBinUp = binIndexDevice(ranges.yMin, ranges.yMax, batch.nBinsY,
                              yCenter + yHalfWidth);

#ifdef ACTS_ENABLE_CUDA_RUNTIME_CHECKS
  if (yBinDown > yBinUp) {
    const int tmp = yBinDown;
    yBinDown = yBinUp;
    yBinUp = tmp;
  }
#endif

  // Those are necessary
  if (yBinDown < 0) {
    yBinDown = 0;
  }

  if (yBinUp >= static_cast<int>(batch.nBinsY)) {
    yBinUp = static_cast<int>(batch.nBinsY) - 1;
  }

  // Top hat add, so 1 for all values
  for (int yBin = yBinDown; yBin <= yBinUp; ++yBin) {
    fillSharedBin(sHits, sLayers, sMask, batch.nBinsX, xBin,
                  static_cast<std::uint32_t>(yBin), layer, weight);
  }
}

__global__ void fillEtaDriftCirclesMdtBatchKernel(
    CudaHoughPlaneBatchArrays batch, CudaMuonSpacePointArrays spacePoints,
    HoughAxisRanges ranges, double widthScale, double maxWidth,
    YieldType weight) {
  const std::uint32_t nCells = batch.nBinsX * batch.nBinsY;

  extern __shared__ unsigned char sharedMemory[];

  auto* sHits = reinterpret_cast<YieldType*>(sharedMemory);
  auto* sLayers = sHits + nCells;

  std::size_t maskOffsetBytes = 2u * nCells * sizeof(YieldType);
  maskOffsetBytes =
      ((maskOffsetBytes + sizeof(LayerMask) - 1u) / sizeof(LayerMask)) *
      sizeof(LayerMask);

  auto* sMask = reinterpret_cast<LayerMask*>(sharedMemory + maskOffsetBytes);

  // Grid-stride loop over buckets.
  for (std::uint32_t bucket = blockIdx.x; bucket < batch.nBuckets; bucket += gridDim.x) {
    // Clear this block's shared-memory accumulator for the current bucket.
    for (std::uint32_t i = threadIdx.x; i < nCells; i += blockDim.x) {
      sHits[i] = 0.0f;
      sLayers[i] = 0.0f;
      sMask[i] = LayerMask{0ull};
    }

    __syncthreads();

    const std::uint32_t bucketStart = spacePoints.bucketStart[bucket];
    const std::uint32_t bucketEnd = spacePoints.bucketEnd[bucket];
    const std::uint32_t nHits = bucketEnd - bucketStart;

    constexpr std::uint32_t nSolutions = 2u;
    const std::uint32_t nTasks = nHits * batch.nBinsX * nSolutions;

    // Thread-stride loop inside one bucket.
    // One logical task is hit x tanTheta-bin x left/right drift-circle solution
    for (std::uint32_t task = threadIdx.x; task < nTasks; task += blockDim.x) {
      const std::uint32_t solution = task % nSolutions;
      const std::uint32_t xBin = (task / nSolutions) % batch.nBinsX;
      const std::uint32_t localHit = task / (nSolutions * batch.nBinsX);
      const std::uint32_t hitIndex = bucketStart + localHit;

      const double tanTheta =
          binCenterDevice(ranges.xMin, ranges.xMax, batch.nBinsX, xBin);

      const double y = spacePoints.localPositionY[hitIndex];
      const double z = spacePoints.localPositionZ[hitIndex];
      const double r = spacePoints.driftRadius[hitIndex];

      const double sign = solution == 0u ? -1.0 : 1.0;
      const double y0 =
          y - tanTheta * z + sign * r * sqrt(1.0 + tanTheta * tanTheta);

      const double cov = spacePoints.covariance1[hitIndex] > 0.0
                             ? spacePoints.covariance1[hitIndex]
                             : 0.0;

      double width = sqrt(cov) * widthScale;

      if (width > maxWidth) {
        width = maxWidth;
      }

      const unsigned layer = layerIndexFromMuonId(spacePoints.muonId[hitIndex]);

      fillSharedYBand(sHits, sLayers, sMask, batch, ranges, xBin, y0, width,
                      layer, weight);
    }

    __syncthreads();

    // Write this bucket's accumulator to global memory.
    const std::uint32_t globalBase = bucket * nCells;

    for (std::uint32_t i = threadIdx.x; i < nCells; i += blockDim.x) {
      batch.nHits[globalBase + i] = sHits[i];
      batch.nLayers[globalBase + i] = sLayers[i];
      batch.layerMask[globalBase + i] = sMask[i];
    }

    __syncthreads();
  }
}

void allocateDeviceData(CudaHoughPlaneBatchArrays& device,
                        std::size_t totalCells) {
  allocateDeviceColumn(device.nHits, totalCells);
  allocateDeviceColumn(device.nLayers, totalCells);
  allocateDeviceColumn(device.layerMask, totalCells);
}

void freeDeviceData(CudaHoughPlaneBatchArrays& device) noexcept {
  freeDeviceColumn(device.nHits);
  freeDeviceColumn(device.nLayers);
  freeDeviceColumn(device.layerMask);
}

void copyHostToDevice(CudaHoughPlaneBatchArrays& device,
                      const std::vector<YieldType>& hits,
                      const std::vector<YieldType>& layers,
                      const std::vector<LayerMask>& layerMask) {
  copyColumnToDevice(device.nHits, hits);
  copyColumnToDevice(device.nLayers, layers);
  copyColumnToDevice(device.layerMask, layerMask);
}

void copyDeviceToHost(std::vector<YieldType>& hits,
                      std::vector<YieldType>& layers,
                      std::vector<LayerMask>& layerMask,
                      const CudaHoughPlaneBatchArrays& device) {
  copyColumnToHost(hits, device.nHits);
  copyColumnToHost(layers, device.nLayers);
  copyColumnToHost(layerMask, device.layerMask);
}

std::size_t sharedBytesForCells(std::size_t nCells) {
  std::size_t bytes = 2u * nCells * sizeof(YieldType);
  std::size_t layerMaskSize = alignof(LayerMask); // for normal types alignof and sizeof should regturn same
  bytes = ((bytes + layerMaskSize - 1u) / layerMaskSize) * layerMaskSize;
  bytes += nCells * sizeof(LayerMask);
  return bytes;
}

}  // namespace

namespace ActsExamples::CudaHoughTransformUtils {

CudaHoughPlaneBatch::CudaHoughPlaneBatch(const HoughPlaneConfig& cfg,
                                         size_type nBuckets)
    : m_cfg{cfg}, m_nBuckets{nBuckets} {
  if (m_cfg.nBinsX == 0 || m_cfg.nBinsY == 0) {
    throw std::invalid_argument(
        "CudaHoughPlaneBatch requires non-zero nBinsX and nBinsY");
  }

  if (m_nBuckets == 0) {
    throw std::invalid_argument(
        "CudaHoughPlaneBatch requires non-zero nBuckets");
  }

  if (m_cfg.nBinsX > std::numeric_limits<std::uint32_t>::max() ||
      m_cfg.nBinsY > std::numeric_limits<std::uint32_t>::max() ||
      m_nBuckets > std::numeric_limits<std::uint32_t>::max()) {
    throw std::overflow_error(
        "CudaHoughPlaneBatch dimensions must fit into std::uint32_t");
  }

  m_hostHits.resize(totalCells(), 0.0f);
  m_hostLayers.resize(totalCells(), 0.0f);
  m_hostLayerMask.resize(totalCells(), LayerMask{0ull});
}

CudaHoughPlaneBatch::CudaHoughPlaneBatch(CudaHoughPlaneBatch&& other) noexcept
    : m_cfg{other.m_cfg},
      m_nBuckets{other.m_nBuckets},
      m_hostHits{std::move(other.m_hostHits)},
      m_hostLayers{std::move(other.m_hostLayers)},
      m_hostLayerMask{std::move(other.m_hostLayerMask)},
      m_device{std::exchange(other.m_device, {})},
      m_onDevice{std::exchange(other.m_onDevice, false)} {
  other.m_cfg = {};
  other.m_nBuckets = 0;
}

CudaHoughPlaneBatch& CudaHoughPlaneBatch::operator=(
    CudaHoughPlaneBatch&& other) noexcept {
  if (this != &other) {
    clearDevice();

    m_cfg = other.m_cfg;
    m_nBuckets = other.m_nBuckets;
    m_hostHits = std::move(other.m_hostHits);
    m_hostLayers = std::move(other.m_hostLayers);
    m_hostLayerMask = std::move(other.m_hostLayerMask);
    m_device = std::exchange(other.m_device, {});
    m_onDevice = std::exchange(other.m_onDevice, false);

    other.m_cfg = {};
    other.m_nBuckets = 0;
  }

  return *this;
}

CudaHoughPlaneBatch::~CudaHoughPlaneBatch() noexcept {
  clearDevice();
}

CudaHoughPlaneBatch::size_type CudaHoughPlaneBatch::globalBin(
    size_type bucket, size_type xBin, size_type yBin) const {
  checkBucket(bucket);
  checkIndices(xBin, yBin);
  return uncheckedGlobalBin(bucket, xBin, yBin);
}

std::pair<std::size_t, std::size_t> CudaHoughPlaneBatch::axisBins(
    size_type globalBin) const {

  const size_type localBin = globalBin % nCellsPerBucket();
  return {localBin % nBinsX(), localBin / nBinsX()};
}

void CudaHoughPlaneBatch::fillBin(size_type bucket, size_type xBin,
                                  size_type yBin, unsigned layer,
                                  YieldType weight) {
  checkBucket(bucket);
  checkIndices(xBin, yBin);

  if (weight == 0.0f) {
    return;
  }

  const size_type bin = uncheckedGlobalBin(bucket, xBin, yBin);

  m_hostHits[bin] += weight;

  const LayerMask bit = layerBitHost(layer);

  if ((m_hostLayerMask[bin] & bit) == LayerMask{0ull}) {
    m_hostLayerMask[bin] |= bit;
    m_hostLayers[bin] += weight;
  }
}

void CudaHoughPlaneBatch::fillEtaDriftCirclesHost(
    const CudaMuonSpacePointContainer& spacePoints,
    const HoughAxisRanges& axisRanges,
    YieldType weight) {

  // No width, however original implementation has this as widthPar:
  // CoordType dy = widthPar(x, measurement);
  // {0, 0} means that for loop will default to single central element
  double widthScale = 0; 
  double maxWidth = 0;

  for (size_type bucket = 0; bucket < nBuckets(); ++bucket) {
    const size_type start = spacePoints.bucketStart(bucket);
    const size_type end = spacePoints.bucketEnd(bucket);

    for (size_type hitIndex = start; hitIndex < end; ++hitIndex) {
      auto sp = spacePoints[hitIndex];

      const double y = sp->localPosition().y();
      const double z = sp->localPosition().z();
      const double r = sp->driftRadius();

      const double cov = std::max(sp->covariance()[1], 0.0);
      const double width = std::min(std::sqrt(cov) * widthScale, maxWidth);

      const unsigned layer = layerIndexFromMuonId(spacePoints.muonId(hitIndex));

      for (size_type xBin = 0; xBin < nBinsX(); ++xBin) {
        const double tanTheta = Acts::HoughTransformUtils::binCenter(
            axisRanges.xMin, axisRanges.xMax, nBinsX(), xBin);

        const double sqrtTerm = std::sqrt(1.0 + tanTheta * tanTheta);

        const double yLeft = y - tanTheta * z - r * sqrtTerm;
        const double yRight = y - tanTheta * z + r * sqrtTerm;

        for (double y0 : {yLeft, yRight}) {
          int yBinDown = Acts::HoughTransformUtils::binIndex(
              axisRanges.yMin, axisRanges.yMax, nBinsY(), y0 - width);
          int yBinUp = Acts::HoughTransformUtils::binIndex(
              axisRanges.yMin, axisRanges.yMax, nBinsY(), y0 + width);

          if (yBinDown > yBinUp) {
            std::swap(yBinDown, yBinUp);
          }

          yBinDown = std::max(yBinDown, 0);
          yBinUp = std::min(yBinUp, static_cast<int>(nBinsY()) - 1);

          for (int yBin = yBinDown; yBin <= yBinUp; ++yBin) {
            fillBin(bucket, xBin, static_cast<size_type>(yBin), layer, weight);
          }
        }
      }
    }
  }
}

void CudaHoughPlaneBatch::fillEtaDriftCirclesOnDevice(
    CudaMuonSpacePointContainer& spacePoints, const HoughAxisRanges& axisRanges,
    YieldType weight,
    std::uint32_t threadsPerBlock, std::uint32_t num_blocks) {

  // No width, however original implementation has this as widthPar:
  // CoordType dy = widthPar(x, measurement);
  // {0, 0} means that for loop will default to single central element
  double widthScale = 0; 
  double maxWidth = 0;

  if (threadsPerBlock == 0) {
    throw std::invalid_argument("threadsPerBlock must be non-zero");
  }

  if (!spacePoints.isOnDevice()) {
    spacePoints.moveToDevice();
  }

  if (!m_onDevice) {
    moveToDevice();
  }

  const std::size_t sharedBytes = sharedBytesForCells(nCellsPerBucket());

  int device = 0;
  cudaCheck(cudaGetDevice(&device));

  int smCount = 0;
  cudaCheck(
      cudaDeviceGetAttribute(&smCount, cudaDevAttrMultiProcessorCount, device));

  int maxSharedMemory = 0;
  cudaCheck(cudaDeviceGetAttribute(&maxSharedMemory,
                                   cudaDevAttrMaxSharedMemoryPerBlock, device));

  int maxThreadsPerBlock = 0;
  cudaCheck(cudaDeviceGetAttribute(&maxThreadsPerBlock,
                                   cudaDevAttrMaxThreadsPerBlock, device));

  if (sharedBytes > static_cast<std::size_t>(maxSharedMemory)) {
    throw std::runtime_error(
        "CudaHoughPlaneBatch MDT shared-memory fill requires too much shared "
        "memory. Use fewer Hough bins or implement a global-memory fallback.");
  }

  if (threadsPerBlock > static_cast<std::uint32_t>(maxThreadsPerBlock)) {
    throw std::runtime_error(
        "threadsPerBlock exceeds cudaDevAttrMaxThreadsPerBlock");
  }

  // If num_blocks == 0, use one block per SM.
  if (num_blocks == 0) {
    num_blocks = static_cast<std::uint32_t>(smCount);
  }

  if (num_blocks == 0) {
    throw std::runtime_error("Resolved num_blocks is zero");
  }

  num_blocks = std::min<std::uint32_t>(num_blocks,
                                       static_cast<std::uint32_t>(nBuckets()));

  fillEtaDriftCirclesMdtBatchKernel<<<static_cast<unsigned>(num_blocks),
                                      static_cast<unsigned>(threadsPerBlock),
                                      sharedBytes>>>(
      m_device, spacePoints.deviceArrays(), axisRanges, widthScale, maxWidth,
      weight);

  cudaCheck(cudaGetLastError());
  cudaCheck(cudaDeviceSynchronize());
}

YieldType CudaHoughPlaneBatch::nHits(size_type bucket, size_type xBin,
                                     size_type yBin) const {
  return m_hostHits[globalBin(bucket, xBin, yBin)];
}

YieldType CudaHoughPlaneBatch::nLayers(size_type bucket, size_type xBin,
                                       size_type yBin) const {
  return m_hostLayers[globalBin(bucket, xBin, yBin)];
}

LayerMask CudaHoughPlaneBatch::layerMask(size_type bucket, size_type xBin,
                                         size_type yBin) const {
  return m_hostLayerMask[globalBin(bucket, xBin, yBin)];
}

bool CudaHoughPlaneBatch::hasLayer(size_type bucket, size_type xBin,
                                   size_type yBin, unsigned layer) const {
  if (layer >= 8u * sizeof(LayerMask)) {
    return false;
  }

  const LayerMask bit = LayerMask{1ull} << layer;
  return (layerMask(bucket, xBin, yBin) & bit) != LayerMask{0ull};
}

std::vector<unsigned> CudaHoughPlaneBatch::layers(size_type bucket,
                                                  size_type xBin,
                                                  size_type yBin) const {
  const LayerMask mask = layerMask(bucket, xBin, yBin);

  std::vector<unsigned> out{};

  for (unsigned layer = 0; layer < 8u * sizeof(LayerMask); ++layer) {
    const LayerMask bit = LayerMask{1ull} << layer;

    if ((mask & bit) != LayerMask{0ull}) {
      out.push_back(layer);
    }
  }

  return out;
}

YieldType CudaHoughPlaneBatch::maxHits(size_type bucket) const {
  checkBucket(bucket);

  const auto begin = m_hostHits.begin() +
                     static_cast<std::ptrdiff_t>(bucket * nCellsPerBucket());
  const auto end = begin + static_cast<std::ptrdiff_t>(nCellsPerBucket());

  return *std::max_element(begin, end);
}

YieldType CudaHoughPlaneBatch::maxLayers(size_type bucket) const {
  checkBucket(bucket);

  const auto begin = m_hostLayers.begin() +
                     static_cast<std::ptrdiff_t>(bucket * nCellsPerBucket());
  const auto end = begin + static_cast<std::ptrdiff_t>(nCellsPerBucket());

  return *std::max_element(begin, end);
}

 std::pair<std::size_t, std::size_t>CudaHoughPlaneBatch::locMaxHits(
    size_type bucket) const {

  checkBucket(bucket);

  const auto begin = m_hostHits.begin() +
                     static_cast<std::ptrdiff_t>(bucket * nCellsPerBucket());
  const auto end = begin + static_cast<std::ptrdiff_t>(nCellsPerBucket());
  const auto iter = std::max_element(begin, end);

  const size_type localBin = static_cast<size_type>(std::distance(begin, iter));

  return {localBin % nBinsX(), localBin / nBinsX()};
}

void CudaHoughPlaneBatch::moveToDevice() {
  clearDevice();

  m_device.nBuckets = static_cast<std::uint32_t>(nBuckets());
  m_device.nBinsX = static_cast<std::uint32_t>(nBinsX());
  m_device.nBinsY = static_cast<std::uint32_t>(nBinsY());

  allocateDeviceData(m_device, totalCells());
  copyHostToDevice(m_device, m_hostHits, m_hostLayers, m_hostLayerMask);

  m_onDevice = true;
}

void CudaHoughPlaneBatch::moveToHost() {
  if (!m_onDevice) {
    return;
  }

  copyDeviceToHost(m_hostHits, m_hostLayers, m_hostLayerMask, m_device);
}

void CudaHoughPlaneBatch::clearDevice() noexcept {
  freeDeviceData(m_device);
  m_device = {};
  m_onDevice = false;
}


void CudaHoughPlaneBatch::checkBucket(size_type bucket) const {
  if (bucket >= nBuckets()) {
    throw std::out_of_range("CudaHoughPlaneBatch bucket index out of range");
  }
}

void CudaHoughPlaneBatch::checkIndices(size_type xBin, size_type yBin) const {
  if (xBin >= nBinsX()) {
    throw std::out_of_range("CudaHoughPlaneBatch x-bin index out of range");
  }

  if (yBin >= nBinsY()) {
    throw std::out_of_range("CudaHoughPlaneBatch y-bin index out of range");
  }
}

}  // namespace ActsExamples::CudaHoughTransformUtils

#endif
