// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/CudaMuonSpacePoint.hpp"
#include "ActsExamples/Utilities/CudaUtilities.hpp"

#include <stdexcept>
#include <utility>

#include <cuda_runtime.h>

namespace {

using namespace ActsExamples;

template <typename T>
void resizeColumn(std::vector<T>& column, std::size_t size) {
  column.resize(size);
}

void resizeHostData(ActsExamples::CudaMuonSpacePointHostData& host,
                    std::size_t size) {
  resizeColumn(host.geometryId, size);
  resizeColumn(host.muonId, size);

  resizeColumn(host.localPositionX, size);
  resizeColumn(host.localPositionY, size);
  resizeColumn(host.localPositionZ, size);

  resizeColumn(host.sensorDirectionX, size);
  resizeColumn(host.sensorDirectionY, size);
  resizeColumn(host.sensorDirectionZ, size);

  resizeColumn(host.toNextSensorX, size);
  resizeColumn(host.toNextSensorY, size);
  resizeColumn(host.toNextSensorZ, size);

  resizeColumn(host.planeNormalX, size);
  resizeColumn(host.planeNormalY, size);
  resizeColumn(host.planeNormalZ, size);

  resizeColumn(host.covariance0, size);
  resizeColumn(host.covariance1, size);
  resizeColumn(host.covariance2, size);

  resizeColumn(host.driftRadius, size);
  resizeColumn(host.time, size);
}

void allocateDeviceData(ActsExamples::CudaMuonSpacePointArrays& device,
                        std::size_t spacePointCount, std::size_t bucketCount) {
  allocateDeviceColumn(device.geometryId, spacePointCount);
  allocateDeviceColumn(device.muonId, spacePointCount);

  allocateDeviceColumn(device.localPositionX, spacePointCount);
  allocateDeviceColumn(device.localPositionY, spacePointCount);
  allocateDeviceColumn(device.localPositionZ, spacePointCount);

  allocateDeviceColumn(device.sensorDirectionX, spacePointCount);
  allocateDeviceColumn(device.sensorDirectionY, spacePointCount);
  allocateDeviceColumn(device.sensorDirectionZ, spacePointCount);

  allocateDeviceColumn(device.toNextSensorX, spacePointCount);
  allocateDeviceColumn(device.toNextSensorY, spacePointCount);
  allocateDeviceColumn(device.toNextSensorZ, spacePointCount);

  allocateDeviceColumn(device.planeNormalX, spacePointCount);
  allocateDeviceColumn(device.planeNormalY, spacePointCount);
  allocateDeviceColumn(device.planeNormalZ, spacePointCount);

  allocateDeviceColumn(device.covariance0, spacePointCount);
  allocateDeviceColumn(device.covariance1, spacePointCount);
  allocateDeviceColumn(device.covariance2, spacePointCount);

  allocateDeviceColumn(device.driftRadius, spacePointCount);
  allocateDeviceColumn(device.time, spacePointCount);

  allocateDeviceColumn(device.bucketStart, bucketCount);
  allocateDeviceColumn(device.bucketEnd, bucketCount);
}

void freeDeviceData(ActsExamples::CudaMuonSpacePointArrays& device) noexcept {
  freeDeviceColumn(device.geometryId);
  freeDeviceColumn(device.muonId);

  freeDeviceColumn(device.localPositionX);
  freeDeviceColumn(device.localPositionY);
  freeDeviceColumn(device.localPositionZ);

  freeDeviceColumn(device.sensorDirectionX);
  freeDeviceColumn(device.sensorDirectionY);
  freeDeviceColumn(device.sensorDirectionZ);

  freeDeviceColumn(device.toNextSensorX);
  freeDeviceColumn(device.toNextSensorY);
  freeDeviceColumn(device.toNextSensorZ);

  freeDeviceColumn(device.planeNormalX);
  freeDeviceColumn(device.planeNormalY);
  freeDeviceColumn(device.planeNormalZ);

  freeDeviceColumn(device.covariance0);
  freeDeviceColumn(device.covariance1);
  freeDeviceColumn(device.covariance2);

  freeDeviceColumn(device.driftRadius);
  freeDeviceColumn(device.time);

  freeDeviceColumn(device.bucketStart);
  freeDeviceColumn(device.bucketEnd);
}

void copyHostToDevice(ActsExamples::CudaMuonSpacePointArrays& device,
                      const ActsExamples::CudaMuonSpacePointHostData& host) {
  copyColumnToDevice(device.geometryId, host.geometryId);
  copyColumnToDevice(device.muonId, host.muonId);

  copyColumnToDevice(device.localPositionX, host.localPositionX);
  copyColumnToDevice(device.localPositionY, host.localPositionY);
  copyColumnToDevice(device.localPositionZ, host.localPositionZ);

  copyColumnToDevice(device.sensorDirectionX, host.sensorDirectionX);
  copyColumnToDevice(device.sensorDirectionY, host.sensorDirectionY);
  copyColumnToDevice(device.sensorDirectionZ, host.sensorDirectionZ);

  copyColumnToDevice(device.toNextSensorX, host.toNextSensorX);
  copyColumnToDevice(device.toNextSensorY, host.toNextSensorY);
  copyColumnToDevice(device.toNextSensorZ, host.toNextSensorZ);

  copyColumnToDevice(device.planeNormalX, host.planeNormalX);
  copyColumnToDevice(device.planeNormalY, host.planeNormalY);
  copyColumnToDevice(device.planeNormalZ, host.planeNormalZ);

  copyColumnToDevice(device.covariance0, host.covariance0);
  copyColumnToDevice(device.covariance1, host.covariance1);
  copyColumnToDevice(device.covariance2, host.covariance2);

  copyColumnToDevice(device.driftRadius, host.driftRadius);
  copyColumnToDevice(device.time, host.time);

  copyColumnToDevice(device.bucketStart, host.bucketStart);
  copyColumnToDevice(device.bucketEnd, host.bucketEnd);
}

void copyDeviceToHost(ActsExamples::CudaMuonSpacePointHostData& host,
                      const ActsExamples::CudaMuonSpacePointArrays& device) {
  copyColumnToHost(host.geometryId, device.geometryId);
  copyColumnToHost(host.muonId, device.muonId);

  copyColumnToHost(host.localPositionX, device.localPositionX);
  copyColumnToHost(host.localPositionY, device.localPositionY);
  copyColumnToHost(host.localPositionZ, device.localPositionZ);

  copyColumnToHost(host.sensorDirectionX, device.sensorDirectionX);
  copyColumnToHost(host.sensorDirectionY, device.sensorDirectionY);
  copyColumnToHost(host.sensorDirectionZ, device.sensorDirectionZ);

  copyColumnToHost(host.toNextSensorX, device.toNextSensorX);
  copyColumnToHost(host.toNextSensorY, device.toNextSensorY);
  copyColumnToHost(host.toNextSensorZ, device.toNextSensorZ);

  copyColumnToHost(host.planeNormalX, device.planeNormalX);
  copyColumnToHost(host.planeNormalY, device.planeNormalY);
  copyColumnToHost(host.planeNormalZ, device.planeNormalZ);

  copyColumnToHost(host.covariance0, device.covariance0);
  copyColumnToHost(host.covariance1, device.covariance1);
  copyColumnToHost(host.covariance2, device.covariance2);

  copyColumnToHost(host.driftRadius, device.driftRadius);
  copyColumnToHost(host.time, device.time);

  copyColumnToHost(host.bucketStart, device.bucketStart);
  copyColumnToHost(host.bucketEnd, device.bucketEnd);
}

std::size_t countSpacePoints(
    const ActsExamples::MuonSpacePointContainer& spacePoints) {
  std::size_t totalSize = 0;

  for (const auto& bucket : spacePoints) {
    totalSize += bucket.size();
  }

  return totalSize;
}

}  // namespace

namespace ActsExamples {

CudaMuonSpacePointProxy::CudaMuonSpacePointProxy(
    const CudaMuonSpacePointContainer& container, std::size_t index) noexcept
    : m_container{&container}, m_index{index} {}

const MuonSpacePoint::MuonId& CudaMuonSpacePointProxy::id() const {
  m_idCache = MuonSpacePoint::MuonId{m_container->m_host.muonId[m_index]};
  return m_idCache;
}

const Acts::GeometryIdentifier& CudaMuonSpacePointProxy::geometryId() const {
  m_geometryIdCache =
      Acts::GeometryIdentifier{m_container->m_host.geometryId[m_index]};
  return m_geometryIdCache;
}

const Acts::Vector3& CudaMuonSpacePointProxy::localPosition() const {
  m_localPositionCache = Acts::Vector3{
      m_container->m_host.localPositionX[m_index],
      m_container->m_host.localPositionY[m_index],
      m_container->m_host.localPositionZ[m_index],
  };

  return m_localPositionCache;
}

const Acts::Vector3& CudaMuonSpacePointProxy::sensorDirection() const {
  m_sensorDirectionCache = Acts::Vector3{
      m_container->m_host.sensorDirectionX[m_index],
      m_container->m_host.sensorDirectionY[m_index],
      m_container->m_host.sensorDirectionZ[m_index],
  };

  return m_sensorDirectionCache;
}

const Acts::Vector3& CudaMuonSpacePointProxy::planeNormal() const {
  m_planeNormalCache = Acts::Vector3{
      m_container->m_host.planeNormalX[m_index],
      m_container->m_host.planeNormalY[m_index],
      m_container->m_host.planeNormalZ[m_index],
  };

  return m_planeNormalCache;
}

const Acts::Vector3& CudaMuonSpacePointProxy::toNextSensor() const {
  m_toNextSensorCache = Acts::Vector3{
      m_container->m_host.toNextSensorX[m_index],
      m_container->m_host.toNextSensorY[m_index],
      m_container->m_host.toNextSensorZ[m_index],
  };

  return m_toNextSensorCache;
}

const std::array<double, 3>& CudaMuonSpacePointProxy::covariance() const {
  m_covarianceCache = {
      m_container->m_host.covariance0[m_index],
      m_container->m_host.covariance1[m_index],
      m_container->m_host.covariance2[m_index],
  };

  return m_covarianceCache;
}

double CudaMuonSpacePointProxy::driftRadius() const {
  return m_container->m_host.driftRadius[m_index];
}

double CudaMuonSpacePointProxy::time() const {
  return m_container->m_host.time[m_index];
}

bool CudaMuonSpacePointProxy::isStraw() const {
  return id().technology() == MuonSpacePoint::MuonId::TechField::Mdt;
}

bool CudaMuonSpacePointProxy::hasTime() const {
  return id().measuresTime();
}

bool CudaMuonSpacePointProxy::measuresLoc0() const {
  return id().measuresPhi();
}

bool CudaMuonSpacePointProxy::measuresLoc1() const {
  return id().measuresEta();
}

CudaMuonSpacePointPtr::CudaMuonSpacePointPtr(
    const CudaMuonSpacePointContainer& container, std::size_t index) noexcept
    : m_proxy{container, index}, m_valid{true} {}

CudaMuonSpacePointContainer::CudaMuonSpacePointContainer(size_type size)
    : m_size{size} {
  resizeHostData(m_host, m_size);
}

CudaMuonSpacePointContainer::CudaMuonSpacePointContainer(
    const MuonSpacePointContainer& spacePoints)
    : CudaMuonSpacePointContainer{countSpacePoints(spacePoints)} {
  m_host.bucketStart.reserve(spacePoints.size());
  m_host.bucketEnd.reserve(spacePoints.size());

  size_type index = 0;

  for (const MuonSpacePointBucket& bucket : spacePoints) {
    const size_type bucketBegin = index;

    for (const MuonSpacePoint& spacePoint : bucket) {
      setGeometryId(index, spacePoint.geometryId().value());
      setId(index, spacePoint.id().toInt());

      defineCoordinates(index, spacePoint.localPosition(),
                        spacePoint.sensorDirection(),
                        spacePoint.toNextSensor());

      setRadius(index, spacePoint.driftRadius());
      setTime(index, spacePoint.time());

      const auto& covariance = spacePoint.covariance();
      setCovariance(index, covariance[0], covariance[1], covariance[2]);

      ++index;
    }

    // Empty buckets are preserved as ranges for which begin == end.
    addBucket(bucketBegin, index);
  }
}

CudaMuonSpacePointContainer::CudaMuonSpacePointContainer(
    CudaMuonSpacePointContainer&& other) noexcept
    : m_size{std::exchange(other.m_size, 0)},
      m_host{std::move(other.m_host)},
      m_device{std::exchange(other.m_device, {})},
      m_onDevice{std::exchange(other.m_onDevice, false)} {}

CudaMuonSpacePointContainer& CudaMuonSpacePointContainer::operator=(
    CudaMuonSpacePointContainer&& other) noexcept {
  if (this != &other) {
    clearDevice();

    m_size = std::exchange(other.m_size, 0);
    m_host = std::move(other.m_host);
    m_device = std::exchange(other.m_device, {});
    m_onDevice = std::exchange(other.m_onDevice, false);
  }

  return *this;
}

CudaMuonSpacePointContainer::~CudaMuonSpacePointContainer() noexcept {
  clearDevice();
}

CudaMuonSpacePointContainer::size_type CudaMuonSpacePointContainer::bucketStart(
    size_type bucket) const {
  checkBucket(bucket);
  return m_host.bucketStart[bucket];
}

CudaMuonSpacePointContainer::size_type CudaMuonSpacePointContainer::bucketEnd(
    size_type bucket) const {
  checkBucket(bucket);
  return m_host.bucketEnd[bucket];
}

void CudaMuonSpacePointContainer::addBucket(size_type start, size_type end) {
  if (start > end || end > m_size) {
    throw std::out_of_range(
        std::format("CudaMuonSpacePointContainer invalid bucket range "
                    "[{:};{:}]. Allowed [0;{:})",
                    start, end, m_size));
  }

  m_host.bucketStart.push_back(static_cast<std::uint32_t>(start));
  m_host.bucketEnd.push_back(static_cast<std::uint32_t>(end));
}

void CudaMuonSpacePointContainer::setGeometryId(
    size_type index, Acts::GeometryIdentifier::Value geometryId) {
  checkIndex(index);
  m_host.geometryId[index] = geometryId;
}

void CudaMuonSpacePointContainer::setId(size_type index, std::uint32_t muonId) {
  checkIndex(index);
  m_host.muonId[index] = muonId;
}

void CudaMuonSpacePointContainer::defineCoordinates(
    size_type index, const Acts::Vector3& position,
    const Acts::Vector3& sensorDirection, const Acts::Vector3& toNextSensor) {
  checkIndex(index);

  m_host.localPositionX[index] = position.x();
  m_host.localPositionY[index] = position.y();
  m_host.localPositionZ[index] = position.z();

  m_host.sensorDirectionX[index] = sensorDirection.x();
  m_host.sensorDirectionY[index] = sensorDirection.y();
  m_host.sensorDirectionZ[index] = sensorDirection.z();

  m_host.toNextSensorX[index] = toNextSensor.x();
  m_host.toNextSensorY[index] = toNextSensor.y();
  m_host.toNextSensorZ[index] = toNextSensor.z();

  // CompositeSpacePoint expects the measurement plane normal
  const Acts::Vector3 planeNormal =
      sensorDirection.cross(toNextSensor).normalized();

  m_host.planeNormalX[index] = planeNormal.x();
  m_host.planeNormalY[index] = planeNormal.y();
  m_host.planeNormalZ[index] = planeNormal.z();
}

void CudaMuonSpacePointContainer::setRadius(size_type index, double radius) {
  checkIndex(index);
  m_host.driftRadius[index] = radius;
}

void CudaMuonSpacePointContainer::setTime(size_type index, double time) {
  checkIndex(index);
  m_host.time[index] = time;
}

void CudaMuonSpacePointContainer::setCovariance(size_type index, double cov0,
                                                double cov1, double cov2) {
  checkIndex(index);

  m_host.covariance0[index] = cov0;
  m_host.covariance1[index] = cov1;
  m_host.covariance2[index] = cov2;
}

void CudaMuonSpacePointContainer::moveToDevice() {
  clearDevice();

  allocateDeviceData(m_device, m_size, bucketCount());
  copyHostToDevice(m_device, m_host);

  m_onDevice = true;
}

void CudaMuonSpacePointContainer::moveToHost() {
  if (!m_onDevice) {
    return;
  }

  copyDeviceToHost(m_host, m_device);
}

void CudaMuonSpacePointContainer::clearDevice() noexcept {
  freeDeviceData(m_device);
  m_device = {};
  m_onDevice = false;
}

CudaMuonSpacePointContainer::value_type CudaMuonSpacePointContainer::at(
    size_type index) const {
  checkIndex(index);
  return value_type{*this, index};
}

inline void CudaMuonSpacePointContainer::checkIndex(size_type index) const {
  if (index >= m_size) {
    throw std::out_of_range(std::format(
        "CudaMuonSpacePointContainer index {:}out of range. Max allowed: {:}",
        index, m_size - 1));
  }
}

void CudaMuonSpacePointContainer::checkBucket(size_type bucket) const {
  if (bucket >= bucketCount()) {
    throw std::out_of_range(std::format(
        "CudaMuonSpacePointContainer bucket {:} out of range. Max allowed: {:}",
        bucket, bucketCount() - 1));
  }
}

// The raw MuonId stores layer - 1 in bits 17..20. Here `layer` is already
// zero-based because CUDA layer masks use layer 0 -> bit 0.
void CudaMuonSpacePointContainer::setLogicalLayer(size_type index,
                                                  std::uint32_t layer) {
  checkIndex(index);

  static constexpr std::uint32_t fourBit = 0xFu;
  static constexpr std::uint32_t layerShift = 17u;

  if (layer > fourBit) {
    throw std::out_of_range(
        "CudaMuonSpacePointContainer logical layer out of range");
  }

  std::uint32_t rawId = m_host.muonId[index];

  // Clear old layer bits 17..20.
  rawId &= ~(fourBit << layerShift);

  // Insert new layer bits 17..20.
  rawId |= ((layer & fourBit) << layerShift);

  m_host.muonId[index] = rawId;
}

}  // namespace ActsExamples
