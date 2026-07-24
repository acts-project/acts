// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/detail/ContainerIterator.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <cuda_runtime.h>

namespace ActsExamples {

/// @brief Decode the detector layer from the raw MuonId representation, taken from MuonSpacePoint.
///
/// This reproduces:
///
///   m_layer = ((rawRep >> 17) & fourBit);
///
/// Bit layout relevant here:
///
///   bits 17..20 = encoded layer - 1
///
/// Therefore:
///
///   stored value 0 -> layer 0
///   stored value 1 -> layer 1
///   ...
///   stored value 15 -> layer 15
/// This function is used and we do not store layers to avoid repetition of
/// memory storage
///
/// @param rawId -> Muon ID as given to MuonSpacePoint
/// @return Layer number in range of 0..31
__host__ __device__ inline unsigned detLayer(std::uint32_t rawId) noexcept {
  static constexpr std::uint32_t fourBit = 0xFu;
  static constexpr std::uint32_t layerShift = 17u;

  return static_cast<unsigned>((rawId >> layerShift) & fourBit);
}

/// @brief Device-side raw structure-of-arrays view.
///
/// This structure does not own memory. It only stores raw device pointers.
/// CUDA kernels should receive this object by value.
struct CudaMuonSpacePointArrays {
  Acts::GeometryIdentifier::Value* geometryId = nullptr;
  std::uint32_t* muonId = nullptr;

  double* localPositionX = nullptr;
  double* localPositionY = nullptr;
  double* localPositionZ = nullptr;

  double* sensorDirectionX = nullptr;
  double* sensorDirectionY = nullptr;
  double* sensorDirectionZ = nullptr;

  double* toNextSensorX = nullptr;
  double* toNextSensorY = nullptr;
  double* toNextSensorZ = nullptr;

  double* planeNormalX = nullptr;
  double* planeNormalY = nullptr;
  double* planeNormalZ = nullptr;

  double* covariance0 = nullptr;
  double* covariance1 = nullptr;
  double* covariance2 = nullptr;

  double* driftRadius = nullptr;
  double* time = nullptr;

  std::uint32_t* bucketStart = nullptr;
  std::uint32_t* bucketEnd = nullptr;
};

/// @brief This is the RAM copy of the data. The container copies this data to VRAM
/// with moveToDevice() and copies it back with moveToHost().
struct CudaMuonSpacePointHostData {
  std::vector<Acts::GeometryIdentifier::Value> geometryId;
  std::vector<std::uint32_t> muonId;

  std::vector<double> localPositionX;
  std::vector<double> localPositionY;
  std::vector<double> localPositionZ;

  std::vector<double> sensorDirectionX;
  std::vector<double> sensorDirectionY;
  std::vector<double> sensorDirectionZ;

  std::vector<double> toNextSensorX;
  std::vector<double> toNextSensorY;
  std::vector<double> toNextSensorZ;

  std::vector<double> planeNormalX;
  std::vector<double> planeNormalY;
  std::vector<double> planeNormalZ;

  std::vector<double> covariance0;
  std::vector<double> covariance1;
  std::vector<double> covariance2;

  std::vector<double> driftRadius;
  std::vector<double> time;

  std::vector<std::uint32_t> bucketStart;
  std::vector<std::uint32_t> bucketEnd;
};

class CudaMuonSpacePointContainer;

/// @brief Host-side proxy representing one logical muon space point.
/// As specified by CompositeSpacePoint concept.
///
/// The proxy does not own data. It stores a pointer to the container and an
/// index into the flat SoA columns.
class CudaMuonSpacePointProxy {
 public:
  /// Empty default constructor.
  CudaMuonSpacePointProxy() = default;

  /// Constructor from container and space point index.
  /// @param container The container owning the data.
  /// @param index The index of the represented space point.
  CudaMuonSpacePointProxy(const CudaMuonSpacePointContainer& container,
                          std::size_t index) noexcept;

  /// Returns the muon identifier.
  const MuonSpacePoint::MuonId& id() const;

  /// Returns the geometry identifier.
  const Acts::GeometryIdentifier& geometryId() const;

  /// Returns the local measurement position.
  const Acts::Vector3& localPosition() const;

  /// Returns the local sensor direction.
  const Acts::Vector3& sensorDirection() const;

  /// Returns the normal vector to the measurement plane.
  const Acts::Vector3& planeNormal() const;

  /// Returns the vector pointing to the next wire / strip.
  const Acts::Vector3& toNextSensor() const;

  /// Returns the space point covariance values.
  const std::array<double, 3>& covariance() const;

  /// Returns the drift radius.
  double driftRadius() const;

  /// Returns the measurement time.
  double time() const;

  /// Returns whether the measurement is a straw measurement.
  bool isStraw() const;

  /// Returns whether the measurement provides time information.
  bool hasTime() const;

  /// Returns whether the measurement constrains local coordinate 0.
  bool measuresLoc0() const;

  /// Returns whether the measurement constrains local coordinate 1.
  bool measuresLoc1() const;

 private:
  const CudaMuonSpacePointContainer* m_container = nullptr;
  std::size_t m_index = 0;

  // Cache is needed because the data is returned as reference
  mutable MuonSpacePoint::MuonId m_idCache{};
  mutable Acts::GeometryIdentifier m_geometryIdCache{};
  mutable Acts::Vector3 m_localPositionCache{Acts::Vector3::Zero()};
  mutable Acts::Vector3 m_sensorDirectionCache{Acts::Vector3::Zero()};
  mutable Acts::Vector3 m_toNextSensorCache{Acts::Vector3::Zero()};
  mutable Acts::Vector3 m_planeNormalCache{Acts::Vector3::Zero()};
  mutable std::array<double, 3> m_covarianceCache{0.0, 0.0, 0.0};
};

/// @brief Pointer-like wrapper required by CompositeSpacePointContainer.
///
/// The ACTS CompositeSpacePointContainer concept expects value_type to be
/// pointer-like. This class provides operator-> and operator* over the proxy.
class CudaMuonSpacePointPtr {
 public:
  /// Type alias for the pointed-to type.
  using element_type = CudaMuonSpacePointProxy;

  /// Empty default constructor.
  CudaMuonSpacePointPtr() = default;

  /// Constructor from container and space point index.
  /// @param container The container owning the data.
  /// @param index The index of the represented space point.
  CudaMuonSpacePointPtr(const CudaMuonSpacePointContainer& container,
                        std::size_t index) noexcept;

  /// Accesses the proxy.
  element_type* operator->() noexcept { return &m_proxy; }

  /// Dereferences the proxy.
  element_type& operator*() noexcept { return m_proxy; }

  /// Checks whether this pointer-like object is valid.
  explicit operator bool() noexcept { return m_valid; }

 private:
  element_type m_proxy{};
  bool m_valid = false;
};

/// @brief CUDA-backed flat SoA container for muon space points.
///
/// The old MuonSpacePointContainer is a jagged vector of buckets. This class
/// keeps all space points in one flat SoA and stores bucket ranges separately.
class CudaMuonSpacePointContainer {
 public:
  /// Type alias for pointer-like value type.
  using value_type = CudaMuonSpacePointPtr;
  /// Type alias for container size type.
  using size_type = std::size_t;

  /// Iterator type over the flat space point container.
  template <bool read_only>
  using Iterator = Acts::detail::ContainerIterator<
      CudaMuonSpacePointContainer, CudaMuonSpacePointPtr, size_type, read_only>;

  /// Mutable iterator.
  using iterator = Iterator<false>;
  /// Const iterator.
  using const_iterator = Iterator<true>;

  /// Empty default constructor.
  CudaMuonSpacePointContainer() = default;

  /// Construct from MuonSpacePointContaier, implemented to avoid reader
  explicit CudaMuonSpacePointContainer(
      const MuonSpacePointContainer& spacePoints);

  /// Constructor with fixed number of space points.
  /// @param size The number of flat space points.
  explicit CudaMuonSpacePointContainer(size_type size);

  /// Deleted because the container owns CUDA device memory.
  CudaMuonSpacePointContainer(const CudaMuonSpacePointContainer&) = delete;
  CudaMuonSpacePointContainer& operator=(const CudaMuonSpacePointContainer&) =
      delete;

  /// Movable to transfer ownership of host/device storage.
  CudaMuonSpacePointContainer(CudaMuonSpacePointContainer&& other) noexcept;
  CudaMuonSpacePointContainer& operator=(
      CudaMuonSpacePointContainer&& other) noexcept;

  /// @brief Destructor releases VRAM pointer.
  ~CudaMuonSpacePointContainer() noexcept;

  /// @brief Returns the number of flat space points.
  size_type size() const noexcept { return m_size; }

  /// @brief Checks whether the container is empty.
  bool empty() const noexcept { return size() == 0; }

  /// @brief Returns the number of buckets.
  size_type bucketCount() const noexcept { return m_host.bucketStart.size(); }

  /// @brief Returns the first space point index of a bucket.
  /// @param bucket The bucket index.
  size_type bucketStart(size_type bucket) const;

  /// @brief Returns the past-the-end space point index of a bucket.
  /// @param bucket The bucket index.
  size_type bucketEnd(size_type bucket) const;

  /// @brief Adds a bucket range.
  /// @param start First space point index in the bucket.
  /// @param end Past-the-end space point index in the bucket.
  void addBucket(size_type start, size_type end);

  /// @brief Sets the geometry identifier.
  /// @param index The space point index.
  /// @param geometryId The raw geometry identifier.
  void setGeometryId(size_type index,
                     Acts::GeometryIdentifier::Value geometryId);

  /// @brief Sets the muon identifier.
  /// @param index The space point index.
  /// @param muonId The raw muon identifier.
  void setId(size_type index, std::uint32_t muonId);

  /// @brief Defines the coordinate system and computes the plane normal.
  /// @param index The space point index.
  /// @param position The local position.
  /// @param sensorDirection The sensor direction.
  /// @param toNextSensor The direction to the next sensor.
  void defineCoordinates(size_type index, const Acts::Vector3& position,
                         const Acts::Vector3& sensorDirection,
                         const Acts::Vector3& toNextSensor);

  /// @brief Sets the drift radius.
  /// @param index The space point index.
  /// @param radius The drift radius.
  void setRadius(size_type index, double radius);

  /// @brief Sets the measurement time.
  /// @param index The space point index.
  /// @param time The measurement time.
  void setTime(size_type index, double time);

  /// @brief Sets covariance values.
  /// @param index The space point index.
  /// @param cov0 First covariance component.
  /// @param cov1 Second covariance component.
  /// @param cov2 Time covariance component.
  void setCovariance(size_type index, double cov0, double cov1, double cov2);

  /// @brief Copies all host columns to device memory.
  void moveToDevice();

  /// @brief Copies all device columns back to host memory.
  void moveToHost();

  /// @brief Releases all device memory.
  void clearDevice() noexcept;

  /// @brief Checks whether device memory is currently allocated.
  bool isOnDevice() const noexcept { return m_onDevice; }

  /// @brief Returns raw device arrays for CUDA kernels.
  CudaMuonSpacePointArrays deviceArrays() const noexcept { return m_device; }

  /// @brief Returns a pointer-like proxy to a space point.
  /// @param index The space point index.
  value_type operator[](size_type index) noexcept {
    return value_type{*this, index};
  }

  /// @brief Returns a const pointer-like proxy to a space point.
  /// @param index The space point index.
  value_type operator[](size_type index) const noexcept {
    return value_type{*this, index};
  }

  value_type at(size_type index) const;

  iterator begin() noexcept { return {*this, 0}; }

  iterator end() noexcept { return {*this, size()}; }

  const_iterator begin() const noexcept { return {*this, 0}; }

  const_iterator end() const noexcept { return {*this, size()}; }

  // Returns raw muonId stored in space point
  std::uint32_t muonId(std::uint32_t idx) const noexcept {
    return m_host.muonId[idx];
  }

  // Write zero based layer into MuonId
  void setLogicalLayer(size_type index, std::uint32_t layer);

 private:
  /// Allow proxy to accesses private data
  friend class CudaMuonSpacePointProxy;

  size_type m_size = 0;
  CudaMuonSpacePointHostData m_host{};
  CudaMuonSpacePointArrays m_device{};
  bool m_onDevice = false;

  void checkIndex(size_type index) const;
  void checkBucket(size_type bucket) const;
};

/// Assert all concepts are fulfilled
static_assert(Acts::Experimental::CompositeSpacePoint<CudaMuonSpacePointProxy>);
static_assert(
    Acts::Experimental::CompositeSpacePointPtr<CudaMuonSpacePointPtr>);
static_assert(Acts::Experimental::CompositeSpacePointContainer<
              CudaMuonSpacePointContainer>);

}  // namespace ActsExamples
