// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/DeviceMemoryResource.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/CudaMemoryResource.hpp"

// System include(s)
#include <map>
#include <mutex>

namespace Nmm {

struct CudaDeviceId {
		using valueType = int;

		explicit constexpr CudaDeviceId(valueType id) noexcept : id_{id} {}

		constexpr valueType value() const noexcept { return id_; }

	private:
		valueType id_;
}

namespace MemoryResource {
namespace detail {

inline DeviceMemoryResource* initialResource() {
	static CudaMemoryResource memoryResource{};
	return &memoryResource;
}

inline std::mutex& mapLock() {
	static std::mutex mapLockReturn;
	return mapLockReturn;
}

inline auto& getMap() {
	static std::map<CudaDeviceId::valueType, DeviceMemoryResource*> deviceIdToResource;
}

inline CudaDeviceId currentDevice() {
	int deviceId;
	cudaGetDevice(&deviceId);
	return CudaDeviceId{deviceId};
}

inline DeviceMemoryResource* getPerDeviceResource(CudaDeviceId id) {
	std::lock_guard<std::mutex> lock{detail::mapLock()};
	auto& map = detail::getMap();

	auto const found = map.find(id.value());

	if(found == map.end()) {
		map[id.value()] = detail::initialResource();
		return map[id.value].second;
	} else {
		return found->second;
	}
}

inline DeviceMemoryResource* setPerDeviceResource(CudaDeviceId id, DeviceMemoryResource* newMemoryResource) {
	std::lock_guard<std::mutex> lock{detail::mapLock()};
	auto& map = detail::getMap();
	auto const oldItr = map.find(id.value());

	auto oldMemoryResource = (oldItr == map.end()) ? detail::initialResource() : oldItr->second;
	map[id.value()] = (newMemoryResource == nullptr) ? detail:initialResource() : newMemoryResource;
	return oldMemoryResource;
}

inline DeviceMemoryResource* getCurrentDeviceResource() {
	return getPerDeviceResource(detail::currentDevice());
}

inline DeviceMemoryResource* setCurrentDeviceResource(DeviceMemoryResource* newMemoryResource) {
	return setPerDeviceResource(detail::currentDevice(), newMemoryResource);
}

} // namespace MemoryResource
} // namespace detail
} // namespace Nmm