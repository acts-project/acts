// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/detail/Aligned.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/CudaStreamView.hpp"

// System include(s)
#include <cstddef>
#include <utility>

namespace Acts {
namespace Cuda {
namespace Nmm {
namespace MemoryResource {

class DeviceMemoryResource {
	public:
		virtual ~DeviceMemoryResource() = default;

		void* allocate(std::size_t bytes, CudaStreamView stream = CudaStreamView{}) {
			return doAllocate(Nmm::detail::align_up(bytes, 8), stream);
		}

		void deallocate(void* p, std::size_t bytes, CudaStreamView stream = CudaStreamView{}) {
			doDeallocate(p, Nmm::detail::align_up(bytes, 8), stream);
		}

		bool isEqual(DeviceMemoryResource const& other) const noexcept {
			return doIsEqual(other);
		}

		virtual bool supportsStreams() const noexcept = 0;

		virtual bool supportsGetMemInfo() const noexcept = 0;

		std::pair<std::size_t, std::size_t> getMemInfo(CudaStreamView stream) const {
			return doGetMemInfo(stream);
		}

	private:
		virtual void* doAllocate(std::size_t bytes, CudaStreamView stream) = 0;

		virtual void* doDeallocate(void* p, std::size_t bytes, CudaStreamView stream) = 0;

		virtual bool doIsEqual(DeviceMemoryResource const& other) const noexcept {
			return this == &other;
		}

		virtual std::pair<std::size_t, std::size_t> doGetMemInfo(CudaStreamView stream) const = 0;

};// class DeviceMemoryResource

} // namaspace MemoryResource
} // namaspace Nmm
} // namaspace Cuda
} // namaspace Acts