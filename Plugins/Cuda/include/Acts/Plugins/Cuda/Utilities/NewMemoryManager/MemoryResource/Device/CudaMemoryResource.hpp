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
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/CudaStreamView.hpp"

namespace Acts {
namespace Cuda {
namespace Nmm {
namespace MemoryResource {

class CudaMemoryResource final : public DeviceMemoryResource {
	public:
		CudaMemoryResource() = default;
		~CudaMemoryResource() = default;
		CudaMemoryResource(CudaMemoryResource const&) = default;
		CudaMemoryResource(CudaMemoryResource&&) = default;
		CudaMemoryResource& operator=(CudaMemoryResource const&) = default;
		CudaMemoryResource& operator=(CudaMemoryResource&&) = default;

		bool supportsStreams() const noexcept override { return false; }

		bool supportsGetMemInfo() const noexcept override { return true; }

	private:
		void* doAllocate(std::size_t bytes, CudaStreamView) override {
			void* p{nullptr};
			cudaMalloc(&p, bytes);
			return p;
		}

		void doDeallocate(void* p, std::size_t, CudaStreamView) override {
			cudaFree(p);
		}

		bool doIsEqual(DeviceMemoryResource const& other) const noexcept override {
			return dynamic_cast<CudaMemoryResource const*>(&other) != nullptr;
		}

		std::pair<size_t, size_t> doGetMemInfo(CudaStreamView) const override {
			std::size_t freeSize;
			std::size_t totalSize;

			cudaMemGetInfo(&freeSize, &totalSize);
			return std::make_pair(freeSize, totalSize);
		}

};// class CudaMemoryResource

} // namespace MemoryResource
} // namespace Nmm
} // namespace Cuda
} // namespace Acts