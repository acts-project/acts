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
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/ArenaCleaner.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/CudaStreamView.hpp"

// CUDA include(s).
#include <cuda_runtime_api.h>

// System include(s)
#include <map>
#include <shared_mutex>
#include <thread>

namespace Acts {
namespace Cuda {
namespace Nmm {
namespace MemoryResource {

template <typename Upstream>
class ArenaMemoryResource final : public DeviceMemoryResource {
	public:
		explicit ArenaMemoryResource(Upstream* upstreamMemoryResource, std::size_t initialSize = globalArena::defaultInitialSize, std::size_t maximumSize = globalArena::defaultMaximumSize)
			: globalArena_{upstreamMemoryResource, initialSize, maximumSize}
			{}

		ArenaMemoryResource(ArenaMemoryResource const&) = delete;
		ArenaMemoryResource& operator=(ArenaMemoryResource const&) = delete;

		bool supportsStreams() const noexcept override { return true; }

		bool supportsGetMemInfo() const noexcept override { return false; }

	private:
		using globalArena = detail::Arena::GlobalArena<Upstream>;
		using arena = detail::Arena::Arena<Upstream>;
		using readLock = std::shared_lock<std::shared_timed_mutex>;
		using writeLock = std::lock_guard<std::shared_timed_mutex>;

		void* doAllocate(std::size_t bytes, CudaStreamView stream) override {
			if(bytes <= 0) return nullptr;

			bytes = detail::Arena::alignUp(bytes);
			return getArena(stream).allocate(bytes);
		}

		void doDeallocate(void* p, std::size_t bytes, CudaStreamView stream) override {
			if(p == nullptr || bytes <= 0) return;

			bytes = detail::Arena::alignUp(bytes);
			if(!getArena(stream).deallocate(p, bytes, stream)) {
				deallocateFromOtherArena(p, bytes, stream);
			}
		}

		void deallocateFromOtherArena(void* p, std::size_t bytes, CudaStreamView stream) {
			stream.synchronize();

			readLock lock(mtx_);

			if(usePerThreadArena(stream)) {
				auto const id = std::this_thread::get_id();
				for(auto& kv : threadArenas_) {
					if(kv.first != id && kv.second>deallocate(p, bytes)) return;
				}
			} else {
				for(auto& kv : streamArenas_){
					if(stream != kv.first && kv.second.deallocate(p, bytes)) return;
				}
			}
			globalArena_.deallocate({p, bytes});
		}

		arena& getArena(CudaStreamView stream) {
			if(usePerThreadArena(stream)){
				return getThreadArena();
			} else {
				return getStreamArena(stream);
			}
		}

		arena& getThreadArena() {
			auto const id = std::this_thread::get_id();
			
			readLock lockRead(mtx_);
			auto const it = threadArenas_.find(id);
			if(it != threadArenas_.end()) {
				return *it->second;
			}

			writeLock lockWrite(mtx_);
			auto a = std::make_shared<arena>(globalArena_);
			threadArenas_.emplace(id, a);
			thread_local detail::Arena::ArenaCleaner<Upstream> cleaner{a};
			return *a;
		}

		arena& getStreamArena(CudaStreamView stream) {
			if(usePerThreadArena(stream)){

				readLock lockRead(mtx_);
				auto const it = streamArenas_.find(stream.value());
				if (it != streamArenas_.end()) { return it->second; }

				writeLock lockWrite(mtx_);
				streamArenas_.emplace(stream.value(), globalArena_);
				return streamArenas_.at(stream.value());
			}
		}

		std::pair<std::size_t, std::size_t> doGetMemInfo(CudaStreamView stream) const override {
			return std::make_pair(0, 0);
		}

		static bool usePerThreadArena(CudaStreamView stream) {
			return stream.is_per_thread_default();
		}

		globalArena globalArena_;
		std::map<std::thread::id, std::shared_ptr<arena>> threadArenas_;
		std::map<cudaStream_t, arena> streamArenas_;
		mutable std::shared_timed_mutex mtx_;
};// class ArenaMemoryResource

} // namaspace MemoryResource
} // namaspace Nmm
} // namaspace Cuda
} // namaspace Acts