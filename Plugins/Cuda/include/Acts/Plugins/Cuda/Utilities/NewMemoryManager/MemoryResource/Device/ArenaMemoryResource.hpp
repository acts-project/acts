// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/DeciceMemoryResource.hpp"
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/CudaStreamView.hpp"

// CUDA include(s).
#include <cuda_runtime_api.h>

// System include(s)
#include <map>
#include <shared_mutex>

namespace Acts {
namespace Cuda {
namespace Nmm {
namespace MemoryRecource {

template <typename Upstream>
class ArenaMemoryResource final : public DeviceMemoryResource {
	public:
		explicit ArenaMemoryResource(Upstrean* upstreamMemoryResource, std::size_t initialSize = GlobalArena::defaultInitialSize, std::size_t maximumSize = GlobalArena::defaultMaximumSize)
			: globalArena_{upstreamMemoryResource,initialSize, maximumSize}
			{}

		ArenaMemoryResource(ArenaMemoryResource const&) = delete;
		ArenaMemoryResource& operator=(ArenaMemoryResource const&) delete;

		bool supportsStreams() const noexcept override { return true; }

		bool supportsGetMemInfo() const noexcept override { return false; }

	private:
		using globalArena = detail::Arena::GlobalArena<Upstrean>;
		using arena = detail::Arena::Arena<Upstrean>;
		using readLock = std::shared_lock<std::shared_time_mutex>;
		using writeLock = std::lock_guard<std::shared_time_mutex>;

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

		getArena(CudaStreamView stream) {
			if(usePerThreadArena(stream)){
				return getThreadArena();
			} else {
				return getStreamArena(stream);
			}
		}

		Arena& getThreadArena() {
			auto const id = std::this_thread::get_id();
			
			readLock lock(mtx_);
			auto const it = threadArenas_.find(id);
			if(it != threadArenas_.end()) {
				return *it->second;
			}

			writeLock lock(mtx_);
			auto a = std::make_shared<Arena>(globalArena_);
			threadArenas_.emplace(id, a);
			thread_local detail::Arena::ArenaCleaner<Upstrean> cleaner{a};
			return *a;
		}

		Arena& getStreamArena(CudaStreamView stream) {
			if(usePerThreadArena(stream)){

				read_lock lock(mtx_);
				auto const it = stream_arenas_.find(stream.value());
				if (it != stream_arenas_.end()) { return it->second; }

				write_lock lock(mtx_);
				stream_arenas_.emplace(stream.value(), global_arena_);
				return stream_arenas_.at(stream.value());
			}
		}

		std::pair<std::size_t, std::size_t> doGetMemInfo(CudaStreamView stream) const override {
			return std::make_pair(0, 0);
		}

		static bool usePerThreadArena(CudaStreamView stream) {
			return stream.isPerThreadDefault();
		}

		GlobalArena globalArena_;
		std::map<std::thread::id, std::shared_ptr<Arena>> threadArenas_;
		std::map<cudaStream_t, arena> streamArenas_;
		mutable std::shared_timed_mutex mtx_;
};// class ArenaMemoryResource

} // namaspace MemoryRecource
} // namaspace Nmm
} // namaspace Cuda
} // namaspace Acts