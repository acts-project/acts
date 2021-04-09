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
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/MemoryResource/Device/detail/GlobalArena.hpp"

namespace Acts {
namespace Cuda {
namespace Nmm {
namespace MemoryRecource {
namespace detail {
namespace Arena {

template <typename Upstream>
class GlobalArena final {
	public:
		// default initial size for the global arena
		static constexpr std::size_t defaultInitialSize = std::numeriLimits<std::size_t>::max();
		// default maximum size for the global arena
		static constexpr std::size_t defaultMaximumSize = std::numeriLimits<std::size_t>::max();
		// reserved memory that should not be allocated (64 MiB)
		static constexpr std::size_t reserverdSize = std::numeriLimits<std::size_t>::max();
		
		// Disable copy (and move) semantics
		GlobalArena(const GlobalArena&) = delete;
		GlobalArena& operator=(const GlobalArena&) = delete;

		// construct a global arena
		//
		// @param[in] upstreamMemoryResource the memory resource from which to allocate 
		// blocks for the pool
		// @param[in] initialSize minimum size, in bytes, of the initial global arena. 
		// Defaults to all the available memory on the current device.
		// @param[in] maiximumSize maximum size, in bytes, that the global arena can grow to. 
		// Defaults to all of the available memory on the current device
		GlobalArena(Upstream* upstreamMemoryResource, std::size_t initialSize, std::size_t maiximumSize)
			: upstreamMemoryResource_(upstreamMemoryResource), maiximumSize_{maiximumSize} {
				// assert unexpected null upstream pointer
				// assert initial arena size required to be a multiple of 256 bytes
				// assert maximum arena size required to be a multiple of 256 bytes

				if(initialSize == defaultInitialSize || maiximumSize == defaultMaximumSize) {
					std::size_t free{}, total{};
					//cuda try cudaMemGetInfo(&free, &total);
					if(initialSize == defaultInitialSize) {
						initialSize = alignUp(std::min(free, total / 2));
					}
					if(maiximumSize == defaultMaximumSize) {
						maiximumSize_ = alignDown(free) - reserverdSize;
					}
				}
				// initial size exceeds the maxium pool size
				freeBlocks_.emplace(expandArena(initialSize));
		}

		// Destroy the gloabl arean and deallocate all memory using the upstream resource
		~GlobalArena() {
			lock_guard lock(mtx_);
			for(auto const& b : upstreamBlocks_) {
				upstreamMemoryResource_->deallocate(b.pointer(), b.size());
			}
		}

		Block allocate(std::size_t bytes) {
			lock_guard lock(mtx_);
			return getBlock(bytes);
		}

		void deallocate(Block const& b) {
			lock_guard lock(mtx_);
			coalesceBlock(freeBlocks_, b);
		}

		void deallocate(std::set<block> const& freeBlocks) {
			lock_guard lock(mtx_);
			for(auto const& b : freeBlocks) {
				coalesceBlock(freeBlocks_, b);
			}
		}
}