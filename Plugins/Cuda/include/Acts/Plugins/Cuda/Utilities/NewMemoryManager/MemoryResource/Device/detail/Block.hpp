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

// CUDA include(s).
#include <cuda_runtime_api.h>

// System include(s).
#include <algorithm>
#include <limits>
#include <memory>
#include <mutex>
#include <set>
#include <unordered_map>

namespace Acts {
namespace Cuda {
namespace Nmm {
namespace MemoryRecource {
namespace detail {
namespace Arena {

// Minimum size of a superblock (256 KiB)
constexpr std::size_t minimumSuperblockSize = 1u << 18u;

class Block {

	public:
		// construct a default block.
		block() = default;

		// construct a block given a pointer and size.
		// 
		// @param[in] pointer the address for the beginning of the block.
		// @param[in] size the size of the block
		block(char* pointer, std::size_t size) : pointer_(pointer), size_(size) {}

		// construct a block given a pointer and size.
		// 
		// @param[in] pointer the address for the beginning of the block.
		// @param[in] size the size of the block
		block(void* pointer, std::size_t size) : pointer_(static_cast<char*>(pointer)), size_(size){}

		// returns the underlying pointer
		void* pointer() const { return pointer_; }

		// returns the size of the block
		std::size_t size() const { return size_; }

		// returns true if this block is valid (non-null), false otherwise
		bool isValid() const { return pointer_ != nullptr; }

		// returns true if this block is a superblock, false otherwise
		bool isSuperblock() const { return size_ >= minimumSuperblockSize; }

		// verifies wheter this block can be merged to the beginning of block b
		//
		// @param[in] b the block to check for contiguity
		// @return true if this block's `pointer` + `size` == `b.ptr`, and `not b.isHead`,
		// false otherwise 
		bool isContiguousBefore(block const& b) const { return pointer_ + size_ == b.pointer_; }

		// is this block large enough to fit that size of bytes?
		//
		// @param[in] sizeOfBytes the size in bytes to check for fit
		// @return true if this block is at least sizeOfBytes
		bool fits(std::size_t sizeOfBytes) const { return size_ >= size; }

		// split this block into two by the given size
		//
		// @param[in] size the size in bytes of the first block
		// @return std::pair<block, block> a pair of blocks split by size
		std::pair<block, block> split(std::size_t size) const {
			//assert condition of size_ >= size
			if(size_ > size) {
				return {{pointer_, size}, {pointer_ + size, size_ - size}};
			} else {
				return {*this, {}};
			}
		}

		// coalesce two contiguos blocks into one, this->isContiguousBefore(b) 
		// must be true
		// 
		// @param[in] b block to merge
		// @return block the merged block
		block merge(block const& b) const {
			//assert condition isContiguousBefore(b)
			return {pointer_, size_ + b.size};
		}

		// used by std::set to compare blocks
		bool operator<(block const& b) const { return pointer_ < b.pointer_; }

	private:
		char* pointer_;      // raw memory pointer
		std::size_t size_{}; // size in bytes
}; // class Block

} // namaspace Arena
} // namaspace detail
} // namaspace MemoryRecource
} // namaspace Nmm
} // namaspace Cuda
} // namaspace Acts