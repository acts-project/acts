// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/BinnedSurfaceMaterialAccumulator.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Material/detail/MaterialSurfaceRegistry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinAdjustment.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"
#include "Acts/Utilities/detail/WorkerPool.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

struct Acts::detail::MaterialSlotCache::Impl {
  /// The accumulated materials, in map order
  std::vector<AccumulatedSurfaceMaterial*> slots;
  /// The position of each geometry ID in `slots`
  std::unordered_map<GeometryIdentifier, std::uint32_t> index;
};

Acts::detail::MaterialSlotCache::MaterialSlotCache() = default;

Acts::detail::MaterialSlotCache::MaterialSlotCache(
    const MaterialSlotCache& /*other*/) {}

Acts::detail::MaterialSlotCache::MaterialSlotCache(
    MaterialSlotCache&& /*other*/) noexcept {}

Acts::detail::MaterialSlotCache& Acts::detail::MaterialSlotCache::operator=(
    const MaterialSlotCache& /*other*/) {
  m_impl.reset();
  return *this;
}

Acts::detail::MaterialSlotCache& Acts::detail::MaterialSlotCache::operator=(
    MaterialSlotCache&& /*other*/) noexcept {
  m_impl.reset();
  return *this;
}

Acts::detail::MaterialSlotCache::~MaterialSlotCache() = default;

Acts::detail::MaterialSlotCache::Impl& Acts::detail::MaterialSlotCache::impl() {
  if (m_impl == nullptr) {
    m_impl = std::make_unique<Impl>();
  }
  return *m_impl;
}

Acts::BinnedSurfaceMaterialAccumulator::BinnedSurfaceMaterialAccumulator(
    const Config& cfg, std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {}

std::unique_ptr<Acts::ISurfaceMaterialAccumulator::State>
Acts::BinnedSurfaceMaterialAccumulator::createState(
    const GeometryContext& gctx) const {
  auto state = std::make_unique<State>();
  state->materialSurfaceRegistry.emplace(m_cfg.materialSurfaces);

  /// Create the surface accumulation
  for (const auto& [registeredId, surface] :
       state->materialSurfaceRegistry->surfaces()) {
    GeometryIdentifier geoID = surface->geometryId();

    // Get the Surface Material
    const ISurfaceMaterial* surfaceMaterial = surface->surfaceMaterial();
    if (surfaceMaterial == nullptr) {
      throw std::invalid_argument(
          "Surface material is not set, inconsistent configuration.");
    }

    // First attempt from ProtoSurfaceMaterial
    auto psm = dynamic_cast<const ProtoSurfaceMaterial*>(surfaceMaterial);
    if (psm != nullptr) {
      auto binUtility = psm->binning();
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - (proto) binning from ProtoSurfaceMateria is "
                 << binUtility);
      // Now adjust to surface type
      binUtility = adjustBinUtility(binUtility, *surface, gctx);
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - adjusted binning is " << binUtility);
      state->accumulatedMaterial[geoID] =
          AccumulatedSurfaceMaterial(binUtility);
      // Material accumulation  is created for this
      continue;
    }
    // Second attempt from ProtoGridSurfaceMaterial
    auto psgm = dynamic_cast<const ProtoGridSurfaceMaterial*>(surfaceMaterial);
    if (psgm != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - (proto) binning from ProtoGridSurfaceMaterial is "
                 << psgm->binning());
      // Resolve the deferred binning against the surface bounds
      BinUtility binUtility(surface->localToGlobalTransform(gctx));
      binUtility += BinUtility(*resolveMultiAxis(psgm->binning(), *surface));
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - resolved binning is " << binUtility);
      state->accumulatedMaterial[geoID] =
          AccumulatedSurfaceMaterial(binUtility);
      // Material accumulation  is created for this
      continue;
    }
    // Third attempt: binned material
    auto bmp = dynamic_cast<const BinnedSurfaceMaterial*>(surfaceMaterial);
    if (bmp != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - binning from BinnedSurfaceMaterial is "
                 << bmp->binUtility());
      state->accumulatedMaterial[geoID] =
          AccumulatedSurfaceMaterial(bmp->binUtility());
      // Material accumulation  is created for this
      continue;
    }
    // Create a homogeneous type of material
    ACTS_DEBUG("       - this is homogeneous material.");
    state->accumulatedMaterial[geoID] = AccumulatedSurfaceMaterial();
  }
  return state;
}

void Acts::BinnedSurfaceMaterialAccumulator::accumulate(
    ISurfaceMaterialAccumulator::State& state, const GeometryContext& gctx,
    const std::vector<MaterialInteraction>& interactions,
    const std::vector<IAssignmentFinder::SurfaceAssignment>&
        surfacesWithoutAssignment) const {
  TrackAccumulation track{&interactions, &surfacesWithoutAssignment};
  accumulateTracks(state, gctx, std::span(&track, 1), nullptr);
}

void Acts::BinnedSurfaceMaterialAccumulator::accumulateTracks(
    ISurfaceMaterialAccumulator::State& state, const GeometryContext& /*gctx*/,
    std::span<const TrackAccumulation> tracks, detail::WorkerPool* pool) const {
  // Cast into the right state object (guaranteed by upstream algorithm)
  State& cState = static_cast<State&>(state);

  // Every accumulated surface material is an independent slot, numbered in
  // map order. The slots are cached in the state, and rebuilt when the
  // surfaces changed. A copied state starts out with an empty cache.
  auto& cache = cState.slotCache.impl();
  if (cache.slots.size() != cState.accumulatedMaterial.size()) {
    cache.slots.clear();
    cache.index.clear();
    cache.slots.reserve(cState.accumulatedMaterial.size());
    cache.index.reserve(cState.accumulatedMaterial.size());
    for (auto& [geoID, accMaterial] : cState.accumulatedMaterial) {
      cache.index.emplace(geoID,
                          static_cast<std::uint32_t>(cache.slots.size()));
      cache.slots.push_back(&accMaterial);
    }
  }
  const auto& slots = cache.slots;
  auto lookup = [&](const Surface& surface) {
    auto it = cache.index.find(surface.geometryId());
    if (it == cache.index.end()) {
      throw std::invalid_argument(
          "Surface material is not found, inconsistent configuration.");
    }
    return it->second;
  };

  // Sharding only pays off with workers and enough operations: every chunk
  // must carry at least as many operations as there are slots, so that the
  // per-chunk slot counts below never cost more than the operations. A
  // single chunk still accumulates the slots concurrently.
  std::size_t nOperations = 0;
  for (const auto& track : tracks) {
    nOperations += track.interactions->size();
    if (m_cfg.emptyBinCorrection) {
      nOperations += track.surfacesWithoutAssignment->size();
    }
  }
  constexpr std::size_t minOperationsPerChunk = 4096;
  const std::size_t nSlots = slots.size();
  const std::size_t nThreads = pool != nullptr ? pool->size() : 1;
  const std::size_t nChunks =
      std::min({tracks.size(), 4 * nThreads,
                nOperations / std::max(minOperationsPerChunk, nSlots)});

  if (nThreads <= 1 || nChunks == 0) {
    // Direct accumulation, one track after the other
    using TouchedBin = std::pair<std::uint32_t, std::array<std::size_t, 3>>;
    std::vector<TouchedBin> touched;
    std::vector<std::array<std::size_t, 3>> bins(1);
    for (const auto& track : tracks) {
      touched.clear();
      for (const auto& mi : *track.interactions) {
        std::uint32_t slot = lookup(*mi.surface);
        touched.emplace_back(
            slot, slots[slot]->accumulate(mi.intersection, mi.materialSlab,
                                          mi.pathCorrection));
      }
      // Average the first touched bin of every slot
      std::ranges::stable_sort(touched, {}, &TouchedBin::first);
      for (std::size_t i = 0; i < touched.size(); ++i) {
        if (i == 0 || touched[i].first != touched[i - 1].first) {
          bins[0] = touched[i].second;
          slots[touched[i].first]->trackAverage(bins, true);
        }
      }
      // Empty bin correction
      if (m_cfg.emptyBinCorrection) {
        for (const auto& missed : *track.surfacesWithoutAssignment) {
          slots[lookup(*missed.surface)]->trackAverage(missed.position, true);
        }
      }
    }
    return;
  }

  // A self-contained accumulation step on one slot: the slot loop below
  // reads nothing but these, sequentially
  struct Operation {
    std::uint32_t track = 0;
    bool empty = false;
    Vector3 position = Vector3::Zero();
    MaterialSlab slab = MaterialSlab::Nothing();
    double pathCorrection = 1.;
  };

  // The tracks are split into contiguous chunks, in order. Each chunk
  // resolves the slot of its operations and counts them per slot.
  auto chunkBegin = [&](std::size_t c) { return c * tracks.size() / nChunks; };
  std::vector<std::vector<std::pair<std::uint32_t, Operation>>> chunkOperations(
      nChunks);
  std::vector<std::size_t> chunkCounts(nChunks * nSlots, 0);
  pool->parallelFor(
      nChunks,
      [&](std::size_t c) {
        auto& ops = chunkOperations[c];
        std::size_t* counts = &chunkCounts[c * nSlots];
        auto add = [&](std::size_t t, const Surface& surface,
                       const Vector3& position, bool empty,
                       const MaterialSlab& slab, double pathCorrection) {
          std::uint32_t slot = lookup(surface);
          ++counts[slot];
          ops.emplace_back(slot, Operation{static_cast<std::uint32_t>(t), empty,
                                           position, slab, pathCorrection});
        };
        for (std::size_t t = chunkBegin(c); t < chunkBegin(c + 1); ++t) {
          const auto& track = tracks[t];
          for (const auto& mi : *track.interactions) {
            add(t, *mi.surface, mi.intersection, false, mi.materialSlab,
                mi.pathCorrection);
          }
          if (m_cfg.emptyBinCorrection) {
            for (const auto& missed : *track.surfacesWithoutAssignment) {
              add(t, *missed.surface, missed.position, true,
                  MaterialSlab::Nothing(), 1.);
            }
          }
        }
      },
      1);

  // Stable counting sort by slot: per slot, the chunks are laid out in
  // order, so each slot sees its operations in the serial order
  std::vector<std::size_t> offsets(nSlots + 1, 0);
  for (std::size_t slot = 0; slot < nSlots; ++slot) {
    std::size_t total = offsets[slot];
    for (std::size_t c = 0; c < nChunks; ++c) {
      std::size_t count = chunkCounts[c * nSlots + slot];
      chunkCounts[c * nSlots + slot] = total;
      total += count;
    }
    offsets[slot + 1] = total;
  }
  // Keep the payloads in their chunks and sort only pointers to them. The
  // chunks are no longer modified, so these pointers remain valid until
  // accumulation has finished.
  std::vector<const Operation*> operations(offsets.back());
  pool->parallelFor(
      nChunks,
      [&](std::size_t c) {
        std::size_t* fill = &chunkCounts[c * nSlots];
        for (const auto& [slot, op] : chunkOperations[c]) {
          operations[fill[slot]++] = &op;
        }
      },
      1);

  // Process the slots concurrently, the busiest ones first
  std::vector<std::uint32_t> active;
  for (std::uint32_t slot = 0; slot < nSlots; ++slot) {
    if (offsets[slot + 1] > offsets[slot]) {
      active.push_back(slot);
    }
  }
  std::ranges::stable_sort(active, [&](auto a, auto b) {
    return offsets[a + 1] - offsets[a] > offsets[b + 1] - offsets[b];
  });

  pool->parallelFor(
      active.size(),
      [&](std::size_t i) {
        const std::uint32_t slot = active[i];
        AccumulatedSurfaceMaterial& accMaterial = *slots[slot];
        // Same as the serial accumulation: only the first touched bin of a
        // track is track-averaged, before any empty bin correction
        std::vector<std::array<std::size_t, 3>> bins;
        std::optional<std::array<std::size_t, 3>> touchedBin;
        auto averageTouched = [&]() {
          if (touchedBin.has_value()) {
            bins.assign(1, *touchedBin);
            accMaterial.trackAverage(bins, true);
            touchedBin.reset();
          }
        };
        std::uint32_t currentTrack = std::numeric_limits<std::uint32_t>::max();
        for (std::size_t o = offsets[slot]; o < offsets[slot + 1]; ++o) {
          const Operation& op = *operations[o];
          if (op.track != currentTrack) {
            averageTouched();
            currentTrack = op.track;
          }
          if (!op.empty) {
            auto bin =
                accMaterial.accumulate(op.position, op.slab, op.pathCorrection);
            if (!touchedBin.has_value()) {
              touchedBin = bin;
            }
          } else {
            averageTouched();
            accMaterial.trackAverage(op.position, true);
          }
        }
        averageTouched();
      },
      1);
}

std::map<Acts::GeometryIdentifier,
         std::shared_ptr<const Acts::ISurfaceMaterial>>
Acts::BinnedSurfaceMaterialAccumulator::finalizeMaterial(
    ISurfaceMaterialAccumulator::State& state,
    const GeometryContext& /*gctx*/) const {
  std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
      sMaterials;

  // Cast into the right state object (guaranteed by upstream algorithm)
  State* cState = static_cast<State*>(&state);
  if (cState == nullptr) {
    throw std::invalid_argument(
        "Invalid state object provided, something is seriously wrong.");
  }

  // iterate over the map to call the total average
  for (auto& accMaterial : cState->accumulatedMaterial) {
    ACTS_DEBUG("Finalizing map for Surface " << accMaterial.first);
    auto sMaterial = accMaterial.second.totalAverage();
    sMaterials[accMaterial.first] = std::move(sMaterial);
  }

  return sMaterials;
}

Acts::TrackingGeometryMaterial
Acts::BinnedSurfaceMaterialAccumulator::finalizeMaps(
    ISurfaceMaterialAccumulator::State& state,
    const GeometryContext& gctx) const {
  auto* concrete = dynamic_cast<State*>(&state);
  if (concrete == nullptr || !concrete->materialSurfaceRegistry) {
    throw std::invalid_argument(
        "State was not created by this material accumulator");
  }
  return concrete->materialSurfaceRegistry->materialMaps(
      finalizeMaterial(state, gctx));
}
