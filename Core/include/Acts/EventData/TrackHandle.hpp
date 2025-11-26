// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <any>
#include <cassert>
#include <memory>

namespace Acts {

/// @brief Type-erased handle for track proxies
///
/// This class provides type erasure for TrackProxy instances without heap
/// allocation per handle. It uses a static handler pattern with virtual
/// functions where each concrete track type has a static handler instance.
/// The handle stores only a pointer to the actual track data without ownership.
///
/// @note The handle does not own the underlying track data. The user must
/// ensure the referenced track outlives the handle.
///
/// Example usage:
/// @code
/// TrackContainer tracks{...};
/// auto track = tracks.getTrack(0);
/// TrackHandle handle(track);  // Type-erased, no heap allocation
/// auto tipIndex = handle.tipIndex();
/// @endcode
class TrackHandle {
 public:
  /// @brief Handler base class for type-specific operations
  ///
  /// Each concrete track type has a static instance of a derived handler that
  /// implements these virtual methods. This provides the vtable for type
  /// erasure.
  class TrackHandlerBase {
   public:
    virtual ~TrackHandlerBase() = default;

    /// Get the tip index of the track
    virtual TrackIndexType tipIndex(const void* instance) const = 0;

    /// Get the stem index of the track
    virtual TrackIndexType stemIndex(const void* instance) const = 0;

    /// Get the particle hypothesis
    virtual const ParticleHypothesis& particleHypothesis(
        const void* instance) const = 0;

    /// Get the reference surface
    virtual const Surface& referenceSurface(const void* instance) const = 0;

    /// Get track parameters
    virtual BoundTrackParameters parameters(const void* instance) const = 0;

    /// Get track covariance
    virtual BoundSquareMatrix covariance(const void* instance) const = 0;

    /// Check if track has reference surface
    virtual bool hasReferenceSurface(const void* instance) const = 0;

    /// Get number of track states
    virtual std::size_t nTrackStates(const void* instance) const = 0;

    /// Get number of measurements
    virtual std::size_t nMeasurements(const void* instance) const = 0;

    /// Get number of outliers
    virtual std::size_t nOutliers(const void* instance) const = 0;

    /// Get number of holes
    virtual std::size_t nHoles(const void* instance) const = 0;

    /// Get number of shared hits
    virtual std::size_t nSharedHits(const void* instance) const = 0;

    /// Get chi2
    virtual double chi2(const void* instance) const = 0;

    /// Get number of degrees of freedom
    virtual std::size_t nDoF(const void* instance) const = 0;

    /// Get component as std::any
    virtual std::any component(const void* instance,
                               HashedString key) const = 0;
  };

  /// @brief Handler implementation for a concrete track type
  ///
  /// This template generates a static handler instance for each unique
  /// track proxy type. The virtual methods static_cast the void pointer
  /// to the concrete type and forward to the actual track proxy methods.
  template <typename track_proxy_t>
  class TrackHandler final : public TrackHandlerBase {
   public:
    /// Get the static handler instance for this track type
    static const TrackHandler& instance() {
      static const TrackHandler s_instance;
      return s_instance;
    }

    TrackIndexType tipIndex(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->tipIndex();
    }

    TrackIndexType stemIndex(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->stemIndex();
    }

    const ParticleHypothesis& particleHypothesis(
        const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->particleHypothesis();
    }

    const Surface& referenceSurface(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->referenceSurface();
    }

    BoundTrackParameters parameters(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->parameters();
    }

    BoundSquareMatrix covariance(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->covariance();
    }

    bool hasReferenceSurface(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->hasReferenceSurface();
    }

    std::size_t nTrackStates(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->nTrackStates();
    }

    std::size_t nMeasurements(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->nMeasurements();
    }

    std::size_t nOutliers(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->nOutliers();
    }

    std::size_t nHoles(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->nHoles();
    }

    std::size_t nSharedHits(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->nSharedHits();
    }

    double chi2(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->chi2();
    }

    std::size_t nDoF(const void* instance) const override {
      return static_cast<const track_proxy_t*>(instance)->nDoF();
    }

    std::any component(const void* instance,
                       HashedString key) const override {
      return static_cast<const track_proxy_t*>(instance)->component(key);
    }

   private:
    TrackHandler() = default;
  };

  /// @brief Construct an empty (invalid) handle
  TrackHandle() : m_ptr(nullptr), m_handler(nullptr) {}

  /// @brief Construct from a concrete track proxy
  ///
  /// @param track Reference to a track proxy (any concrete type)
  ///
  /// @note The handle does not take ownership. The track must outlive the
  /// handle.
  template <typename track_container_t, typename trajectory_t,
            template <typename> class holder_t, bool read_only>
  TrackHandle(
      TrackProxy<track_container_t, trajectory_t, holder_t, read_only>& track)
      : m_ptr(static_cast<void*>(&track)),
        m_handler(&TrackHandler<TrackProxy<track_container_t, trajectory_t,
                                           holder_t, read_only>>::instance()) {}

  /// @brief Construct from a const track proxy
  template <typename track_container_t, typename trajectory_t,
            template <typename> class holder_t, bool read_only>
  TrackHandle(const TrackProxy<track_container_t, trajectory_t, holder_t,
                                read_only>& track)
      : m_ptr(const_cast<void*>(static_cast<const void*>(&track))),
        m_handler(&TrackHandler<TrackProxy<track_container_t, trajectory_t,
                                           holder_t, read_only>>::instance()) {}

  /// @brief Check if the handle is valid (non-null)
  bool isValid() const { return m_ptr != nullptr && m_handler != nullptr; }

  /// @brief Explicit conversion to bool
  explicit operator bool() const { return isValid(); }

  /// @brief Get the tip index of the track
  TrackIndexType tipIndex() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->tipIndex(m_ptr);
  }

  /// @brief Get the stem index of the track
  TrackIndexType stemIndex() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->stemIndex(m_ptr);
  }

  /// @brief Get the particle hypothesis
  const ParticleHypothesis& particleHypothesis() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->particleHypothesis(m_ptr);
  }

  /// @brief Get the reference surface
  const Surface& referenceSurface() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->referenceSurface(m_ptr);
  }

  /// @brief Get track parameters
  BoundTrackParameters parameters() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->parameters(m_ptr);
  }

  /// @brief Get track covariance
  BoundSquareMatrix covariance() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->covariance(m_ptr);
  }

  /// @brief Check if the track has a reference surface
  bool hasReferenceSurface() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->hasReferenceSurface(m_ptr);
  }

  /// @brief Get the number of track states
  std::size_t nTrackStates() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->nTrackStates(m_ptr);
  }

  /// @brief Get the number of measurements
  std::size_t nMeasurements() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->nMeasurements(m_ptr);
  }

  /// @brief Get the number of outliers
  std::size_t nOutliers() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->nOutliers(m_ptr);
  }

  /// @brief Get the number of holes
  std::size_t nHoles() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->nHoles(m_ptr);
  }

  /// @brief Get the number of shared hits
  std::size_t nSharedHits() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->nSharedHits(m_ptr);
  }

  /// @brief Get the chi2 of the track
  double chi2() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->chi2(m_ptr);
  }

  /// @brief Get the number of degrees of freedom
  std::size_t nDoF() const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->nDoF(m_ptr);
  }

  /// @brief Get a dynamic component by key
  std::any component(HashedString key) const {
    assert(isValid() && "TrackHandle is not valid");
    return m_handler->component(m_ptr, key);
  }

  /// @brief Get a typed dynamic component by key
  template <typename T>
  T componentAs(HashedString key) const {
    return std::any_cast<T>(component(key));
  }

 private:
  void* m_ptr;  // Pointer to the actual track (no ownership)
  const TrackHandlerBase* m_handler;  // Static handler for type-specific ops
};

}  // namespace Acts
