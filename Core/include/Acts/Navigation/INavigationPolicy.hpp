// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationDelegate.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Utilities/Any.hpp"

#include <type_traits>

namespace Acts {

class TrackingVolume;
class INavigationPolicy;
class Surface;
class Navigator;

class NavigationPolicyStateManager;

/// Small-buffer size (in bytes) for a type-erased navigation policy state.
/// Sized so the largest policy state is stored inline, keeping the
/// per-volume-entry pushState() on the navigation hot path free of heap
/// allocations. A static_assert in pushState() guards against silent heap
/// fallback: bump this value if a policy needs a larger state.
static constexpr std::size_t kNavigationPolicyStateSbo = 128;

/// Type-erased storage for a single navigation policy state, with small-buffer
/// optimization so the common (small) states avoid heap allocation.
using NavigationPolicyStateStorage = AnyBase<kNavigationPolicyStateSbo>;

/// Wrapper class for a navigation policy state stored in the state manager.
/// This class provides type-safe access to the underlying state through
/// the `as<T>()` method. The state is stored as a type-erased, small-buffer
/// optimized value (@ref NavigationPolicyStateStorage) in the
/// NavigationPolicyStateManager.
class NavigationPolicyState {
 public:
  /// Cast the state to a specific type
  /// @tparam T The type to cast to
  /// @return Reference to the state as type T
  template <typename T>
  T& as() {
    return payload().as<T>();
  }

  /// Cast the state to a specific type (const version)
  /// @tparam T The type to cast to
  /// @return Const reference to the state as type T
  template <typename T>
  const T& as() const {
    return payload().as<T>();
  }

  NavigationPolicyState() = default;

  /// Check if the state is empty (not associated with a manager)
  /// @return True if the state is empty, false otherwise
  bool empty() const { return m_manager == nullptr; }

  /// Whether this state is the default sentinel, i.e. either no state at all or
  /// the @c INavigationPolicy::EmptyState pushed by policies that do not
  /// override createState(). Such policies also use the default (always-true)
  /// isValid(), so the navigator treats a default state as "no validity
  /// constraint" and skips the per-step isValid() call. A policy with a
  /// meaningful isValid() necessarily pushes a real state to check against, so
  /// this signal cannot be forgotten.
  /// @return True if this is the default/empty state
  bool isDefault() const;

  /// Absolute index of this state within its manager's state stack.
  /// @return the stack index
  std::size_t index() const { return m_index; }

  /// Access another state managed by the same manager, by absolute stack index.
  /// Composite policies use this to reach their child states, which are pushed
  /// contiguously onto the stack directly below the composite's own state.
  /// @param index absolute stack index of the state to access
  /// @return wrapper for the state at @p index in the same manager
  NavigationPolicyState atIndex(std::size_t index) const {
    return NavigationPolicyState{*m_manager, index};
  }

 private:
  NavigationPolicyState(NavigationPolicyStateManager& manager,
                        std::size_t index)
      : m_manager(&manager), m_index(index) {}

  NavigationPolicyStateStorage& payload();
  const NavigationPolicyStateStorage& payload() const;

  NavigationPolicyStateManager* m_manager = nullptr;
  std::size_t m_index = 0;

  friend class NavigationPolicyStateManager;
};

/// Manager class for navigation policy states. Maintains a stack of
/// type-erased states and provides access to the current state.
/// Navigation policies use this manager to push and pop their states
/// during navigation.
class NavigationPolicyStateManager {
 public:
  /// Push a new state onto the stack and construct it in-place
  /// @tparam T The type of state to create
  /// @tparam Args The types of the constructor arguments
  /// @param args Arguments to forward to the state constructor
  /// @return Reference to the newly created state
  template <typename T, typename... Args>
  T& pushState(Args&&... args) {
    static_assert(
        sizeof(T) <= kNavigationPolicyStateSbo,
        "Navigation policy state does not fit the small buffer and would be "
        "heap-allocated on the navigation hot path. Increase "
        "kNavigationPolicyStateSbo to accommodate it.");
    // Reserve a typical stack depth on first use only. This keeps the stack
    // from reallocating during ordinary navigation, while not allocating at all
    // for navigators that never push a state (e.g. the Gen1 navigator, which
    // carries this manager but does not use it). Once cleared by reset() the
    // capacity is retained, so this stays a no-op afterwards.
    if (m_stateStack.empty()) {
      m_stateStack.reserve(16);
    }
    NavigationPolicyStateStorage& state = m_stateStack.emplace_back();
    return state.emplace<T>(std::forward<Args>(args)...);
  }

  friend class Navigator;

  /// Get the current (top) state from the stack
  /// @return NavigationPolicyState wrapper for the current state, or empty state if stack is empty
  NavigationPolicyState currentState() {
    if (m_stateStack.empty()) {
      return {};  // Empty state as sentinel
    }
    return NavigationPolicyState{*this, m_stateStack.size() - 1};
  }

  /// Remove the current state from the stack
  /// @throws std::runtime_error if the stack is empty
  void popState() {
    if (m_stateStack.empty()) {
      throw std::runtime_error(
          "NavigationPolicyStateManager: Attempt to pop from empty stack");
    }
    m_stateStack.pop_back();
  }

  /// Completely reset the state manager by clearing the stack
  void reset() { m_stateStack.clear(); }

 private:
  friend class NavigationPolicyState;

  // We might want to extend this to a stack
  std::vector<NavigationPolicyStateStorage> m_stateStack;
};

inline NavigationPolicyStateStorage& NavigationPolicyState::payload() {
  if (m_manager == nullptr) {
    throw std::runtime_error(
        "NavigationPolicyState: Attempt to access empty payload");
  }
  return m_manager->m_stateStack.at(m_index);
}

inline const NavigationPolicyStateStorage& NavigationPolicyState::payload()
    const {
  if (m_manager == nullptr) {
    throw std::runtime_error(
        "NavigationPolicyState: Attempt to access empty payload");
  }
  return m_manager->m_stateStack.at(m_index);
}

namespace detail {
template <typename T>
concept HasOldInitializeCandidates = requires {
  requires requires(T policy, const GeometryContext& gctx,
                    const NavigationArguments& args,
                    AppendOnlyNavigationStream& stream, const Logger& logger) {
    policy.initializeCandidates(gctx, args, stream, logger);
  };
};

template <typename T>
concept HasNewInitializeCandidates = requires {
  requires requires(T policy, const GeometryContext& gctx,
                    const NavigationArguments& args,
                    NavigationPolicyState& state,
                    AppendOnlyNavigationStream& stream, const Logger& logger) {
    policy.initializeCandidates(gctx, args, state, stream, logger);
  };
};

template <detail::HasOldInitializeCandidates T>
void oldToNewSignatureAdapter(const void* instance, const GeometryContext& gctx,
                              const NavigationArguments& args,
                              NavigationPolicyState& /*state*/,
                              AppendOnlyNavigationStream& stream,
                              const Logger& logger) {
  const auto* policy = static_cast<const T*>(instance);
  policy->initializeCandidates(gctx, args, stream, logger);
}
}  // namespace detail

/// Concept for a navigation policy
/// This exists so `initializeCandidates` can be a non-virtual method and we
/// still have a way to enforce it exists.
template <typename T>
concept NavigationPolicyConcept = requires {
  requires std::is_base_of_v<INavigationPolicy, T>;

  // Require either of the signatures to allow backwards compatibility
  requires(detail::HasOldInitializeCandidates<T> ||
           detail::HasNewInitializeCandidates<T>);
};

/// Base class for all navigation policies. The policy needs to be *connected*
/// to a delegate via a virtual method for it to become active. The update
/// method is not part of the class interface. The conventional `updateState`
/// method is only required for use with the navigation policy factory,
/// otherwise `connect` is free to connect any function.
class INavigationPolicy {
 public:
  /// Noop update function that is suitable as a default for default navigation
  /// delegates.
  static void noopInitializeCandidates(
      const GeometryContext& /*unused*/, const NavigationArguments& /*unused*/,
      NavigationPolicyState& /*unused*/,
      const AppendOnlyNavigationStream& /*unused*/, const Logger& /*unused*/) {
    // This is a noop
  }

  /// Virtual destructor so policies can be held through this base class.
  virtual ~INavigationPolicy() = default;

  /// Connect a policy with a delegate (usually a member of a volume).
  /// This method exists to allow a policy to ensure a non-virtual function is
  /// registered with the delegate.
  /// @param delegate The delegate to connect to
  virtual void connect(NavigationDelegate& delegate) const = 0;

  /// Convenience function to walk over all navigation policies. The default
  /// implementation just calls this on itself, while the @ref
  /// MultiNavigationPolicy will call it on all it's children.
  /// @param visitor The visitor function to call for each policy
  virtual void visit(
      const std::function<void(const INavigationPolicy&)>& visitor) const {
    visitor(*this);
  }

  /// Check if the policy is in a valid state for navigation
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param state The navigation policy state to check
  /// @param logger Logger for debug output
  /// @return True if the policy state is valid, false otherwise
  virtual bool isValid([[maybe_unused]] const GeometryContext& gctx,
                       [[maybe_unused]] const NavigationArguments& args,
                       [[maybe_unused]] NavigationPolicyState& state,
                       const Logger& logger) const {
    ACTS_VERBOSE("Default navigation policy isValid check. (always true)");
    return true;
  }

  struct EmptyState {};

  /// Create and initialize the state for this policy
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param stateManager The state manager to push the new state onto
  /// @param logger Logger for debug output
  virtual void createState([[maybe_unused]] const GeometryContext& gctx,
                           [[maybe_unused]] const NavigationArguments& args,
                           NavigationPolicyStateManager& stateManager,
                           const Logger& logger) const {
    ACTS_VERBOSE(
        "Default navigation policy state initialization. (empty state)");
    stateManager.pushState<EmptyState>();
  }

  /// One-time initialization step that probes whether this policy pushes only
  /// default (empty) states and caches the result (mutating internal state),
  /// allowing the navigator to skip the per-volume-entry state creation and the
  /// matching pop entirely for volumes with such policies. This is an eager
  /// initialization (deliberately not a lazy first-call latch): it is called
  /// for all attached policies at the end of Gen3 geometry construction, once
  /// every policy is in place. It must be called again should a policy be
  /// attached after that point; an uninitialized policy conservatively reports
  /// stateful via @ref isStateless, so skipping this step is always safe.
  ///
  /// @note This relies on the contract that the *defaultness* of the state a
  ///       navigation policy pushes does not depend on the geometry context
  ///       or the navigation arguments: a policy either always pushes a real
  ///       state (it has something to validate) or always pushes the default
  ///       empty state.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param logger Logger for debug output
  void initializeStatelessCache(const GeometryContext& gctx,
                                const Logger& logger = getDummyLogger()) {
    m_stateless = false;
    // Create the policy state once into a scratch manager and observe whether
    // the result is the default sentinel. Per the contract above, the
    // defaultness must not depend on the arguments, so probing with arbitrary
    // ones is valid.
    NavigationPolicyStateManager scratchManager;
    NavigationArguments args{.position = Vector3::Zero(),
                             .direction = Vector3::UnitZ()};
    createState(gctx, args, scratchManager, logger);
    m_stateless = scratchManager.currentState().isDefault();
  }

  /// Whether this policy is known to push only default (empty) states, see
  /// @ref initializeStatelessCache.
  /// @return true if the policy was probed and found stateless
  bool isStateless() const { return m_stateless; }

  /// Remove the state for this policy from the state manager
  /// @param stateManager The state manager to pop the state from
  /// @param logger Logger for debug output
  virtual void popState(NavigationPolicyStateManager& stateManager,
                        const Logger& logger) const {
    // By default, we didn't push anything, so we don't need to poop anything
    ACTS_VERBOSE("Default navigation policy pop state. (pops empty state)");
    stateManager.popState();  // This will be correct regardless of if a derived
                              // class pushed a concrete state type that's
                              // different from the default EmptyState.
  }

 protected:
  /// Internal helper function for derived classes that conform to the concept
  /// and have a conventional `updateState` method. Mainly used to save some
  /// boilerplate.
  /// @tparam T The type of the navigation policy
  /// @param delegate The delegate to connect to
  template <NavigationPolicyConcept T>
  void connectDefault(NavigationDelegate& delegate) const {
    // This cannot be a concept because we use it in CRTP below
    const auto* self = static_cast<const T*>(this);

    if constexpr (detail::HasNewInitializeCandidates<T>) {
      delegate.template connect<&T::initializeCandidates>(self);
    } else {
      // @TODO: Remove navigation policy signature compatibility eventually
      delegate.connect(&detail::oldToNewSignatureAdapter<T>, self);
    }
  }

 private:
  /// Whether this policy was probed and found to push only default states
  /// (see initializeStatelessCache). Conservative default: an unprobed
  /// policy is treated as stateful.
  bool m_stateless = false;
};

inline bool NavigationPolicyState::isDefault() const {
  if (empty()) {
    return true;  // no state at all -> nothing to validate
  }
  const NavigationPolicyStateStorage& stored = payload();
  return !stored || stored.is<INavigationPolicy::EmptyState>();
}

}  // namespace Acts
