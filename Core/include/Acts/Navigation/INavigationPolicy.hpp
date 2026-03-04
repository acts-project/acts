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

/// Wrapper class for a navigation policy state stored in the state manager.
/// This class provides type-safe access to the underlying state through
/// the `as<T>()` method. The state is stored as a type-erased std::any
/// in the NavigationPolicyStateManager.
class NavigationPolicyState {
 public:
  /// Cast the state to a specific type
  /// @tparam T The type to cast to
  /// @return Reference to the state as type T
  template <typename T>
  T& as() {
    return std::any_cast<T&>(payload());
  }

  /// Cast the state to a specific type (const version)
  /// @tparam T The type to cast to
  /// @return Const reference to the state as type T
  template <typename T>
  const T& as() const {
    return std::any_cast<T&>(payload());
  }

  NavigationPolicyState() = default;

  /// Check if the state is empty (not associated with a manager)
  /// @return True if the state is empty, false otherwise
  bool empty() const { return m_manager == nullptr; }

 private:
  NavigationPolicyState(NavigationPolicyStateManager& manager,
                        std::size_t index)
      : m_manager(&manager), m_index(index) {}

  std::any& payload();
  const std::any& payload() const;

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
    std::any& state = m_stateStack.emplace_back();
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
  std::vector<std::any> m_stateStack;
};

inline std::any& NavigationPolicyState::payload() {
  if (m_manager == nullptr) {
    throw std::runtime_error(
        "NavigationPolicyState: Attempt to access empty payload");
  }
  return m_manager->m_stateStack.at(m_index);
}

inline const std::any& NavigationPolicyState::payload() const {
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
};

}  // namespace Acts
