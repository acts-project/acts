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

class NavigationPolicyState {
 public:
  template <typename T>
  T& as() {
    return std::any_cast<T&>(payload());
  }

  template <typename T>
  const T& as() const {
    return std::any_cast<T&>(payload());
  }

  NavigationPolicyState() = default;

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

class NavigationPolicyStateManager {
 public:
  template <typename T, typename... Args>
  std::pair<T&, NavigationPolicyState> pushState(Args&&... args) {
    std::any& state = m_stateStack.emplace_back();
    T& content = state.emplace<T>(std::forward<Args>(args)...);
    return std::pair<T&, NavigationPolicyState>{
        content, {*this, m_stateStack.size() - 1}};
  }

  friend class Navigator;

  NavigationPolicyState currentState() {
    if (m_stateStack.empty()) {
      return {};  // Emtpy state as sentinel
    }
    return NavigationPolicyState{*this, m_stateStack.size() - 1};
  }

  void popState() {
    if (m_stateStack.empty()) {
      throw std::runtime_error(
          "NavigationPolicyStateManager: Attempt to pop from empty stack");
    }
    m_stateStack.pop_back();
  }

 private:
  // We might want to extend this to a stack
  std::vector<std::any> m_stateStack;

  friend class NavigationPolicyState;
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

  virtual bool isValid(const GeometryContext& /*gctx*/,
                       const NavigationArguments /*args*/,
                       NavigationPolicyState& /*state*/,
                       const Logger& logger) const {
    ACTS_VERBOSE("Default navigation policy isValid check. (always true)");
    return true;
  }

  // MAKE SURE to also implement the pop method if you implement this one!
  virtual NavigationPolicyState createState(
      const GeometryContext& /*gctx*/, const NavigationArguments /*args*/,
      NavigationPolicyStateManager& /*stateManager*/,
      const Logger& logger) const {
    ACTS_VERBOSE("Default navigation policy state initialization. (noop)");
    return {};
  }

  virtual void popState(NavigationPolicyStateManager& /*stateManager*/,
                        const Logger& logger) const {
    // By default, we didn't push anything, so we don't need to poop anything
    ACTS_VERBOSE("Default navigation policy pop state. (noop)");
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
