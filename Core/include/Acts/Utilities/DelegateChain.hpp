// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/TypeList.hpp"

#include <concepts>
#include <tuple>
#include <utility>

namespace Acts {

template <typename, typename payload_types = TypeList<>>
class DelegateChain;

/// This class is a thin wrapper around a vector of delegates that are supposed
/// to be called in order.
template <typename R, typename... payload_types, typename... Args>
class DelegateChain<R(Args...), TypeList<payload_types...>> {
 public:
  using return_type =
      std::conditional_t<std::is_same_v<R, void>, void,
                         std::array<R, sizeof...(payload_types)>>;
  using tuple_type = std::tuple<payload_types...>;
  using delegate_type = Delegate<return_type(Args...)>;

  DelegateChain(std::unique_ptr<tuple_type> payloads)
      : m_payloads(std::move(payloads)) {}

  delegate_type& delegate() { return m_delegate; }

  /// The call operator that exposes the functionality of the @c DelegateChain type.
  /// @param args The arguments to call the contained functions with
  template <typename... Ts>
  auto operator()(Ts&&... args) const {
    return m_delegate(std::forward<Ts>(args)...);
  }

 private:
  delegate_type m_delegate;
  std::unique_ptr<tuple_type> m_payloads{};
};

template <typename, typename payload_types = TypeList<>, auto... Callables>
class DelegateChainFactory;

struct DelegateNoPayloadTag {};

// @TODO: Maybe add concept requirement for default initialization of R
template <typename R, typename... payload_types, auto... callables,
          typename... callable_args>
class DelegateChainFactory<R(callable_args...), TypeList<payload_types...>,
                           callables...> {
  using return_type = R;
  using chain_type =
      DelegateChain<R(callable_args...), TypeList<payload_types...>>;
  using tuple_type = typename chain_type::tuple_type;

 public:
  DelegateChainFactory() = default;
  DelegateChainFactory(std::tuple<payload_types...> payloads)
      : m_payloads(payloads) {}

  template <auto Callable, typename payload_type>
  auto add(payload_type&& payload) {
    std::tuple<payload_types..., payload_type> payloads =
        std::tuple_cat(m_payloads, std::make_tuple(payload));

    return DelegateChainFactory<R(callable_args...),
                                TypeList<payload_types..., payload_type>,
                                callables..., Callable>{payloads};
  }

  template <auto Callable>
  auto add() {
    std::tuple<payload_types..., DelegateNoPayloadTag> payloads =
        std::tuple_cat(m_payloads, std::make_tuple(DelegateNoPayloadTag{}));

    return DelegateChainFactory<
        R(callable_args...), TypeList<payload_types..., DelegateNoPayloadTag>,
        callables..., Callable>{payloads};
  }

  template <std::size_t I, std::size_t J, auto head, auto... tail>
  static constexpr auto findCallable() {
    if constexpr (I == J) {
      return head;
    } else {
      return findCallable<I, J + 1, tail...>();
    }
  }

  template <std::size_t I = 0, typename result_ptr>
  static constexpr auto invoke(result_ptr result, const tuple_type* payloads,
                               callable_args... args) {
    const auto& callable = findCallable<I, 0, callables...>();

    if constexpr (!std::is_same_v<std::tuple_element_t<I, tuple_type>,
                                  DelegateNoPayloadTag>) {
      auto payload = std::get<I>(*payloads);

      if constexpr (!std::is_same_v<result_ptr, std::nullptr_t>) {
        std::get<I>(*result) = std::invoke(
            callable, payload, std::forward<callable_args>(args)...);
      } else {
        std::invoke(callable, payload, std::forward<callable_args>(args)...);
      }

    } else {
      if constexpr (!std::is_same_v<result_ptr, std::nullptr_t>) {
        std::get<I>(*result) =
            std::invoke(callable, std::forward<callable_args>(args)...);
      } else {
        std::invoke(callable, std::forward<callable_args>(args)...);
      }
    }

    if constexpr (I < sizeof...(payload_types) - 1) {
      invoke<I + 1>(result, payloads, std::forward<callable_args>(args)...);
    }
  }

  chain_type build() {
    auto payloads = std::make_unique<tuple_type>(m_payloads);
    const tuple_type* payloadsPtr = payloads.get();
    chain_type chain{std::move(payloads)};

    typename chain_type::delegate_type::function_type function =
        [](const void* payloadPack, callable_args... args) ->
        typename chain_type::return_type {
          const auto* tuplePtr = static_cast<const tuple_type*>(payloadPack);
          if constexpr (std::is_same_v<return_type, void>) {
            invoke(nullptr, tuplePtr, std::forward<callable_args>(args)...);
          } else {
            typename chain_type::return_type result{};
            invoke(&result, tuplePtr, std::forward<callable_args>(args)...);
            return result;
          }
        };

    chain.delegate().connect(function, payloadsPtr);
    return chain;
  }

 private:
  tuple_type m_payloads{};
};

}  // namespace Acts
