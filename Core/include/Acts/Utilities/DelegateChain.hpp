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

template <typename, typename payload_types = TypeList<>, auto... Callables>
class DelegateChainFactory;

// @TODO: Maybe add concept requirement for default initialization of R
template <typename R, typename... payload_types, auto... callables,
          typename... callable_args>
class DelegateChainFactory<R(callable_args...), TypeList<payload_types...>,
                           callables...> {
  using return_type =
      std::conditional_t<std::is_same_v<R, void>, void,
                         std::array<R, sizeof...(payload_types)>>;
  using delegate_type =
      Delegate<return_type(callable_args...), void, DelegateType::Owning>;
  using tuple_type = std::tuple<payload_types...>;

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
    std::tuple<payload_types..., std::nullptr_t> payloads =
        std::tuple_cat(m_payloads, std::make_tuple(std::nullptr_t{}));

    return DelegateChainFactory<R(callable_args...),
                                TypeList<payload_types..., std::nullptr_t>,
                                callables..., Callable>{payloads};
  }

  struct DispatchBlock {
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
                                    std::nullptr_t>) {
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

    DispatchBlock(tuple_type payloads) : m_payloads(std::move(payloads)) {}

    tuple_type m_payloads{};

    auto dispatch(callable_args... args) const {
      if constexpr (std::is_same_v<R, void>) {
        invoke(nullptr, &m_payloads, std::forward<callable_args>(args)...);
      } else {
        return_type result{};
        invoke(&result, &m_payloads, std::forward<callable_args>(args)...);
        return result;
      }
    }
  };

  delegate_type build() {
    auto block = std::make_unique<const DispatchBlock>(m_payloads);
    delegate_type delegate;
    delegate.template connect<&DispatchBlock::dispatch>(std::move(block));
    return delegate;
  }

  void store(delegate_type& delegate) {
    auto block = std::make_unique<const DispatchBlock>(m_payloads);
    delegate.template connect<&DispatchBlock::dispatch>(std::move(block));
  }

 private:
  tuple_type m_payloads{};
};

}  // namespace Acts
