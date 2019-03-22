// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <optional>
#include <variant>

namespace Acts {

template <typename Derived, typename... States>
class FiniteStateMachine
{
public:
  struct Terminated
  {
  };
  using StateVariant = std::variant<Terminated, States...>;

protected:
  using self_type = FiniteStateMachine<Derived, States...>;

  using event_return = std::optional<StateVariant>;

public:
  FiniteStateMachine()
    : m_state(typename std::tuple_element<0, std::tuple<States...>>::type{}){};

  FiniteStateMachine(StateVariant state) : m_state(std::move(state)){};

  const StateVariant&
  getState() const
  {
    return m_state;
  }

  StateVariant&
  getState()
  {
    return m_state;
  }

  void
  setState(StateVariant state)
  {
    Derived& child = static_cast<Derived&>(*this);

    // call on exit function
    std::visit([&](auto& s) { child.on_exit(s); }, m_state);

    // no change state
    m_state = std::move(state);

    // call on enter function
    std::visit([&](auto& s) { child.on_enter(s); }, m_state);
  }

  template <typename S>
  bool
  is(const S& /*state*/) const noexcept
  {
    if (std::get_if<S>(&m_state)) { return true; }
    return false;
  }

  bool
  terminated() const noexcept
  {
    return is(Terminated{});
  }

  template <typename Event, typename... Args>
  event_return
  process_event(Event&& event, Args&&... args)
  {

    Derived& child     = static_cast<Derived&>(*this);
    auto     new_state = std::visit(
        [&](auto& s) -> std::optional<StateVariant> {
          return child.on_event(
              s, std::forward<Event>(event), std::forward<Args...>(args)...);
        },
        m_state);
    return std::move(new_state);
  }

  template <typename Event, typename... Args>
  void
  dispatch(Event&& event, Args&&... args)
  {
    auto new_state = process_event(std::forward<Event>(event),
                                    std::forward<Args...>(args)...);
    if (new_state) { setState(std::move(*new_state)); }
  }

private:
  StateVariant m_state;
};

}  // namespace Acts
