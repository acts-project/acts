// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <system_error>
#include <type_traits>
#include <utility>
#include <variant>

namespace Acts {

/**
 * Class which encapsulates either a valid result, or an error
 * @tparam T The valid result value
 * @tparam E The error, defaults to `std::error_code`
 */
template <typename T, typename E = std::error_code>
class Result
{
  /**
   * Private constructor which accepts an external variant.
   * This is used by the factory static methods to set up
   * the variant unambiguously in all cases.
   */
  Result(std::variant<T, E>&& var) : m_var(var) {}

public:
  /**
   * Default construction is disallowed.
   */
  Result() = delete;

  /**
   * Copy construction is disallowed
   */
  Result(const Result<T, E>& other) = delete;

  /**
   * Assignment is disallowed
   */
  Result<T, E>&
  operator=(const Result<T, E>& other)
      = delete;

  /**
   * Move construction is allowed
   */
  Result(Result<T, E>&& other) : m_var(std::move(other.m_var)){};

  /**
   * Move assignment is allowed
   */
  Result<T, E>&
  operator=(Result<T, E>&& other)
  {
    m_var = std::move(other.m_var);
    return *this;
  }

  /**
   * @brief Constructor from arbitrary value
   * This constructor allows construction from any value. This constructor is
   * only enabled if T and E are unambiguous, meaning the cannot be implicitly
   * converted and there is T cannot be constructed from E and vice-versa.
   * This means that when this is invoked, the value can be propagated to the
   * underlying variant, and the assignment will be correct, and error will be
   * an error, and a value will be a value.
   * @note If T and E are ambigious, use the `success` and `failure` static
   * factory methods.
   * @tparam T2 Type of the potential assignment
   * @param value The potential value, could be an actual valid value or an
   * error.
   */
  template <
      typename T2,
      typename _E = E,
      typename _T = T,
      typename    = std::
          enable_if_t<!std::
                          is_same_v<_T,
                                    _E> && !std::is_constructible_v<_T, _E> && !std::is_constructible_v<_T, _E>>>
  Result(T2 value) noexcept : m_var(std::move(value))
  {
  }

  /**
   * Static helper factory which forces assignment as valid value.
   * @param value The valid value to assign. Will not be converted to E.
   * @return Initialized result object
   */
  static Result<T, E>
  success(T value)
  {
    return Result<T, E>(
        std::variant<T, E>{std::in_place_index<0>, std::move(value)});
  }

  /**
   * Static helper factory which forces assignment as an error.
   * @param value The error to assign. Will not be converted to T.
   * @return Initialized result object
   */
  static Result<T, E>
  failure(E error)
  {
    return Result<T, E>(
        std::variant<T, E>{std::in_place_index<1>, std::move(error)});
  }

  /**
   * Checks whether this result contains a valid value, and no error.
   * @return bool Whether result contains an error or not.
   */
  bool
  ok() const noexcept
  {
    return m_var.index() == 0;
  }

  /**
   * Returns a reference into the variant to the valid value.
   * @note If `!res.ok()`, this method will abort (noexcept)
   * @return Reference to value stored in the variant.
   */
  T& operator*() noexcept { return std::get<T>(m_var); }

  /**
   * Returns a reference to the error stored in the variant.
   * @note If `res.ok()` this method will abort (noexcept)
   * @return Reference to the error
   */
  E&
      error()
      & noexcept
  {
    return std::get<E>(m_var);
  }

  /**
   * Returns the error by-value.
   * @note If `res.ok()` this method will abort (noexcept)
   * @return The error
   */
  E
      error()
      && noexcept
  {
    return std::move(std::get<E>(m_var));
  }

  /**
   * Retrieves the valid value from the result object.
   * @note This is the lvalue version, returns a reference to the value
   * @return The valid value as a reference
   */
  T&
  value() &
  {
    if (m_var.index() != 0) {
      if
        constexpr(std::is_same_v<E, std::error_code>)
        {
          std::stringstream ss;
          const auto&       e = std::get<E>(m_var);
          ss << "Value called on error value: " << e.category().name() << ": "
             << e.message() << " [" << e.value() << "]";
          throw std::runtime_error(ss.str());
        }
      else {
        throw std::runtime_error("Value called on error value");
      }
    }

    return std::get<T>(m_var);
  }

  /**
   * Retrieves the valid value from the result object.
   * @note This is the rvalue version, returns the value
   * by-value and moves out of the variant.
   * @return The valid value by value, moved out of the variant.
   */
  T
  value() &&
  {
    if (m_var.index() != 0) {
      if
        constexpr(std::is_same_v<E, std::error_code>)
        {
          std::stringstream ss;
          const auto&       e = std::get<E>(m_var);
          ss << "Value called on error value: " << e.category().name() << ": "
             << e.message() << " [" << e.value() << "]";
          throw std::runtime_error(ss.str());
        }
      else {
        throw std::runtime_error("Value called on error value");
      }
    }

    return std::move(std::get<T>(m_var));
  }

private:
  std::variant<T, E> m_var;
};
}  // namespace Acts
