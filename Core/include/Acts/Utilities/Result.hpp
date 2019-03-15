#pragma once

#include <system_error>
#include <type_traits>
#include <utility>
#include <variant>

namespace Acts {

template <typename T, typename E = std::error_code>
class Result
{
  static constexpr bool unambiguous
      = (!std::
             is_same_v<T,
                       E> && !std::is_constructible_v<T, E> && !std::is_constructible_v<T, E>);

  Result(std::variant<T, E>&& var) : m_var(var) {}

public:
  Result() = delete;

  Result(const Result<T, E>& other) = delete;
  Result<T, E>&
  operator=(const Result<T, E>& other)
      = delete;

  Result(Result<T, E>&& other) = default;
  Result<T, E>&
  operator=(Result<T, E>&& other)
      = default;

  template <typename T2, typename = std::enable_if_t<unambiguous>>
  Result(T2 value) noexcept : m_var(std::move(value))
  {
  }

  static Result<T, E>
  success(T value)
  {
    return Result<T, E>(
        std::variant<T, E>{std::in_place_index<0>, std::move(value)});
  }

  static Result<T, E>
  failure(E error)
  {
    return Result<T, E>(
        std::variant<T, E>{std::in_place_index<1>, std::move(error)});
  }

  bool
  ok() const noexcept
  {
    return m_var.index() == 0;
  }

  T& operator*() noexcept { return std::get<T>(m_var); }

  E&
  error() noexcept
  {
    return std::get<E>(m_var);
  }

  T
  unwrap()
  {
    if (m_var.index() != 0) {
      throw std::runtime_error("Unwrap called on error value");
    }

    return std::move(std::get<T>(m_var));
  }

private:
  std::variant<T, E> m_var;
};
}
