// SPDX-License-Identifier: MIT
// Copyright 2018 Moritz Kiehn
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/// \file
/// \brief   Command dispatcher to register functions and call them by name
/// \author  Moritz Kiehn <msmk@cern.ch>
/// \date    2018-02-20

#pragma once

#include <cassert>
#include <cstdint>
#include <functional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace dfe {

/// Variable-type value object a.k.a. a poor mans std::variant.
class Variable final {
public:
  /// Supported value types.
  enum class Type { Empty, Boolean, Integer, Float, String };

  Variable() : m_type(Type::Empty) {}
  Variable(Variable&& v) { *this = std::move(v); }
  Variable(const Variable& v) { *this = v; }
  explicit Variable(std::string&& s)
    : m_string(std::move(s)), m_type(Type::String) {}
  explicit Variable(const std::string& s) : Variable(std::string(s)) {}
  explicit Variable(const char* s) : Variable(std::string(s)) {}
  // suppport all possible integer types
  template<typename I, typename = std::enable_if_t<std::is_integral<I>::value>>
  explicit Variable(I integer)
    : m_integer(static_cast<int64_t>(integer)), m_type(Type::Integer) {}
  explicit Variable(double d) : m_float(d), m_type(Type::Float) {}
  explicit Variable(float f) : Variable(static_cast<double>(f)) {}
  explicit Variable(bool b) : m_boolean(b), m_type(Type::Boolean) {}
  ~Variable() = default;

  Variable& operator=(Variable&& v);
  Variable& operator=(const Variable& v);

  /// Parse a string into a value of the requested type.
  static Variable parse_as(const std::string& str, Type type);

  /// In a boolean context a variable is false if it does not contain a value.
  ///
  /// \warning This is not the value of the stored boolean.
  constexpr bool operator!() const { return m_type == Type::Empty; }
  /// \see operator!()
  constexpr explicit operator bool() const { return !!(*this); }
  /// The type of the currently stored value.
  constexpr Type type() const { return m_type; }
  /// Get value of the variable as a specific type.
  ///
  /// \exception std::invalid_argument if the requested type is incompatible
  template<typename T>
  auto as() const;

private:
  template<typename T>
  struct Converter;
  template<typename I>
  struct IntegerConverter;

  union {
    int64_t m_integer;
    double m_float;
    bool m_boolean;
  };
  // std::string has non-trivial constructor; cannot store by value in union
  // TODO 2018-11-29 msmk: more space-efficient string storage
  std::string m_string;
  Type m_type;

  friend std::ostream& operator<<(std::ostream& os, const Variable& v);
};

/// A simple command dispatcher.
///
/// You can register commands and call them by name.
class Dispatcher {
public:
  /// The native dispatcher function interface.
  using Interface = std::function<Variable(const std::vector<Variable>&)>;

  /// Register a native dispatcher function.
  ///
  /// \param name      Unique function name
  /// \param func      Function object
  /// \param arg_types Arguments types
  /// \param help      Optional help text
  void add(
    std::string name, Interface&& func, std::vector<Variable::Type>&& arg_types,
    std::string help = std::string());
  /// Register a function with arbitrary arguments.
  ///
  /// The return type and the argument types must be compatible with `Variable`.
  template<typename R, typename... Args>
  void add(
    std::string name, std::function<R(Args...)>&& func,
    std::string help = std::string());
  template<typename R, typename... Args>
  void add(
    std::string name, R (*func)(Args...), std::string help = std::string());
  template<typename T, typename R, typename... Args>
  void add(
    std::string name, R (T::*member_func)(Args...), T* t,
    std::string help = std::string());

  /// Call a command with arbitrary arguments.
  template<typename... Args>
  Variable call(const std::string& name, Args&&... args);
  /// Call a command with arguments parsed from strings into the expected types.
  Variable call_parsed(
    const std::string& name, const std::vector<std::string>& args);
  /// Call a command using the native argument encoding.
  Variable call_native(
    const std::string& name, const std::vector<Variable>& args);

  /// Return a list of registered commands.
  std::vector<std::string> commands() const;
  /// Return the help text for the command.
  const std::string& help(const std::string& name) const;

private:
  struct Command {
    Interface func;
    std::vector<Variable::Type> argument_types;
    std::string help;
  };
  std::unordered_map<std::string, Command> m_commands;
};

// implementation Variable

inline Variable
Variable::parse_as(const std::string& str, Type type) {
  if (type == Type::Boolean) {
    return Variable((str == "true"));
  } else if (type == Type::Integer) {
    return Variable(std::stoll(str));
  } else if (type == Type::Float) {
    return Variable(std::stod(str));
  } else if (type == Type::String) {
    return Variable(str);
  } else {
    return Variable();
  }
}

inline std::ostream&
operator<<(std::ostream& os, const Variable& v) {
  if (v.type() == Variable::Type::Boolean) {
    os << (v.m_boolean ? "true" : "false");
  } else if (v.m_type == Variable::Type::Integer) {
    os << v.m_integer;
  } else if (v.m_type == Variable::Type::Float) {
    os << v.m_float;
  } else if (v.m_type == Variable::Type::String) {
    os << v.m_string;
  }
  return os;
}

inline Variable&
Variable::operator=(Variable&& v) {
  // handle `x = std::move(x)`
  if (this == &v) {
    return *this;
  }
  if (v.m_type == Type::Boolean) {
    m_boolean = v.m_boolean;
  } else if (v.m_type == Type::Integer) {
    m_integer = v.m_integer;
  } else if (v.m_type == Type::Float) {
    m_float = v.m_float;
  } else if (v.m_type == Type::String) {
    m_string = std::move(v.m_string);
  }
  m_type = v.m_type;
  return *this;
}

inline Variable&
Variable::operator=(const Variable& v) {
  if (v.m_type == Type::Boolean) {
    m_boolean = v.m_boolean;
  } else if (v.m_type == Type::Integer) {
    m_integer = v.m_integer;
  } else if (v.m_type == Type::Float) {
    m_float = v.m_float;
  } else if (v.m_type == Type::String) {
    m_string = v.m_string;
  }
  m_type = v.m_type;
  return *this;
}

template<>
struct Variable::Converter<bool> {
  static constexpr Type type() { return Type::Boolean; }
  static constexpr bool as_t(const Variable& v) { return v.m_boolean; }
};
template<>
struct Variable::Converter<float> {
  static constexpr Type type() { return Type::Float; }
  static constexpr float as_t(const Variable& v) {
    return static_cast<float>(v.m_float);
  }
};
template<>
struct Variable::Converter<double> {
  static constexpr Type type() { return Type::Float; }
  static constexpr double as_t(const Variable& v) { return v.m_float; }
};
template<>
struct Variable::Converter<std::string> {
  static constexpr Type type() { return Type::String; }
  static constexpr const std::string& as_t(const Variable& v) {
    return v.m_string;
  }
};
template<typename I>
struct Variable::IntegerConverter {
  static constexpr Type type() { return Type::Integer; }
  static constexpr I as_t(const Variable& v) {
    return static_cast<I>(v.m_integer);
  }
};
template<>
struct Variable::Converter<int8_t> : Variable::IntegerConverter<int8_t> {};
template<>
struct Variable::Converter<int16_t> : Variable::IntegerConverter<int16_t> {};
template<>
struct Variable::Converter<int32_t> : Variable::IntegerConverter<int32_t> {};
template<>
struct Variable::Converter<int64_t> : Variable::IntegerConverter<int64_t> {};
template<>
struct Variable::Converter<uint8_t> : Variable::IntegerConverter<uint8_t> {};
template<>
struct Variable::Converter<uint16_t> : Variable::IntegerConverter<uint16_t> {};
template<>
struct Variable::Converter<uint32_t> : Variable::IntegerConverter<uint32_t> {};
template<>
struct Variable::Converter<uint64_t> : Variable::IntegerConverter<uint64_t> {};

template<typename T>
inline auto
Variable::as() const {
  if (m_type != Variable::Converter<T>::type()) {
    throw std::invalid_argument(
      "Requested type is incompatible with stored type");
  }
  return Variable::Converter<T>::as_t(*this);
}

// implementation Dispatcher

namespace dispatcher_impl {
namespace {

// Wrap a function that returns a value
template<typename R, typename... Args>
struct InterfaceWrappper {
  std::function<R(Args...)> func;

  Variable operator()(const std::vector<Variable>& args) {
    return call(args, std::index_sequence_for<Args...>());
  }
  template<std::size_t... I>
  Variable call(const std::vector<Variable>& args, std::index_sequence<I...>) {
    return Variable(func(args.at(I).as<typename std::decay_t<Args>>()...));
  }
};

// Wrap a function that does not return anything
template<typename... Args>
struct InterfaceWrappper<void, Args...> {
  std::function<void(Args...)> func;

  Variable operator()(const std::vector<Variable>& args) {
    return call(args, std::index_sequence_for<Args...>());
  }
  template<std::size_t... I>
  Variable call(const std::vector<Variable>& args, std::index_sequence<I...>) {
    func(args.at(I).as<typename std::decay_t<Args>>()...);
    return Variable();
  }
};

template<typename R, typename... Args>
inline Dispatcher::Interface
make_wrapper(std::function<R(Args...)>&& function) {
  return InterfaceWrappper<R, Args...>{std::move(function)};
}

template<typename R, typename... Args>
std::vector<Variable::Type>
make_types(const std::function<R(Args...)>&) {
  return {Variable(std::decay_t<Args>()).type()...};
}

} // namespace
} // namespace dispatcher_impl

inline void
Dispatcher::add(
  std::string name, Dispatcher::Interface&& func,
  std::vector<Variable::Type>&& arg_types, std::string help) {
  if (name.empty()) {
    throw std::invalid_argument("Can not register command with empty name");
  }
  if (m_commands.count(name)) {
    throw std::invalid_argument(
      "Can not register command '" + name + "' more than once");
  }
  m_commands[std::move(name)] =
    Command{std::move(func), std::move(arg_types), std::move(help)};
}

template<typename R, typename... Args>
inline void
Dispatcher::add(
  std::string name, std::function<R(Args...)>&& func, std::string help) {
  auto args = dispatcher_impl::make_types(func);
  add(
    std::move(name), dispatcher_impl::make_wrapper(std::move(func)),
    std::move(args), std::move(help));
}

template<typename R, typename... Args>
inline void
Dispatcher::add(std::string name, R (*func)(Args...), std::string help) {
  assert(func && "Function pointer must be non-null");
  add(std::move(name), std::function<R(Args...)>(func), std::move(help));
}

template<typename T, typename R, typename... Args>
inline void
Dispatcher::add(
  std::string name, R (T::*member_func)(Args...), T* t, std::string help) {
  assert(member_func && "Member function pointer must be non-null");
  assert(t && "Object pointer must be non-null");
  add(
    std::move(name), std::function<R(Args...)>([=](Args... args) {
      return (t->*member_func)(args...);
    }),
    std::move(help));
}

inline Variable
Dispatcher::call_native(
  const std::string& name, const std::vector<Variable>& args) {
  auto cmd = m_commands.find(name);
  if (cmd == m_commands.end()) {
    throw std::invalid_argument("Unknown command '" + name + "'");
  }
  if (args.size() != cmd->second.argument_types.size()) {
    throw std::invalid_argument("Invalid number of arguments");
  }
  return cmd->second.func(args);
}

inline Variable
Dispatcher::call_parsed(
  const std::string& name, const std::vector<std::string>& args) {
  // dont reuse call_native since we need to have access to the command anyways
  auto cmd = m_commands.find(name);
  if (cmd == m_commands.end()) {
    throw std::invalid_argument("Unknown command '" + name + "'");
  }
  if (args.size() != cmd->second.argument_types.size()) {
    throw std::invalid_argument("Invalid number of arguments");
  }
  // convert string arguments into Variable values
  std::vector<Variable> vargs;
  for (std::size_t i = 0; i < args.size(); ++i) {
    vargs.push_back(Variable::parse_as(args[i], cmd->second.argument_types[i]));
  }
  return cmd->second.func(vargs);
}

template<typename... Args>
inline Variable
Dispatcher::call(const std::string& name, Args&&... args) {
  return call_native(
    name, std::vector<Variable>{Variable(std::forward<Args>(args))...});
}

inline std::vector<std::string>
Dispatcher::commands() const {
  std::vector<std::string> cmds;

  for (const auto& cmd : m_commands) {
    cmds.emplace_back(cmd.first);
  }
  return cmds;
}

inline const std::string&
Dispatcher::help(const std::string& name) const {
  return m_commands.at(name).help;
}

} // namespace dfe
