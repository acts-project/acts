// SPDX-License-Identifier: MIT

/// \file
/// \brief Demonstrate dfe::Dispatcher functionality

#include <cstdlib>
#include <iostream>

#include <dfe/dfe_dispatcher.hpp>

// example functions

void
func_noreturn(int x, double f) {
  std::cout << "free function w/o return: x=" << x << " f=" << f << '\n';
}

double
func_return(int x, double f) {
  std::cout << "free function w/ return: x=" << x << " f=" << f << '\n';
  return x + f;
}

dfe::Variable
func_native(const std::vector<dfe::Variable>& args) {
  std::cout << "native w/ " << args.size() << " arguments\n";
  std::string ret;
  for (const auto& arg : args) {
    ret += arg.as<std::string>();
  }
  return dfe::Variable(ret);
}

struct WithFunctions {
  float x;

  float member_add(float y) {
    std::cout << "member add x=" << x << " y=" << y << '\n';
    return x + y;
  }
  static float static_add(float a, float b) {
    std::cout << "static add a=" << a << " b=" << b << '\n';
    return a + b;
  }
};

int
main() {
  using Type = dfe::Variable::Type;
  using dfe::Dispatcher;

  Dispatcher dispatch;

  // add functions to dispatcher
  dispatch.add("noreturn", func_noreturn, "does something");
  dispatch.add("return", func_return, "returns something");
  // native function different number of arguments
  dispatch.add("native1", func_native, {Type::String});
  dispatch.add(
    "native3", func_native, {Type::String, Type::String, Type::String});
  dispatch.add(
    "native5", func_native,
    {Type::String, Type::String, Type::String, Type::String, Type::String});
  WithFunctions adder = {5.5};
  dispatch.add("member_add", &WithFunctions::member_add, &adder);
  dispatch.add("static_add", WithFunctions::static_add, "a static member");

  // list registered functions
  std::cout << "registered commands:\n";
  for (const auto& cmd : dispatch.commands()) {
    std::cout << "  " << cmd << '\n';
    std::cout << "    " << dispatch.help(cmd) << '\n';
  }

  // call functions by name
  std::cout << dispatch.call("noreturn", 1, 1.24) << '\n';
  std::cout << dispatch.call("return", 1, 1.24) << '\n';
  std::cout << dispatch.call("native1", "x") << '\n';
  std::cout << dispatch.call("native3", "x", "y", "z") << '\n';
  std::cout << dispatch.call_parsed("native3", {"x", "y", "z"}) << '\n';
  std::cout << dispatch.call_parsed("native5", {{"x", "y", "z", "1", "d"}})
            << '\n';
  std::cout << dispatch.call("member_add", 1.2) << '\n';
  std::cout << dispatch.call("static_add", 4.2, 2.3) << '\n';

  return EXIT_SUCCESS;
}
