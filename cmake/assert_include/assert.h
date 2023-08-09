#pragma once

#include <stdexcept>
namespace {
inline void __assert_fail(const char* expression, const char* file, int line,
                          const char* function) {
  throw std::runtime_error{"Condition failed: " + std::string{expression} +
                           " at " + std::string{file} + ":" +
                           std::to_string(line) + ": " + function};
}
}  // namespace

#define assert(expr)       \
  (static_cast<bool>(expr) \
       ? (void)(0)         \
       : __assert_fail(#expr, __FILE__, __LINE__, __PRETTY_FUNCTION__))
