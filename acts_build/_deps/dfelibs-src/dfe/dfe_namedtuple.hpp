// SPDX-License-Identifier: MIT
// Copyright 2015,2018-2020 Moritz Kiehn
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
/// \brief   Allow structs to be accessed like std::tuple
/// \author  Moritz Kiehn <msmk@cern.ch>
/// \date    2015-06-00, Initial version
/// \date    2018-02-09, Major rework
/// \date    2019-09-09, Split i/o components into separate libraries
/// \date    2020-10-13, Add support for direct use of `get<I>(named_tuple)`

#pragma once

#include <array>
#include <cassert>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>

/// Enable tuple-like access and conversion for selected class/struct members.
///
/// This allows access to the selected members via `.get<I>()` or `get<I>(...)`,
/// conversion to equivalent `std::tuple<...>` via implicit conversion or
/// explicitely via `.tuple()`,  and assignment from equivalent tuples.
/// The names can be accessed via `::names()`.
#define DFE_NAMEDTUPLE(name, ...) \
  using Tuple = decltype(::std::make_tuple(__VA_ARGS__)); \
  static ::std::array<::std::string, ::std::tuple_size<Tuple>::value> \
  names() { \
    return ::dfe::namedtuple_impl::unstringify< \
      ::std::tuple_size<Tuple>::value>((#__VA_ARGS__)); \
  } \
  template<typename... U> \
  name& operator=(const ::std::tuple<U...>& other) { \
    ::std::tie(__VA_ARGS__) = other; \
    return *this; \
  } \
  template<typename... U> \
  name& operator=(::std::tuple<U...>&& other) { \
    ::std::tie(__VA_ARGS__) = ::std::forward<std::tuple<U...>>(other); \
    return *this; \
  } \
  operator Tuple() const { return ::std::make_tuple(__VA_ARGS__); } \
  Tuple tuple() const { return ::std::make_tuple(__VA_ARGS__); } \
  template<std::size_t I> \
  constexpr ::std::tuple_element_t<I, Tuple>& get() { \
    return ::std::get<I>(std::tie(__VA_ARGS__)); \
  } \
  template<::std::size_t I> \
  constexpr const ::std::tuple_element_t<I, Tuple>& get() const { \
    return ::std::get<I>(std::tie(__VA_ARGS__)); \
  } \
  template<::std::size_t I> \
  friend constexpr ::std::tuple_element_t<I, Tuple>& get(name& nt) { \
    return nt.template get<I>(); \
  } \
  template<::std::size_t I> \
  friend constexpr const ::std::tuple_element_t<I, Tuple>& get( \
    const name& nt) { \
    return nt.template get<I>(); \
  } \
  friend inline ::std::ostream& operator<<(::std::ostream& os, const name& nt) \
    __attribute__((unused)) { \
    return ::dfe::namedtuple_impl::print_tuple( \
      os, nt.names(), nt.tuple(), \
      ::std::make_index_sequence<::std::tuple_size<Tuple>::value>{}); \
  }

// implementation helpers
namespace dfe {
namespace namedtuple_impl {

// Reverse macro stringification.
//
// Splits a string of the form `a, b, c` into components a, b, and c.
template<std::size_t N>
constexpr std::array<std::string, N>
unstringify(const char* str) {
  assert(str and "Input string must be non-null");

  std::array<std::string, N> out;

  for (std::size_t idx = 0; idx < N; ++idx) {
    // skip leading whitespace
    while ((*str != '\0') and (*str == ' ')) {
      ++str;
    }
    // find the next separator or end-of-string
    const char* sep = str;
    while ((*sep != '\0') and (*sep != ',')) {
      ++sep;
    }
    // store component w/o the separator
    out[idx].assign(str, sep - str);
    // we can quit as soon as we reached the end of the input
    if (*sep == '\0') {
      break;
    }
    // start search for next component after the separator
    str = ++sep;
  }
  // TODO handle inconsistent number of entries? can it occur in expected use?
  return out;
}

// modified from http://stackoverflow.com/a/6245777
template<typename Names, typename Values, std::size_t... I>
inline std::ostream&
print_tuple(
  std::ostream& os, const Names& n, const Values& v,
  std::index_sequence<I...>) {
  // we want to execute some expression for every entry in the index pack. this
  // requires a construction that can take a variable number of arguments into
  // which we can unpack the indices. inside a function, constructing an
  // array with an initializer list will do the job, i.e. we will effectively
  // create the following statement
  //
  //     int x[] = {...};
  //
  // since we do not care about the actual values within the array, the
  // initializer list is cast twice: once to the array type and then to void.
  // this ignores the actual values and silences warnings about unused
  // variables. to get the correct initializer list syntax, the array type
  // must be typedef'd as a single type. what we get is
  //
  //     using Vacuum = int[];
  //     (void)Vacuum{...};
  //
  // in order for this to work, the expression that we want to instantiate
  // needs to evaluate to the element type of the array (here: `int`). this can
  // be done with the comma operator (yep, `,` is a weird but helpful operator)
  // for arbitrary expressions. `(<expr1>, <expr2>)` executes both expressions
  // but evaluates only to the return value of the second expression. thus,
  // `(<expr>, 0)` executes `<expr>` but always evalutes to an integer of value
  // zero. if <expr> uses the index pack variable `I` in the following setup
  //
  //     (void)Vacuum{(<expr>, 0)...};
  //
  // it is instantiatied for each element within the pack (with appropriate ,
  // placements). thus, effectively looping over every entry in the pack and
  // calling <expr> for each (here: printing to os);
  using std::get;
  using Vacuum = int[];
  (void)Vacuum{
    (os << ((0 < I) ? " " : "") << get<I>(n) << "=" << get<I>(v), 0)...};
  return os;
}

} // namespace namedtuple_impl
} // namespace dfe
