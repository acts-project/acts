/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace traccc::opts {

template <typename TYPE, std::size_t SIZE>
value_array<TYPE, SIZE>& value_array<TYPE, SIZE>::operator*=(TYPE factor) {

    for (auto& v : *this) {
        v *= factor;
    }
    return *this;
}

template <typename TYPE, std::size_t SIZE>
value_array<TYPE, SIZE>& value_array<TYPE, SIZE>::operator/=(TYPE factor) {

    for (auto& v : *this) {
        v /= factor;
    }
    return *this;
}

namespace details {

/// Separator character for the float arguments
static constexpr std::string_view float_array_separator = ":";

}  // namespace details

template <typename TYPE, std::size_t SIZE>
std::istream& operator>>(std::istream& is, value_array<TYPE, SIZE>& values) {

    // Read the argument as a single string.
    std::string arg;
    is >> arg;

    // Split it into the individual values. This is admittedly a bit
    // complicated, but I didn't want to use Boost just for this tiny thing.
    std::size_t token_pos = 0, array_index = 0;
    while ((token_pos = arg.find(details::float_array_separator)) !=
           std::string::npos) {
        if (array_index >= SIZE) {
            throw std::invalid_argument("Too many values in the input string");
        }
        std::string token = arg.substr(0, token_pos);
        std::istringstream(token) >> values[array_index];
        arg.erase(0, token_pos + details::float_array_separator.length());
        ++array_index;
    }
    if (arg.empty() || (array_index != (SIZE - 1))) {
        throw std::invalid_argument("Too few values in the input string");
    }
    std::istringstream(arg) >> values[array_index];

    // Return the input stream.
    return is;
}

template <typename TYPE, std::size_t SIZE>
std::ostream& operator<<(std::ostream& os,
                         const value_array<TYPE, SIZE>& values) {

    for (std::size_t i = 0; i < SIZE; ++i) {
        os << values[i];
        if (i < (SIZE - 1)) {
            os << details::float_array_separator;
        }
    }
    return os;
}

template <typename TYPE, std::size_t SIZE>
value_array<TYPE, SIZE> operator*(const value_array<TYPE, SIZE>& values,
                                  TYPE factor) {

    value_array<TYPE, SIZE> result = values;
    result *= factor;
    return result;
}

template <typename TYPE, std::size_t SIZE>
value_array<TYPE, SIZE> operator/(const value_array<TYPE, SIZE>& values,
                                  TYPE factor) {

    value_array<TYPE, SIZE> result = values;
    result /= factor;
    return result;
}

}  // namespace traccc::opts
