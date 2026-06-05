/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/data_format.hpp"

// System include(s).
#include <iostream>

namespace traccc {

std::ostream& operator<<(std::ostream& out, data_format format) {

    switch (format) {
        case data_format::csv:
            out << "csv";
            break;
        case data_format::binary:
            out << "binary";
            break;
        case data_format::json:
            out << "json";
            break;
        case data_format::obj:
            out << "wavefront obj";
            break;
        default:
            out << "?!?unknown?!?";
            break;
    }
    return out;
}

}  // namespace traccc
