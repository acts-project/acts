/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/io/data_format.hpp"
#include "traccc/options/details/interface.hpp"

// Project include(s).
#include "traccc/definitions/common.hpp"

// System include(s).
#include <string>

namespace traccc::opts {

/// Options for the used magnetic field
struct magnetic_field : public interface {

    /// @name Options
    /// @{

    /// Use an inhomogeneous magnetic field from an input file
    bool read_from_file = false;
    /// The file containing the magnetic field description
    std::string file = "geometries/odd/odd-bfield.cvf";
    /// The data format of the magnetic field file
    traccc::data_format format = data_format::binary;
    /// Magnetic field value when not reading from a file
    float value = 2.f * unit<float>::T;

    /// @}

    /// Constructor
    magnetic_field();

    /// @name Functions implemented from @c traccc::opts::interface
    /// @{

    /// Read/process the command line options
    ///
    /// @param vm The command line options to interpret/read
    ///
    void read(const boost::program_options::variables_map& vm) override;

    /// Get a printable representation of the configuration
    std::unique_ptr<configuration_printable> as_printable() const override;

    /// @}

};  // class magnetic_field

}  // namespace traccc::opts
