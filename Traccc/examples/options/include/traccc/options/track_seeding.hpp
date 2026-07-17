/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"
#include "traccc/options/details/value_array.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"

namespace traccc::opts {

/// Command line options used to configure track seeding
class track_seeding : public interface,
                      public config_provider<seedfinder_config>,
                      public config_provider<seedfilter_config>,
                      public config_provider<spacepoint_grid_config>,
                      public config_provider<vector3> {

    public:
    /// Constructor
    track_seeding();

    /// @name Configuration conversion operators
    /// @{

    /// Configuration provider for the seed-finder
    explicit operator seedfinder_config() const override;
    /// Configuration provider for the seed-filter
    explicit operator seedfilter_config() const override;
    /// Configuration provider for the spacepoint grid
    explicit operator spacepoint_grid_config() const override;
    /// Configuration provider for the constant magnetic field assumed
    explicit operator vector3() const override;

    /// @}

    /// Read/process the command line options
    ///
    /// @param vm The command line options to interpret/read
    ///
    void read(const boost::program_options::variables_map& vm) override;

    /// Get a printable representation of the configuration
    std::unique_ptr<configuration_printable> as_printable() const override;

    private:
    /// Configuration for the seed-finding
    seedfinder_config m_seedfinder;
    /// Configuration for the seed filtering
    seedfilter_config m_seedfilter;

    /// Z range for the used spacepoints
    opts::value_array<float, 2> m_z_range{m_seedfinder.zMin / unit<float>::mm,
                                          m_seedfinder.zMax / unit<float>::mm};
    /// R range for the used spacepoints
    opts::value_array<float, 2> m_r_range{m_seedfinder.rMin / unit<float>::mm,
                                          m_seedfinder.rMax / unit<float>::mm};
    /// Z range for the collision region
    opts::value_array<float, 2> m_vertex_range{
        m_seedfinder.collisionRegionMin / unit<float>::mm,
        m_seedfinder.collisionRegionMax / unit<float>::mm};
    /// Delta R range for the used spacepoints
    opts::value_array<float, 2> m_delta_r_range{
        m_seedfinder.deltaRMin / unit<float>::mm,
        m_seedfinder.deltaRMax / unit<float>::mm};

};  // struct track_seeding

}  // namespace traccc::opts
