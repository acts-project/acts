/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/utils/event_data.hpp"

#include "traccc/edm/measurement_collection.hpp"
#include "traccc/io/csv/make_cell_reader.hpp"
#include "traccc/io/csv/make_hit_reader.hpp"
#include "traccc/io/csv/make_measurement_edm.hpp"
#include "traccc/io/csv/make_measurement_hit_id_reader.hpp"
#include "traccc/io/csv/make_measurement_reader.hpp"
#include "traccc/io/csv/make_particle_reader.hpp"
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/utils/particle.hpp"

// System include(s).
#include <algorithm>
#include <filesystem>

namespace traccc {

event_data::event_data(const std::string& event_dir, const std::size_t event_id,
                       vecmem::memory_resource& resource,
                       bool use_acts_geom_source, const host_detector* det,
                       data_format format, bool include_silicon_cells)
    : m_measurements(resource),
      m_event_dir(event_dir),
      m_event_id(event_id),
      m_mr(resource) {

    // Currently, we only support csv type for event data
    assert(format == data_format::csv);
    if (format == data_format::csv) {
        setup_csv(use_acts_geom_source, det, include_silicon_cells);
    }
}

void event_data::setup_csv(bool use_acts_geom_source, const host_detector* det,
                           bool include_silicon_cells) {

    /********************
     *  Read Csv files  *
     ********************/

    // CSV IO EDM containers
    std::vector<io::csv::cell> csv_cells;
    std::vector<io::csv::hit> csv_hits;
    std::vector<io::csv::measurement> csv_measurements;
    std::vector<io::csv::measurement_hit_id> csv_meas_hit_ids;
    std::vector<io::csv::particle> csv_particles;

    if (include_silicon_cells) {
        // Read the cells from the relevant event file
        std::string io_cells_file =
            io::get_absolute_path((std::filesystem::path(m_event_dir) /
                                   std::filesystem::path(io::get_event_filename(
                                       m_event_id, "-cells.csv")))
                                      .native());

        auto creader = io::csv::make_cell_reader(io_cells_file);
        io::csv::cell iocell;

        while (creader.read(iocell)) {
            csv_cells.push_back(iocell);
        }
    }

    // Read the hits from the relevant event file
    std::string io_hits_file = io::get_absolute_path(
        (std::filesystem::path(m_event_dir) /
         std::filesystem::path(io::get_event_filename(m_event_id, "-hits.csv")))
            .native());

    auto hreader = io::csv::make_hit_reader(io_hits_file);
    {
        io::csv::hit iohit;
        while (hreader.read(iohit)) {
            csv_hits.push_back(iohit);
        }
    }

    // Read the measurements from the relevant event file
    std::string io_measurements_file =
        io::get_absolute_path((std::filesystem::path(m_event_dir) /
                               std::filesystem::path(io::get_event_filename(
                                   m_event_id, "-measurements.csv")))
                                  .native());

    auto mreader = io::csv::make_measurement_reader(io_measurements_file);
    {
        io::csv::measurement iomeas;
        while (mreader.read(iomeas)) {
            csv_measurements.push_back(iomeas);
        }
    }

    // Read the measurement hit id from the relevant event file
    std::string io_measurement_hit_id_file =
        io::get_absolute_path((std::filesystem::path(m_event_dir) /
                               std::filesystem::path(io::get_event_filename(
                                   m_event_id, "-measurement-simhit-map.csv")))
                                  .native());

    auto mhid_reader =
        io::csv::make_measurement_hit_id_reader(io_measurement_hit_id_file);
    {
        io::csv::measurement_hit_id mh_id;
        while (mhid_reader.read(mh_id)) {
            csv_meas_hit_ids.push_back(mh_id);
        }
    }

    // Read the particles from the relevant event file
    std::string io_particles_file =
        io::get_absolute_path((std::filesystem::path(m_event_dir) /
                               std::filesystem::path(io::get_event_filename(
                                   m_event_id, "-particles_initial.csv")))
                                  .native());

    auto preader = io::csv::make_particle_reader(io_particles_file);
    {
        io::csv::particle ioptc;
        while (preader.read(ioptc)) {
            csv_particles.push_back(ioptc);
        }
    }

    /********************
     * Make geom_id map *
     ********************/

    // For Acts data, build a map of acts->detray geometry IDs
    std::map<geometry_id, geometry_id> acts_to_detray_id;
    if (use_acts_geom_source) {
        host_detector_visitor<detector_type_list>(
            *det, [&acts_to_detray_id]<typename detector_traits_t>(
                      const typename detector_traits_t::host& d) {
                for (const auto& surface_desc : d.surfaces()) {
                    acts_to_detray_id[surface_desc.source] =
                        surface_desc.identifier().value();
                }
            });
    }

    /********************
     * Make EDM maps    *
     ********************/

    bool csv_measurements_have_time = false;

    for (const auto& iomeas : csv_measurements) {
        if (std::fabs(iomeas.time) != 0.f) {
            csv_measurements_have_time = true;
            break;
        }
    }

    // TODO: Put some proper logging here in a future PR
    if (csv_measurements_have_time) {
        std::cout << "Using measurement time" << std::endl;
    } else {
        std::cout << "Using hit time" << std::endl;
    }

    // Measurement map
    for (const auto& iomeas : csv_measurements) {
        // Construct the measurement object.
        m_measurements.resize(m_measurements.size() + 1u);
        auto meas = m_measurements.at(m_measurements.size() - 1u);
        if (use_acts_geom_source) {
            traccc::io::csv::make_measurement_edm(iomeas, meas,
                                                  &acts_to_detray_id);
        } else {
            traccc::io::csv::make_measurement_edm(iomeas, meas, nullptr);
        }

        if (!csv_measurements_have_time) {
            const auto hid = csv_meas_hit_ids.at(iomeas.measurement_id).hit_id;
            const auto& iohit = csv_hits.at(hid);

            meas.time() = iohit.tt;
        }

        if (iomeas.measurement_id <=
            std::numeric_limits<measurement_id_type>::max()) {
            m_measurement_map[static_cast<measurement_id_type>(
                iomeas.measurement_id)] = meas;
        } else {
            throw std::runtime_error("Measurement ID exceeds the bound");
        }
    }

    // Particle map
    for (const auto& ioptc : csv_particles) {

        point3 pos{ioptc.vx, ioptc.vy, ioptc.vz};
        vector3 mom{ioptc.px, ioptc.py, ioptc.pz};

        m_particle_map[ioptc.particle_id] =
            traccc::particle{ioptc.particle_id, ioptc.particle_type,
                             ioptc.process,     pos,
                             ioptc.vt,          mom,
                             ioptc.m,           ioptc.q};
    }

    /************************************
     * Make measurement to particle map *
     ************************************/

    // When including silicon cells
    if (include_silicon_cells) {

        std::map<measurement_proxy, std::vector<io::csv::cell>>
            meas_to_cluster_map;

        for (const auto& iocell : csv_cells) {

            // Fill the measurement_to_cluster_map
            auto meas_id = iocell.measurement_id;
            auto hid = csv_meas_hit_ids[meas_id].hit_id;
            const auto& iohit = csv_hits[hid];

            measurement_proxy meas =
                m_measurement_map.at(static_cast<measurement_id_type>(meas_id));
            meas_to_cluster_map[meas].push_back(iocell);

            const auto& ptc = m_particle_map.at(iohit.particle_id);
            m_cell_to_particle_map[iocell] = ptc;
        }

        // Fill the meas_to_particle_map
        for (auto const& [ms, cluster] : meas_to_cluster_map) {
            for (const auto& cell : cluster) {
                const auto& ptc = m_cell_to_particle_map.at(cell);
                m_meas_to_ptc_map[ms][ptc]++;
            }
        }
    }

    for (const auto& iomeas : csv_measurements) {

        // Hit index
        const auto hid = csv_meas_hit_ids.at(iomeas.measurement_id).hit_id;

        // Make spacepoint
        const auto& iohit = csv_hits.at(hid);
        point3 global_pos{iohit.tx, iohit.ty, iohit.tz};
        point3 global_mom{iohit.tpx, iohit.tpy, iohit.tpz};

        // Make particle
        const auto& ptc = m_particle_map.at(iohit.particle_id);

        // Construct the measurement object.
        measurement_proxy meas = m_measurement_map.at(
            static_cast<measurement_id_type>(iomeas.measurement_id));

        // Fill measurement to truth global position and momentum map
        m_meas_to_param_map[meas] = std::make_pair(global_pos, global_mom);

        // Fill particle to measurement map
        auto& meas_vec = m_ptc_to_meas_map[ptc];

        meas_vec.insert(std::upper_bound(meas_vec.begin(), meas_vec.end(), meas,
                                         [](const measurement_proxy& val,
                                            const measurement_proxy& old) {
                                             return old.time() > val.time();
                                         }),
                        meas);

        if (!include_silicon_cells) {
            auto insert_return = m_meas_to_ptc_map.insert({meas, {}});
            if (insert_return.second == false) {
                // TODO: Put some logging here when that's ready
            } else {
                (*(insert_return.first)).second[ptc] = 1u;
            }
        }
    }
}

void event_data::fill_cca_result(
    const edm::silicon_cell_collection::host& cells,
    const edm::silicon_cluster_collection::host& cca_clusters,
    const edm::measurement_collection::host& cca_measurements,
    const detector_conditions_description::host& det_cond) {

    const std::size_t n_cca_clusters = cca_measurements.size();

    std::map<measurement_proxy, std::vector<io::csv::cell>>
        found_meas_to_cluster_map;

    for (std::size_t i = 0; i < n_cca_clusters; i++) {
        const auto meas = cca_measurements.at(i);
        const auto cluster = cca_clusters[i];

        std::vector<io::csv::cell> iocells;
        for (const unsigned int cell_idx : cluster.cell_indices()) {

            const auto cell = cells.at(cell_idx);
            io::csv::cell iocell{
                det_cond.acts_geometry_id().at(cell.module_index()),
                0u,
                cell.channel0(),
                cell.channel1(),
                static_cast<float>(cell.time()),
                static_cast<float>(cell.activation())};

            iocells.push_back(iocell);
        }
        found_meas_to_cluster_map[meas] = iocells;
    }

    // Construct a mapping of geometry IDs to cells and the associated
    // particles, in order to get much faster lookups later.
    std::map<geometry_id, std::vector<std::reference_wrapper<const std::decay_t<
                              decltype(m_cell_to_particle_map)>::value_type>>>
        geo_id_to_cell_to_particle_map;

    for (auto const& v : m_cell_to_particle_map) {
        geo_id_to_cell_to_particle_map[v.first.geometry_id].emplace_back(v);
    }

    for (auto const& [ms, cluster] : found_meas_to_cluster_map) {

        std::map<uint64_t, std::size_t> meas_counts;

        // Cells from CCL
        for (const auto& cell1 : cluster) {

            // Cells from truth [cell, particle] map
            for (auto const& v :
                 geo_id_to_cell_to_particle_map[cell1.geometry_id]) {
                const auto& cell2 = v.get().first;
                const auto& ptc = v.get().second;
                assert(cell1.geometry_id == cell2.geometry_id);
                // Increase the particle number if the cell is the same
                if (cell1.channel0 == cell2.channel0 &&
                    cell1.channel1 == cell2.channel1) {
                    m_found_meas_to_ptc_map[ms][ptc]++;
                    meas_counts[cell2.measurement_id]++;
                }
            }
        }

        // Find most contributing measurement and its corresponding hit
        using pair_type = decltype(meas_counts)::value_type;
        auto pr =
            std::max_element(std::begin(meas_counts), std::end(meas_counts),
                             [](const pair_type& p1, const pair_type& p2) {
                                 return p1.second < p2.second;
                             });

        const auto& meas_id = pr->first;

        if (meas_id <= std::numeric_limits<measurement_id_type>::max()) {
            m_found_meas_to_param_map[ms] = m_meas_to_param_map
                [m_measurement_map[static_cast<measurement_id_type>(meas_id)]];

        } else {
            throw std::runtime_error("Measurement ID exceeds the bound");
        }
    }
}
}  // namespace traccc
