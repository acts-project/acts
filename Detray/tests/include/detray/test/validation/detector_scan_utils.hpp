// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// Detray plugin include(s)
#include "detray/plugins/svgtools/illustrator.hpp"

// Detray test include(s)
#include "detray/test/validation/svg_display.hpp"

// System include(s)
#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <ranges>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace detray::detector_scanner {

/// @brief Filter the intersection trace for overlapping surfaces and
/// remove duplicates
///
/// ACTS geometries can produce multiple overlapping portals, since some portals
/// are larger than the volumes they belong to. This leads to duplicate
/// intersections in the traces.
template <typename record_container>
inline dindex_range overlaps_removal(record_container &intersection_records,
                                     const float tol = 1e-4f *
                                                       unit<float>::mm) {
  // Document any overlaps on the range
  constexpr auto inv_idx{detray::detail::invalid_value<dindex>()};
  dindex_range overlap_idx{inv_idx, inv_idx};

  const std::size_t n_rec{intersection_records.size()};
  std::size_t n_eq_intrs{1u};  //< Number of consecutive overlapping inters.
  for (std::size_t i = 0u; i < n_rec - 1u; ++i) {
    const auto &rec = intersection_records.at(i);
    const auto &next_rec = intersection_records.at(i + 1u);

    // This record and the following are overlapping: Count until we reach
    // the end of the overlapping surfaces or end of trace
    if (math::fabs(next_rec.intersection.path() - rec.intersection.path()) <
            tol &&
        i != (n_rec - 2u)) {
      ++n_eq_intrs;
      continue;
    }

    // No overlap: continue
    if (n_eq_intrs == 1u) {
      continue;
    }
    // Found two overlapping surfaces
    if (n_eq_intrs == 2u) {
      // Two overlapping portals form a valid, connected volume boundary
      if (const auto &prev_rec = intersection_records.at(i - 1u);
          !(rec.intersection.surface().is_portal() &&
            prev_rec.intersection.surface().is_portal())) {
        auto prev_sf_desc = prev_rec.intersection.surface();
        auto sf_desc = rec.intersection.surface();

        // Other types of surfaces must not overlap!
        std::stringstream err_stream;
        err_stream << "The following surfaces overlap at\n"
                   << "POS:\n"
                   << "glob: " << prev_rec.pos
                   << ", loc: " << prev_rec.intersection.local()
                   << "\nvs.\nglob: " << rec.pos
                   << ", loc: " << rec.intersection.local() << std::endl;
        err_stream << "SURFACES:\n -> " << prev_sf_desc << std::endl;
        err_stream << " -> " << sf_desc << std::endl;

        // TODO: Fix wire_chamber geometry
        // throw std::invalid_argument(err_stream.str());
        DETRAY_ERROR_HOST(err_stream.str());
        overlap_idx = {static_cast<dindex>(i - 1u), static_cast<dindex>(i)};
      }
      // Reset to find the next range
      n_eq_intrs = 1u;
      continue;
    }

    // Found more than two overlapping surfaces: Remove oversized portal
    // intersections
    std::size_t first{(i - (n_eq_intrs - 1u))};
    std::size_t last{i};

    // Map of the volume index the portals link to and the portal
    // intersection indices in the intersection trace
    std::multimap<std::size_t, std::size_t> pt_buckets{};

    bool is_all_portals{true};
    for (std::size_t j = first; j <= last; ++j) {
      const auto &intr = intersection_records.at(j).intersection;
      if (!intr.surface().is_portal()) {
        is_all_portals = false;
      }
      pt_buckets.insert({static_cast<std::size_t>(intr.volume_link()), j});
    }

    // Remove buckets that contain only one portal: This is the valid exit
    // portal that the overlapping portals need to be matched against
    auto exit_idx{std::numeric_limits<std::size_t>::max()};
    const auto n_erased =
        std::erase_if(pt_buckets, [&pt_buckets, &exit_idx](const auto &item) {
          const auto &[vol_link, rec_idx] = item;
          if (pt_buckets.count(vol_link) == 1u) {
            exit_idx = rec_idx;
            return true;
          }
          return false;
        });

    // This is not a case of oversized portals, treat as actual overlaps
    if (n_erased != 1u || pt_buckets.empty() || !is_all_portals) {
      DETRAY_ERROR_HOST("Could not resolve exit portal in overlap correction");

      overlap_idx = {static_cast<dindex>(first), static_cast<dindex>(last)};
    } else {
      // Keep only the portal that forms the correct volume crossing
      // with the portal intersection at 'exit_idx'
      const auto &exit_rec = intersection_records.at(exit_idx);
      auto itr = intersection_records.cbegin();

      // Remove elements indiceated by ordered map with higher index first
      for (const auto &[v_link, idx] : pt_buckets | std::views::reverse) {
        // All these portals link against the exit portal
        assert(v_link == exit_rec.vol_idx);

        auto test_elem = std::next(itr, static_cast<int>(idx));
        if (exit_rec.intersection.volume_link() != test_elem->vol_idx) {
          intersection_records.erase(test_elem);
        }
      }
    }

    // Reset
    n_eq_intrs = 1u;
  }

  return overlap_idx;
}

/// Check if a set of volume indices from portal intersections from a path
/// (works even if the pairs are not sorted). That the volume links of the
/// portals at a single boundary crossing match was already checked by the
/// @c trace_intersections function at this point.
///
/// For example, with (a, b) modeling a valid portal crossing between volume 'a'
/// and 'b' (leaving out the portal indices that are also in those pairs), then:
///     (0, 1)(1, 16)(16, 17)  is valid, (0, 1)(1, 15)(16, 17) is not
///         |__|   |__|    |_                |__|   | x |   |_
///
/// Also checks that the first portal that was found, lies in the start volume.
///
/// @tparam invalid_value to keep the implementation simple, all indices are of
///                       type @c dindex here, but that does not need to be the
///                       case in the intersection code, so we need to know what
///                       the invalid value is interpreted as a @c dindex
/// @tparam check_sorted_trace if the trace is not sorted, perform a search for
///                            a fitting record across the trace.
/// @tparam entry_type the record entry, which must contain the portal index
///                    and the volume index the portal was discovered in.
///
/// @param trace the recorded portal crossings between volumes
/// @param start_volume where the ray started
///
/// @return true if the volumes indices form a connected chain.
template <dindex invalid_value = dindex_invalid, bool check_sorted_trace = true,
          typename entry_type = std::pair<dindex, dindex>>
inline bool check_connectivity(
    std::vector<std::pair<entry_type, entry_type>> trace,
    dindex start_volume = 0u) {
  /// Error messages
  std::stringstream err_stream;

  /// Print errors of this function
  auto print_err = [](const std::stringstream &stream) {
    std::cerr << "\n<<<<<<<<<<<<<<< ERROR in connectivity check\n" << std::endl;
    std::cerr << stream.str() << std::endl;
    std::cerr << "\n>>>>>>>>>>>>>>>\n" << std::endl;
  };

  // There must always be portals!
  if (trace.empty()) {
    err_stream << "Trace empty!";
    print_err(err_stream);

    return false;
  }
  // Keep record of leftovers
  std::stringstream record_stream;

  // Where are we on the trace?
  dindex current_volume = start_volume;

  // If the intersection trace comes from the ray gun/trace intersections
  // function it should be sorted, which is the stronger constraint
  using vector_t = decltype(trace);
  using records_iterator_t = typename vector_t::iterator;
  using index_t = std::iter_difference_t<records_iterator_t>;
  std::function<records_iterator_t(index_t)> get_connected_record;
  if constexpr (check_sorted_trace) {
    // Get the next record
    get_connected_record =
        [&trace, &current_volume](index_t next) -> records_iterator_t {
      // Make sure that the record contains the volume that is currently
      // being checked for connectivity
      if (auto rec = trace.begin() + next;
          rec != trace.end() &&
          ((std::get<1>(rec->first) == current_volume) ||
           (std::get<1>(rec->second) == current_volume))) {
        return rec;
      }
      return trace.end();
    };
  } else {
    // Search for the existence of a fitting record over the entire trace
    get_connected_record =
        [&trace, &current_volume](index_t /*next*/) -> records_iterator_t {
      return std::ranges::find_if(
          trace, [&](const std::pair<entry_type, entry_type> &rec) -> bool {
            return (std::get<1>(rec.first) == current_volume) ||
                   (std::get<1>(rec.second) == current_volume);
          });
    };
  }

  // Init chain search
  index_t i{0};
  auto record = get_connected_record(i);

  // Check first volume index, which has no partner otherwise
  if (std::get<1>(record->first) != start_volume) {
    err_stream << "First record does not start at given initial volume: "
               << std::get<1>(record->first) << " vs. " << start_volume;

    print_err(err_stream);

    return false;
  }

  // Walk along the trace as long as a connection is found
  while (record != trace.end()) {
    auto first_vol = std::get<1>(record->first);
    auto second_vol = std::get<1>(record->second);

    record_stream << "On volume: " << current_volume << " and record ("
                  << first_vol << ", " << second_vol << ")";

    // update to next volume
    current_volume = (current_volume == first_vol ? second_vol : first_vol);

    record_stream << " -> next volume: " << current_volume << std::endl;

    // Don't search this key again -> only one potential key with current
    // index left
    if constexpr (!check_sorted_trace) {
      trace.erase(record);
    }

    // find connected record for the current volume
    record = get_connected_record(++i);
  }

  // There are unconnected elements left (we didn't leave world before
  // termination)
  if (current_volume != invalid_value) {
    err_stream << "Didn't leave world or unconnected elements left in trace:"
               << "\n\nValid connections that were found:" << std::endl;
    err_stream << record_stream.str();

    err_stream << "\nPairs left to match:" << std::endl;
    for (auto j = static_cast<std::size_t>(i); j < trace.size(); ++j) {
      auto first_vol = std::get<1>(trace[j].first);
      auto second_vol = std::get<1>(trace[j].second);

      err_stream << "(" << first_vol << ", " << second_vol << ")" << std::endl;
    }

    print_err(err_stream);

    return false;
  }

  return true;
}

/// Check if a recording of portal/module intersections form a coherent
/// trace through a geometry.
///
/// Various consistency checks are performed on the trace:
///     - Was a surface intersected as part of the wrong volume?
///     - Do the volume links at adjacent portals point at each other?
///     - Is a portal crossing happening within a volume, as opposed to at its
///       boundary surfaces?
///     - Do portals always appear in pairs (exit + entry portal)?
///
/// @tparam record_container contains volume indices and intersections
///
/// @param volume_record the recorded portal crossings between volumes
/// @param start_volume where the ray started
///
/// @note the input record needs to be sorted according to the distance from the
///       ray origin
///
/// @return a set of volume connections that were found by portal intersection
///         of a ray.
template <dindex invalid_value = dindex_invalid, typename record_container>
inline auto trace_intersections(const record_container &intersection_records,
                                dindex start_volume = 0u) {
  /// surface index and index of the volume the intersection was found in
  using trace_entry = std::pair<dindex, dindex>;
  /// Pairs of adjacent portals along the ray
  std::vector<std::pair<trace_entry, trace_entry>> portal_trace = {};
  /// Trace of module surfaces
  std::vector<trace_entry> module_trace = {};
  /// Debug output if an error in the trace is discovered
  std::stringstream record_stream;
  /// Error messages
  std::stringstream err_stream;
  bool error_code{true};

  /// Readable access to the data of a recorded intersection
  struct record {
    const typename record_container::value_type &entry;

    /// getter
    /// @{
    inline bool is_invalid() const {
      return entry.intersection.surface().identifier().is_invalid();
    }
    inline auto surface_idx() const {
      return entry.intersection.surface().index();
    }
    inline auto surface_volume_idx() const {
      return entry.intersection.surface().volume();
    }
    inline auto &inters() const { return entry.intersection; }
    inline auto volume_idx() const { return entry.vol_idx; }
    inline auto volume_link() const { return entry.intersection.volume_link(); }
    inline auto dist() const { return entry.intersection.path(); }
    inline bool is_portal() const {
      return entry.intersection.surface().is_portal();
    }
    inline bool is_sensitive() const {
      return entry.intersection.surface().is_sensitive();
    }
    inline bool is_passive() const {
      return entry.intersection.surface().is_passive();
    }
    /// @}
  };

  /// Print errors of this function
  auto print_err = [&error_code](const std::stringstream &stream) {
    std::cerr << "\n<<<<<<<<<<<<<<< ERROR intersection trace\n" << std::endl;
    std::cerr << stream.str() << std::endl;
    std::cerr << "\n>>>>>>>>>>>>>>>\n" << std::endl;
    error_code = false;
  };

  // No intersections found by ray
  if (intersection_records.empty()) {
    err_stream << "No surfaces found in detector!";
    print_err(err_stream);

    return std::make_tuple(portal_trace, module_trace, error_code);
  }

  // If there is only one surface in the trace, it must be a portal
  if (intersection_records.size() == 1u) {
    const record rec{intersection_records.at(0u)};

    // No exit potal
    if (!rec.is_portal()) {
      if (rec.is_invalid()) {
        err_stream << "No surfaces found in detector!";
        print_err(err_stream);
      } else {
        const std::string sf_type{rec.is_sensitive() ? "sensitive" : "passive"};

        err_stream << "We don't leave the detector by portal!" << std::endl;
        err_stream << "Only found single " << sf_type
                   << " surface: portal(s) missing!";

        print_err(err_stream);
      }

      return std::make_tuple(portal_trace, module_trace, error_code);
    }
  }

  // Go through recorded intersection (two at a time)
  dindex current_vol = start_volume;
  // The first entry is the dummy record that preserves initial track
  // parameters, skip it
  const std::size_t start_idx{
      record{intersection_records.at(0)}.is_invalid() ? 1u : 0u};
  for (std::size_t rec = start_idx; rec < (intersection_records.size() - 1u);) {
    const record current_rec = record{intersection_records.at(rec)};
    const record next_rec = record{intersection_records.at(rec + 1u)};

    // Keep a debug stream
    record_stream << current_rec.volume_idx() << "\t" << current_rec.inters()
                  << std::endl;

    // If the current record is not a portal, add an entry to the module
    // trace and continue in more fine-grained steps (sensitive/passive
    // surfaces do not come in pairs)
    if (!current_rec.is_portal()) {
      // Check that the surface was found in the volume it claims to
      // belong to
      const bool is_in_volume =
          (current_rec.volume_idx() == current_rec.surface_volume_idx()) &&
          (current_rec.surface_volume_idx() == current_vol);
      if (is_in_volume) {
        module_trace.emplace_back(current_rec.surface_idx(),
                                  current_rec.volume_idx());
      } else {
        err_stream << "\n(!!) Surface " << current_rec.surface_idx()
                   << " outside of its volume (Found in: " << current_vol
                   << ", belongs in: " << current_rec.surface_volume_idx()
                   << ")\n";

        err_stream << record_stream.str();

        print_err(err_stream);

        return std::make_tuple(portal_trace, module_trace, error_code);
      }
      ++rec;
      continue;
    }
    // If the record is a portal, the current volume switches
    // the portals in the pair may be unordered, so check both
    else if (current_vol == current_rec.volume_idx()) {
      current_vol = current_rec.volume_link();
    } else if (current_vol == next_rec.volume_idx()) {
      current_vol = next_rec.volume_link();
    }

    record_stream << next_rec.volume_idx() << "\t" << next_rec.inters()
                  << std::endl;

    // Check that also the second surface was found in the volume it claims
    // to belong to
    const bool is_in_volume =
        next_rec.volume_idx() == next_rec.surface_volume_idx();
    // Is this doublet connected via a valid portal intersection?
    const bool is_valid =
        (current_rec.inters() == next_rec.inters()) &&
        (current_rec.volume_idx() == next_rec.volume_link()) &&
        (next_rec.volume_idx() == current_rec.volume_link());
    // Is this indeed a portal crossing, i.e. changing volumes?
    const bool is_self_link = current_rec.volume_idx() == next_rec.volume_idx();
    // Is the record doublet we picked made up of a portal and a module?
    const bool is_mixed = (current_rec.is_portal() && !next_rec.is_portal()) ||
                          (next_rec.is_portal() && !current_rec.is_portal());

    if (!is_in_volume) {
      record_stream << "\n(!!) Surface outside of its volume (Found: "
                    << next_rec.volume_idx()
                    << ", belongs in: " << next_rec.surface_volume_idx() << ")"
                    << std::endl;
    }
    if (!is_valid) {
      record_stream << "\n(!!) Not a valid portal crossing ("
                    << current_rec.volume_idx() << " <-> "
                    << next_rec.volume_idx() << "):\nPortals are not "
                    << "connected, either geometrically or by linking!"
                    << std::endl;
    }
    if (is_self_link) {
      record_stream << "\n(!!) Found portal crossing inside volume ("
                    << current_rec.volume_idx() << ")!" << std::endl;
    }
    if (is_mixed) {
      record_stream << "\n(!!) Portal crossing involves module surface ("
                    << current_rec.volume_idx() << " <-> "
                    << next_rec.volume_idx() << ")! The second surface in"
                    << " this portal crossing is not a portal!" << std::endl;
    }
    if (is_in_volume && is_valid && !is_self_link && !is_mixed) {
      // Insert into portal trace
      trace_entry lower{current_rec.surface_idx(), current_rec.volume_idx()};
      trace_entry upper{next_rec.surface_idx(), next_rec.volume_idx()};
      portal_trace.emplace_back(lower, upper);
    }
    // Something went wrong
    else {
      // Print search log
      err_stream << "\nError in portal matching:\n" << std::endl;
      err_stream << "volume id\t(intersection info)" << std::endl;

      err_stream << record_stream.str() << std::endl;

      err_stream << "-----\nINFO: Ray terminated at portal x-ing "
                 << (rec + 1) / 2 << ":\n"
                 << current_rec.inters() << " <-> " << next_rec.inters()
                 << std::endl;

      const record rec_front{record{intersection_records.front()}.is_invalid()
                                 ? intersection_records[1u]
                                 : intersection_records.front()};
      const record rec_back{intersection_records.back()};
      err_stream << "Start volume : " << start_volume << std::endl;
      err_stream << "- first recorded intersection: (sf id:"
                 << rec_front.surface_idx() << ", dist:" << rec_front.dist()
                 << ")," << std::endl;
      err_stream << "- last recorded intersection:  (sf id:"
                 << rec_back.surface_idx() << ", dist:" << rec_back.dist()
                 << "),";

      print_err(err_stream);

      return std::make_tuple(portal_trace, module_trace, error_code);
    }

    // Advance to inspect next pair
    rec += 2u;
  }

  // Look at the last entry, which is a single portal
  if (const record rec_back{intersection_records.back()};
      !rec_back.is_portal()) {
    err_stream << "We don't leave the detector by portal!";
    print_err(err_stream);
  } else {
    trace_entry lower(rec_back.surface_idx(), rec_back.volume_idx());
    trace_entry upper(rec_back.surface_idx(), rec_back.volume_link());
    portal_trace.emplace_back(lower, upper);
  }

  return std::make_tuple(portal_trace, module_trace, error_code);
}

/// Build an adjacency list from intersection traces.
///
/// @tparam portal_trace_type container of portal link pairs
/// @tparam module_trace_type container of module surface links
///
/// @param portal_trace the portal indices and their volume links (in adjacent
///                     portal pairs)
/// @param module_trace the module indices and their volume links
/// @param obj_hashes record which modules/portals were already added
///
/// @return an adjacency list from the traced ray scan of a given geometry.
template <dindex invalid_value = dindex_invalid, typename portal_trace_type,
          typename module_trace_type,
          typename entry_type = std::pair<dindex, dindex>>
  requires std::is_same_v<typename portal_trace_type::value_type,
                          std::pair<entry_type, entry_type>> &&
           std::is_same_v<typename module_trace_type::value_type, entry_type>
inline auto build_adjacency(
    const portal_trace_type &portal_trace,
    const module_trace_type &module_trace,
    std::map<dindex, std::map<dindex, dindex>> &adj_list,
    std::unordered_set<dindex> &obj_hashes) {
  // Every module that was recorded adds a link to the mother volume
  for (const auto &record : module_trace) {
    const auto sf_index = std::get<0>(record);
    const auto vol_index = std::get<1>(record);
    // Check whether we have seen this module in this volume before
    if (obj_hashes.find(sf_index) == obj_hashes.end()) {
      adj_list[vol_index][vol_index]++;
      obj_hashes.insert(sf_index);
    }
  }

  // Portal in first volume links to second volume in the record
  for (const auto &record : portal_trace) {
    const auto pt_index_1 = std::get<0>(record.first);
    const auto vol_index_1 = std::get<1>(record.first);
    const auto pt_index_2 = std::get<0>(record.second);
    const auto vol_index_2 = std::get<1>(record.second);

    if (obj_hashes.find(pt_index_1) == obj_hashes.end()) {
      adj_list[vol_index_1][vol_index_2]++;
      obj_hashes.insert(pt_index_1);
    }
    // Assume the return link for now (filter out portal that leaves world)
    if (vol_index_2 != invalid_value &&
        obj_hashes.find(pt_index_2) == obj_hashes.end()) {
      adj_list[vol_index_2][vol_index_1]++;
      obj_hashes.insert(pt_index_2);
    }
  }

  return adj_list;
}

/// Build an adjacency list from intersection traces.
///
/// @tparam portal_trace_type container of portal link pairs
/// @tparam module_trace_type container of module links
///
/// @param portal_trace the portal indices and their volume links (in adjacent
///                     portal pairs)
/// @param module_trace the module indices and their volume links
/// @param obj_hashes record which modules/portals were already added
///
/// @return an adjacency list from the traced ray scan of a given geometry.
template <dindex invalid_value = dindex_invalid, typename portal_trace_type,
          typename module_trace_type,
          typename entry_type = std::pair<dindex, dindex>>
  requires std::is_same_v<typename portal_trace_type::value_type,
                          std::pair<entry_type, entry_type>> &&
           std::is_same_v<typename module_trace_type::value_type, entry_type>
inline auto build_adjacency(const portal_trace_type &portal_trace,
                            const module_trace_type &module_trace,
                            dvector<dindex> &adj_matrix,
                            std::unordered_set<dindex> &obj_hashes) {
  const auto dim = static_cast<dindex>(math::sqrt(adj_matrix.size()));

  // Every module that was recorded adds a link to the mother volume
  for (const auto &record : module_trace) {
    const auto sf_index = std::get<0>(record);
    const auto vol_index = std::get<1>(record);
    // Check whether we have seen this module in this volume before
    if (obj_hashes.find(sf_index) == obj_hashes.end()) {
      adj_matrix[dim * vol_index + vol_index]++;
      obj_hashes.insert(sf_index);
    }
  }

  // Portal in first volume links to second volume in the record
  for (const auto &record : portal_trace) {
    const auto pt_index_1 = std::get<0>(record.first);
    const auto vol_index_1 = std::get<1>(record.first);
    const auto pt_index_2 = std::get<0>(record.second);
    const auto vol_index_2 = std::get<1>(record.second);

    if (obj_hashes.find(pt_index_1) == obj_hashes.end()) {
      dindex mat_elem_vol1{dindex_invalid};
      // Assume the return link for now (filtering out portals that leave
      // world)
      if (vol_index_2 != invalid_value) {
        mat_elem_vol1 = dim * vol_index_1 + vol_index_2;

        if (obj_hashes.find(pt_index_2) == obj_hashes.end()) {
          adj_matrix[dim * vol_index_2 + vol_index_1]++;
          obj_hashes.insert(pt_index_2);
        }
      } else {
        mat_elem_vol1 = dim * vol_index_1 + dim - 1;
      }
      adj_matrix[mat_elem_vol1]++;
      obj_hashes.insert(pt_index_1);
    }
  }

  return adj_matrix;
}

/// Run all checks on an intersection trace.
///
/// @param[in] intersection_trace the intersection records along the track
/// @param[in] start_index the index of the intended start volume
/// @param[out] adj_mat_scan adjacency matrix to be filled for the detector
/// @param[out] obj_hashes objects in a volume that were already visisted
///
/// @return true if the checks were successful
template <typename detector_t, typename record_t>
inline bool check_trace(const std::vector<record_t> &intersection_trace,
                        const dindex start_index, dvector<dindex> &adj_mat_scan,
                        std::unordered_set<dindex> &obj_hashes) {
  using nav_link_t = typename detector_t::surface_type::navigation_link;
  static constexpr auto leaving_world{
      detray::detail::invalid_value<nav_link_t>()};

  // Create a trace of the volume indices that were encountered
  // and check that portal intersections are connected
  auto [portal_trace, surface_trace, err_code] =
      detector_scanner::trace_intersections<leaving_world>(intersection_trace,
                                                           start_index);

  // Is the succession of volumes consistent ?
  err_code = err_code &&
             detector_scanner::check_connectivity<leaving_world>(portal_trace);

  if (!adj_mat_scan.empty()) {
    // Build an adjacency matrix from this trace that can be checked
    // against the geometry hash (see 'track_geometry_changes')
    detector_scanner::build_adjacency<leaving_world>(
        portal_trace, surface_trace, adj_mat_scan, obj_hashes);
  }

  return err_code;
}

/// Print the failed intersection trace as svg (dumped to files)
///
/// @param gctx current geometry context
/// @param det the detector object
/// @param vol_names the volume name map of the detector
/// @param test_name the name of the test for which to print the error
/// @param test_track trajectory that was used for the scan (ray or helix)
/// @param truth_trace the intersection records along the test track
/// @param svg_style svgtools style for the detector display
/// @param i_track index of the test track
/// @param n_track total number of test tracks
template <typename detector_t, typename trajectory_t, typename truth_trace_t,
          typename recorded_trace_t>
inline void display_error(
    const typename detector_t::geometry_context gctx, const detector_t &det,
    const typename detector_t::name_map &vol_names,
    const std::string &test_name, const trajectory_t &test_track,
    const truth_trace_t &truth_trace,
    const detray::svgtools::styling::style &svg_style,
    const std::size_t i_track, [[maybe_unused]] const std::size_t n_tracks,
    const recorded_trace_t &recorded_trace = {},
    const dindex_range overlap_idx = {detray::detail::invalid_value<dindex>(),
                                      detray::detail::invalid_value<dindex>()},
    const bool verbose = true) {
  // Creating the svg generator for the detector.
  detray::svgtools::illustrator il{det, vol_names, svg_style};
  il.show_info(true);
  il.hide_eta_lines(true);
  il.hide_portals(false);
  il.hide_passives(false);

  std::string track_type{};
  static constexpr auto is_ray{
      std::is_same_v<trajectory_t,
                     detray::detail::ray<typename detector_t::algebra_type>>};
  if constexpr (is_ray) {
    track_type = "ray";
  } else {
    track_type = "helix";
  }

  if (verbose) {
    DETRAY_ERROR_HOST("\nFailed on " << track_type << ": " << i_track << "/"
                                     << n_tracks << "\n"
                                     << test_track);
  }

  detray::detail::svg_display(gctx, il, truth_trace, test_track,
                              track_type + "_" + std::to_string(i_track),
                              test_name, recorded_trace, overlap_idx, verbose);
}

/// Print an intersection trace
template <typename truth_trace_t>
inline std::string print_trace(const truth_trace_t &truth_trace,
                               std::size_t n) {
  std::stringstream out_stream{};
  out_stream << "TRACE NO. " << n << std::endl;

  for (const auto &[idx, record] : detray::views::enumerate(truth_trace)) {
    out_stream << "\nRecord " << idx << std::endl;

    out_stream << " -> volume " << record.vol_idx << std::endl;

    const auto pos = record.pos;
    const auto dir = record.dir;
    out_stream << " -> track pos: [" << pos[0] << ", " << pos[1] << ", "
               << pos[2] << std::endl;
    out_stream << " -> track dir: [" << dir[0] << ", " << dir[1] << ", "
               << dir[2] << std::endl;

    out_stream << " -> intersection " << record.intersection << std::endl;
  }

  return out_stream.str();
}

/// Print an adjacency list
inline std::string print_adj(const dvector<dindex> &adjacency_matrix) {
  std::size_t dim = static_cast<dindex>(math::sqrt(adjacency_matrix.size()));
  std::stringstream out_stream{};

  for (std::size_t i = 0u; i < dim - 1; ++i) {
    out_stream << "[>>] Node with index " << i << std::endl;
    out_stream << " -> edges: " << std::endl;
    for (std::size_t j = 0u; j < dim; ++j) {
      const auto degr = adjacency_matrix[dim * i + j];
      if (degr == 0) {
        continue;
      }
      std::string n_occur =
          degr > 1 ? "\t\t\t\t(" + std::to_string(degr) + "x)" : "";

      // Edge that leads out of the detector world
      if (j == dim - 1 && degr != 0) {
        out_stream << "    -> leaving world " + n_occur << std::endl;
      } else {
        out_stream << "    -> " << std::to_string(j) + "\t" + n_occur
                   << std::endl;
      }
    }
  }

  return out_stream.str();
}

}  // namespace detray::detector_scanner
