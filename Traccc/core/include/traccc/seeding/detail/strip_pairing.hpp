/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/seeding/detail/strip_geometry.hpp"
#include "traccc/seeding/detail/strip_pair.hpp"

// Detray include(s).
#include <detray/geometry/identifier.hpp>

// System include(s).
#include <cstdint>

namespace traccc::details {

struct measurement_surface_before {
  template <typename measurement_t>
  TRACCC_HOST_DEVICE bool operator()(const measurement_t& measurement,
                                     const unsigned int surface_index) const {
    return measurement.surface_link().index() < surface_index;
  }
};

struct measurement_surface_before_or_equal {
  template <typename measurement_t>
  TRACCC_HOST_DEVICE bool operator()(const measurement_t& measurement,
                                     const unsigned int surface_index) const {
    return measurement.surface_link().index() <= surface_index;
  }
};

struct strip_rule_before {
  TRACCC_HOST_DEVICE bool operator()(const strip_pairing_rule& rule,
                                     const std::uint64_t surface_link) const {
    return rule.reference_surface_link < surface_link;
  }
};

template <typename measurement_collection_t>
TRACCC_HOST_DEVICE inline unsigned int lower_measurement_bound(
    const measurement_collection_t& measurements,
    const std::uint64_t surface_link) {
  // Measurement sorting compares detray barcodes by surface index, not by
  // their full encoded value. Use the same key for this binary search.
  const auto surface_index = detray::geometry::identifier{surface_link}.index();
  return collection_bound(measurements, surface_index,
                          measurement_surface_before{});
}

template <typename measurement_collection_t>
TRACCC_HOST_DEVICE inline unsigned int upper_measurement_bound(
    const measurement_collection_t& measurements,
    const std::uint64_t surface_link) {
  const auto surface_index = detray::geometry::identifier{surface_link}.index();
  return collection_bound(measurements, surface_index,
                          measurement_surface_before_or_equal{});
}

/// Rules are sorted by the full reference barcode, independently of measurement
/// sorting (which uses the surface index).
template <typename rules_t>
TRACCC_HOST_DEVICE inline unsigned int lower_strip_rule_bound(
    const rules_t& rules, const std::uint64_t surface_link) {
  return collection_bound(rules, surface_link, strip_rule_before{});
}

template <typename measurement_t>
TRACCC_HOST_DEVICE inline bool strip_pairing_value(
    const measurement_t& measurement,
    const strip_measurement_surface_info& info,
    strip_pairing_coordinate coordinate, scalar& value) {
  switch (coordinate) {
    case strip_pairing_coordinate::local0:
      value = measurement.local_position()[0];
      return true;
    case strip_pairing_coordinate::local1:
      value = measurement.local_position()[1];
      return true;
    case strip_pairing_coordinate::strip_index: {
      const int index = measured_strip_index(measurement, info.index_mapping);
      if (index < 0) {
        return false;
      }
      value = static_cast<scalar>(index);
      return true;
    }
    default:
      return false;
  }
}

template <typename measurement_t>
TRACCC_HOST_DEVICE inline bool match_strip_pair(
    const measurement_t& reference, const measurement_t& candidate,
    const strip_measurement_surface_info& reference_info,
    const strip_measurement_surface_info& candidate_info,
    const strip_pairing_rule& rule) {
  if (reference.dimensions() != 1u || candidate.dimensions() != 1u ||
      reference.surface_link().value() != rule.reference_surface_link ||
      candidate.surface_link().value() != rule.candidate_surface_link) {
    return false;
  }
  scalar reference_value{}, candidate_value{};
  if (!strip_pairing_value(reference, reference_info, rule.reference_coordinate,
                           reference_value) ||
      !strip_pairing_value(candidate, candidate_info, rule.candidate_coordinate,
                           candidate_value)) {
    return false;
  }
  switch (rule.mode) {
    case strip_pairing_mode::difference: {
      const scalar difference =
          (candidate_value - reference_value) * rule.difference_scale;
      return difference >= rule.difference_min &&
             difference <= rule.difference_max;
    }
    case strip_pairing_mode::windows:
      return reference_value >= rule.reference_min &&
             reference_value <= rule.reference_max &&
             candidate_value >= rule.candidate_min &&
             candidate_value <= rule.candidate_max;
    default:
      return false;
  }
}

/// Shared enumeration keeps the count and write passes' pair selection logic in
/// one place. Their callers must nevertheless guard against a floating-point
/// cut-boundary disagreement between the separately compiled kernels. Geometry
/// and rules must be immutable for both passes. Rules must contain exactly one
/// entry per directed surface pair; reverse pairs are not implicit.
template <typename measurements_t, typename surfaces_t, typename rules_t,
          typename visitor_t>
TRACCC_HOST_DEVICE inline void visit_strip_pairs(
    unsigned int reference_index, const measurements_t& measurements,
    const surfaces_t& surfaces, const rules_t& rules, const point3& beam_spot,
    visitor_t& visitor) {
  if (reference_index >= measurements.size()) {
    return;
  }
  const auto reference = measurements.at(reference_index);
  const auto link = reference.surface_link().value();
  const auto rule_begin = lower_strip_rule_bound(rules, link);
  if (reference.dimensions() != 1u || rule_begin == rules.size() ||
      rules.at(rule_begin).reference_surface_link != link) {
    return;
  }
  strip_measurement_surface_info reference_info{};
  if (!find_strip_surface_info(surfaces, link, reference_info) ||
      make_strip_material(reference, reference_info, beam_spot).valid == 0u) {
    return;
  }
  for (auto i = rule_begin; i < rules.size(); ++i) {
    const auto& rule = rules.at(i);
    if (rule.reference_surface_link != link) {
      break;
    }
    if (rule.category != strip_pair_category::standard &&
        rule.category != strip_pair_category::overlap) {
      continue;
    }
    strip_measurement_surface_info candidate_info{};
    if (!find_strip_surface_info(surfaces, rule.candidate_surface_link,
                                 candidate_info)) {
      continue;
    }
    const auto begin =
        lower_measurement_bound(measurements, rule.candidate_surface_link);
    const auto end =
        upper_measurement_bound(measurements, rule.candidate_surface_link);
    for (auto candidate_index = begin; candidate_index < end;
         ++candidate_index) {
      if (candidate_index == reference_index) {
        continue;
      }
      const auto candidate = measurements.at(candidate_index);
      if (match_strip_pair(reference, candidate, reference_info, candidate_info,
                           rule) &&
          make_strip_material(candidate, candidate_info, beam_spot).valid !=
              0u) {
        visitor(candidate_index, rule);
      }
    }
  }
}

}  // namespace traccc::details
