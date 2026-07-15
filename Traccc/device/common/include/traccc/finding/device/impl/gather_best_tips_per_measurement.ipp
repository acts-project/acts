/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/array_insertion_mutex.hpp"

// Project include(s).
#include "traccc/utils/prob.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

// System include(s).
#include <cassert>
#include <limits>

namespace traccc::device {

template <typename algebra_t, concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void gather_best_tips_per_measurement(
    global_index_t thread_id, const barrier_t& barrier,
    const gather_best_tips_per_measurement_payload<algebra_t>& payload) {

    const vecmem::device_vector<const unsigned int> tips(payload.tips);
    const vecmem::device_vector<const candidate_link> links(payload.links);
    const edm::measurement_collection::const_device measurements(
        payload.measurements);
    vecmem::device_vector<unsigned long long int> insertion_mutex(
        payload.insertion_mutex);
    vecmem::device_vector<unsigned int> tip_index(payload.tip_index);
    vecmem::device_vector<typename algebra_t::scalar> tip_pval(
        payload.tip_pval);
    const unsigned int n_meas = measurements.size();

    typename algebra_t::scalar pval = 0.f;
    unsigned int link_idx = 0;
    unsigned int num_states = 0;

    bool need_to_write = true;
    candidate_link L{
        .step = std::numeric_limits<unsigned int>::max(),
        .previous_candidate_idx = std::numeric_limits<unsigned int>::max(),
        .meas_idx = std::numeric_limits<unsigned int>::max(),
        .seed_idx = std::numeric_limits<unsigned int>::max(),
        .n_skipped = std::numeric_limits<unsigned int>::max(),
        .n_consecutive_skipped = std::numeric_limits<unsigned int>::max(),
        .chi2 = std::numeric_limits<traccc::scalar>::max(),
        .chi2_sum = std::numeric_limits<traccc::scalar>::max(),
        .ndf_sum = std::numeric_limits<unsigned int>::max()};

    if (thread_id < tips.size()) {
        link_idx = tips.at(thread_id);
        const candidate_link link = links.at(link_idx);
        pval =
            prob(link.chi2_sum,
                 static_cast<typename algebra_t::scalar>(link.ndf_sum) - 5.f);
        num_states = link.step + 1 - link.n_skipped;

        L = link;

        // Skip any holes at the start; there shouldn't be any.
        while (L.meas_idx >= n_meas && L.step != 0u) {
            L = links.at(L.previous_candidate_idx);
        }
    } else {
        need_to_write = false;
    }

    unsigned int current_state = 0;

    while (barrier.blockOr(current_state < num_states || need_to_write)) {
        if (current_state < num_states || need_to_write) {
            assert(L.meas_idx < n_meas);

            if (need_to_write) {
                vecmem::device_atomic_ref<unsigned long long int> mutex(
                    insertion_mutex.at(L.meas_idx));

                unsigned long long int assumed = mutex.load();
                unsigned long long int desired_set;
                auto [locked, size, worst] = decode_insertion_mutex(assumed);

                if (need_to_write &&
                    size >= payload.max_num_tracks_per_measurement &&
                    pval <= worst) {
                    need_to_write = false;
                }

                bool holds_lock = false;

                if (need_to_write && !locked) {
                    desired_set = encode_insertion_mutex(true, size, worst);

                    if (mutex.compare_exchange_strong(assumed, desired_set)) {
                        holds_lock = true;
                    }
                }

                if (holds_lock) {
                    unsigned int new_size;
                    unsigned int offset =
                        L.meas_idx * payload.max_num_tracks_per_measurement;
                    unsigned int out_idx;

                    if (size == payload.max_num_tracks_per_measurement) {
                        new_size = size;

                        typename algebra_t::scalar worst_pval =
                            std::numeric_limits<
                                typename algebra_t::scalar>::max();

                        out_idx = std::numeric_limits<unsigned int>::max();
                        for (unsigned int i = 0; i < size; ++i) {
                            if (tip_pval.at(offset + i) < worst_pval) {
                                worst_pval = tip_pval.at(offset + i);
                                out_idx = i;
                            }
                        }
                        assert(out_idx !=
                               std::numeric_limits<unsigned int>::max());
                    } else {
                        new_size = size + 1;
                        out_idx = size;
                    }

                    tip_index.at(offset + out_idx) = thread_id;
                    tip_pval.at(offset + out_idx) = pval;

                    float new_worst = std::numeric_limits<float>::max();

                    for (unsigned int i = 0; i < new_size; ++i) {
                        new_worst = std::min(
                            new_worst,
                            static_cast<float>(tip_pval.at(offset + i)));
                    }

                    [[maybe_unused]] bool cas_result =
                        mutex.compare_exchange_strong(
                            desired_set,
                            encode_insertion_mutex(false, new_size, new_worst));

                    assert(cas_result);

                    need_to_write = false;
                }
            }

            if (!need_to_write) {
                if (current_state < num_states - 1) {
                    L = links.at(L.previous_candidate_idx);
                    while (L.meas_idx >= n_meas && L.step != 0u) {
                        L = links.at(L.previous_candidate_idx);
                    }
                    need_to_write = true;
                } else {
#ifndef NDEBUG
                    if (L.step != 0) {
                        do {
                            L = links.at(L.previous_candidate_idx);
                        } while (L.meas_idx >= n_meas && L.step != 0u);
                        assert(L.meas_idx >= n_meas);
                    }
                    assert(L.step == 0);
#endif
                }

                current_state++;
            }
        }
    }
}

}  // namespace traccc::device
