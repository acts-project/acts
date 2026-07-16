/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <vecmem/containers/data/vector_view.hpp>

#include "traccc/device/global_index.hpp"
#include "traccc/edm/device/sort_key.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/candidate_link.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/utils/prob.hpp"

namespace traccc::device {

/**
 * @brief Payload for the duplicate-removal kernel.
 */
struct remove_duplicates_payload {
    /**
     * @brief View of the links vector.
     */
    const vecmem::data::vector_view<const candidate_link> links_view;

    /**
     * @brief View of the last measurement in each link of the current step.
     */
    const vecmem::data::vector_view<const unsigned int>
        link_last_measurement_view;

    /**
     * @brief View of the parameter IDs, sorted in order of last measurement
     * index.
     */
    const vecmem::data::vector_view<const unsigned int> param_ids_view;

    /**
     * @brief View to the vector of parameter livenesses, which is used to
     * remove duplicate tracks.
     */
    vecmem::data::vector_view<unsigned int> param_liveness_view;

    /**
     * @brief The total number of links in this step.
     */
    const unsigned int n_links;

    /**
     * @brief The starting index of the links belonging to this step.
     */
    const unsigned int curr_links_idx;

    /**
     * @brief The total number of measurements, used to find holes.
     */
    const unsigned int n_measurements;

    /**
     * @brief The current step.
     */
    const unsigned int step;
};

/**
 * @brief The CKF duplicate removal kernel.
 *
 * This kernel is designed to remove duplicate tracks from the CKF algorithm.
 * As a motivating example, consider a truth track consisting of spacepoints
 * A, B, C, D, and E. During seed finding, we find seeds 1: (A, B, C) and
 * 2: (A, B, D). Although the seeds are different, they will eventually end up
 * with exactly the same spacepoints, and thus one of them can be removed.
 *
 * This kernel detects tracks with the same spacepoints and keeps only the one
 * with the lowest $\chi^2$ score.
 *
 * @note For performance reasons, this kernel can only deduplicate tracks with
 * the last measurement the same. Thus, while a track (A, B, C, D) might be
 * able to replace a shorter track (A, B, C), this algorithm will not do so.
 * This is not a large problem, however, as tracks with a different initial
 * spacepoint are much more likely to have significantly different initial
 * track parameters.
 *
 * @warning If this kernel is run very early in the CKF process, it might end
 * up removing tracks that are the same before they begin the branch. The very
 * minimum length at which tracks can be removed is 3 spacepoints, but the
 * minimum length should be configured conform physics performance
 * requirements. Intuitively, longer tracks with the same spacepoints have
 * undergone more Kálmán gain to bring them closer together, and thus it is
 * increasingly less likely that they will diverge afterwards.
 */
TRACCC_HOST_DEVICE inline void remove_duplicates(
    global_index_t tid, const finding_config& cfg,
    const remove_duplicates_payload& payload) {

    const vecmem::device_vector<const candidate_link> links(payload.links_view);
    vecmem::device_vector<unsigned int> param_liveness(
        payload.param_liveness_view);
    const vecmem::device_vector<const unsigned int> link_last_measurement(
        payload.link_last_measurement_view);
    const vecmem::device_vector<const unsigned int> param_ids(
        payload.param_ids_view);

    /*
     * As is standard fare, we ignore tracks that are out of bounds or that
     * have already been marked as "dead". Since this kernel contains no
     * synchronization points, early returns are safe.
     */
    if (tid >= param_ids.size()) {
        return;
    }

    const unsigned int this_param_id = param_ids.at(tid);

    if (param_liveness.at(this_param_id) == 0u) {
        return;
    }

    const unsigned int last_measurement = link_last_measurement.at(tid);

    /*
     * We will now gather the range of tracks with the same last measurement
     * as the track belonging to this thread. Because the tracks are sorted by
     * the last measurement, we can perform a simple linear search.
     */
    unsigned int min_tid = tid;
    unsigned int max_tid = tid;

    while (min_tid > 0 &&
           link_last_measurement.at(min_tid - 1) == last_measurement) {
        min_tid--;
    }

    while (max_tid < link_last_measurement.size() - 1 &&
           link_last_measurement.at(max_tid + 1) == last_measurement) {
        max_tid++;
    }

    /*
     * Throughout this function, we refer to "this" as the track uniquely
     * belonging to the executing thread and "that" as any other track
     * considered.
     */
    const candidate_link& Lthisbase =
        links.at(payload.curr_links_idx + param_ids.at(tid));

    /*
     * A core design of this function is that a thread can only remove itself,
     * and a track cannot be removed if it has too few measurements. Thus, we
     * return early if the thread under study is too short.
     */
    if (payload.step + 1 - Lthisbase.n_skipped <=
            cfg.duplicate_removal_minimum_length ||
        Lthisbase.ndf_sum <= 5) {
        return;
    }

    const scalar prob_this =
        prob(Lthisbase.chi2_sum, static_cast<scalar>(Lthisbase.ndf_sum - 5));

    /*
     * We now compare the current track with every other track that has the
     * same final measurement.
     */
    for (unsigned int i = min_tid; i <= max_tid; ++i) {
        /*
         * Obviously, a thread cannot be compared to itself.
         */
        if (i == tid)
            continue;

        /*
         * We will consider a track dominated if every single one of its
         * track states is also found in another track and the other way
         * around.
         */
        bool this_is_dominated = true;

        candidate_link Lthis = Lthisbase;
        candidate_link Lthat =
            links.at(payload.curr_links_idx + param_ids.at(i));

        if (payload.step + 1 - Lthat.n_skipped <=
                cfg.duplicate_removal_minimum_length ||
            Lthis.ndf_sum <= 5 || Lthat.ndf_sum <= 5) {
            continue;
        }

        const scalar prob_that =
            prob(Lthat.chi2_sum, static_cast<scalar>(Lthat.ndf_sum - 5));

        /*
         * This loop is the main workhorse of the algorithm, comparing the
         * track states of the two tracks.
         */
        while (true) {
            /*
             * First, we eliminate any holes in either track, to ensure that
             * we only compare actual measurements.
             */
            while (Lthis.meas_idx >= payload.n_measurements &&
                   Lthis.step != 0u) {
                Lthis = links.at(Lthis.previous_candidate_idx);
            }

            while (Lthat.meas_idx >= payload.n_measurements &&
                   Lthat.step != 0u) {
                Lthat = links.at(Lthat.previous_candidate_idx);
            }

            /*
             * If the two measurements at the top of the recursively examined
             * subtracks are the same, our search continues.
             */
            if (Lthis.meas_idx == Lthat.meas_idx) {
                /*
                 * If we have reached the last step in the "this" track, then
                 * we know for certain that the tracks are the same and thus
                 * we break out of the loop. Note that not setting the
                 * domination flag to false leaves it true, thus marking this
                 * track as dominated.
                 */
                if (Lthis.step == 0) {
                    break;
                }
                /*
                 * If, instead, the other track has run out then we know the
                 * current track has a state that the other does not, so it is
                 * definitionally not dominated.
                 */
                else if (Lthat.step == 0) {
                    this_is_dominated = false;
                    break;
                }
                /*
                 * Finally, in all other cases we still have some track to go
                 * and we recurse further.
                 */
                else {
                    Lthis = links.at(Lthis.previous_candidate_idx);
                    Lthat = links.at(Lthat.previous_candidate_idx);
                }
            }
            /*
             * If the two measurements at the top are different, we can be
             * certain that the tracks are different, so the current track is
             * _not_ dominated by the other.
             */
            else {
                this_is_dominated = false;
                break;
            }
        }

        /*
         * We now decide which of the tracks is better. Note that this block
         * of code does nothing if the domination flag is false. If it is true
         * it will remain true only if either the other track has a better
         * p-value or, if the p-value is equal, if the other track has a lower
         * identifier (this should be a very uncommon tie-breaker).
         */
        if (prob_this != prob_that) {
            this_is_dominated &= prob_that >= prob_this;
        } else {
            this_is_dominated &= param_ids.at(i) < param_ids.at(tid);
        }

        /*
         * Finally, if we determine the track to be dominated we mark it as
         * such in global memory and we exit out of the main loop, thus
         * exiting the kernel.
         */
        if (this_is_dominated) {
            param_liveness.at(param_ids.at(tid)) = 0u;
            break;
        }
    }
}

}  // namespace traccc::device
