# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

#
# Set up the "throughput applications".
#
traccc_add_library( traccc_examples_cpu examples_cpu SHARED
   "include/traccc/examples/cpu/full_chain_algorithm.hpp"
   "src/full_chain_algorithm.cpp" )
target_link_libraries( traccc_examples_cpu
   PUBLIC vecmem::core detray::core detray::detectors traccc::core
   traccc::examples_common )

traccc_add_executable( throughput_st "apps/throughput_st.cpp"
   LINK_LIBRARIES vecmem::core detray::detectors
                  detray::io traccc::core traccc::io traccc::performance
                  traccc::options traccc::examples_cpu )

traccc_add_executable( throughput_mt "apps/throughput_mt.cpp"
   LINK_LIBRARIES TBB::tbb vecmem::core detray::detectors
                  detray::io traccc::core traccc::io traccc::performance
                  traccc::options traccc::examples_cpu )

traccc_add_executable( seeding_example "apps/seeding_example.cpp"
   LINK_LIBRARIES vecmem::core traccc::core traccc::io
   traccc::performance traccc::options detray::detectors detray::io
   traccc::examples_common )

traccc_add_executable( seq_example "apps/seq_example.cpp"
   LINK_LIBRARIES vecmem::core traccc::core traccc::io
   traccc::performance traccc::options detray::detectors detray::io
   traccc::examples_common )

traccc_add_executable( truth_finding_example "apps/truth_finding_example.cpp"
   LINK_LIBRARIES vecmem::core detray::detectors traccc::core traccc::io
   traccc::performance traccc::options traccc::examples_common )

traccc_add_executable( truth_fitting_example "apps/truth_fitting_example.cpp"
   LINK_LIBRARIES vecmem::core detray::io detray::detectors traccc::core
   traccc::io traccc::performance traccc::options traccc::examples_common )

traccc_add_executable( misaligned_truth_fitting_example "apps/misaligned_truth_fitting_example.cpp"
   LINK_LIBRARIES vecmem::core detray::io detray::detectors traccc::core
   traccc::io traccc::performance traccc::options traccc::examples_common )

traccc_add_executable( kalman_filter_validation "apps/kalman_filter_validation.cpp"
   LINK_LIBRARIES vecmem::core detray::io detray::detectors traccc::core
   traccc::io traccc::simulation traccc::performance traccc::options traccc_examples_common )
