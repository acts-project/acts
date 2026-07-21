/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/options/threading.hpp"

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <stdexcept>

namespace traccc::opts {

/// Type alias for concurrent slots option
using concurrent_slots_type = std::size_t;
/// Name of the concurrent slots option
static const char* concurrent_slots_option = "concurrent-slots";

threading::threading() : interface("Multi-Threading Options") {
  m_desc.add_options()(
      "cpu-threads",
      boost::program_options::value(&threads)->default_value(threads),
      "The number of CPU threads to use")(
      concurrent_slots_option,
      boost::program_options::value<concurrent_slots_type>(),
      "The number of events that can be "
      "processed concurrently, be default equal to cpu-threads");
}

void threading::read(const boost::program_options::variables_map& vm) {
  if (threads == 0) {
    throw std::invalid_argument{"Must use threads>0"};
  }
  if (!vm.count(concurrent_slots_option)) {
    concurrent_slots = threads;
  } else {
    concurrent_slots = vm[concurrent_slots_option].as<concurrent_slots_type>();
    if (concurrent_slots == 0) {
      throw std::invalid_argument{"Must use concurrent-slots>0"};
    }
  }
}

std::unique_ptr<configuration_printable> threading::as_printable() const {
  auto cat = std::make_unique<configuration_category>(m_description);

  cat->add_child(std::make_unique<configuration_kv_pair>(
      "Number of CPU thread", std::to_string(threads)));
  cat->add_child(std::make_unique<configuration_kv_pair>(
      "Number of concurrent event slots", std::to_string(concurrent_slots)));
  return cat;
}

}  // namespace traccc::opts
