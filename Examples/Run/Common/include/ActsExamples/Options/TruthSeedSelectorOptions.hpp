#pragma once

#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {
namespace Options {

void addTruthSeedSelectorOptions(Options::Description& desc);

ActsExamples::TruthSeedSelector::Config
readTruthSeedSelectorConfig(const Options::Variables& vars);

}
}
