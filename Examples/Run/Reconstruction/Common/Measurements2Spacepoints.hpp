#pragma once

#include <memory>
namespace ActsExamples {
class IBaseDetector;
}

int runMeasurements2SP(int argc, char* argv[],
    std::shared_ptr<ActsExamples::IBaseDetector> detector);