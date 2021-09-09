#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"
#include "Measurements2Spacepoints.hpp"

int main(int argc, char* argv[]) {
    return runMeasurements2SP(argc, argv, std::make_shared<ActsExamples::TGeoDetector>());
}