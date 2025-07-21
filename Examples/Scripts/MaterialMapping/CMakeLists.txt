add_executable(ActsAnalysisMaterialComposition MaterialComposition.cpp)
target_link_libraries(
    ActsAnalysisMaterialComposition
    Acts::ExamplesFramework
    ROOT::Core
    ROOT::Hist
    ROOT::Tree
    ROOT::TreePlayer
    Boost::program_options
    nlohmann_json::nlohmann_json
)
