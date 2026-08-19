target_sources(
    ActsCore
    PRIVATE Common.cpp Direction.cpp ParticleData.cpp ParticleIdHelper.cpp
)

include(ActsCodegen)

acts_code_generation(
    ADD_TO_TARGET ActsCore
    KEY Core/src/Definitions/ParticleDataTable.hpp
)
