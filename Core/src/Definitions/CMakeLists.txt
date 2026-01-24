target_sources(
    ActsCore
    PRIVATE Common.cpp Direction.cpp ParticleData.cpp ParticleIdHelper.cpp
)

include(ActsCodegen)

acts_code_generation(
    ADD_TO_TARGET ActsCore
    PYTHON ${PROJECT_SOURCE_DIR}/Fatras/scripts/generate_particle_data_table.py
    OUTPUT ParticleDataTable.hpp
)
