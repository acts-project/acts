target_sources(
    ActsCore
    PRIVATE Common.cpp Direction.cpp ParticleData.cpp ParticleIdHelper.cpp
)

include(ActsCodegen)

set(PARTICLE_TABLE_CODEGEN_PYTHON_ARGS "")

if(DEFINED ACTS_ADDITIONAL_PARTICLES_CSV AND NOT ACTS_ADDITIONAL_PARTICLES_CSV STREQUAL "")
    list(APPEND PARTICLE_TABLE_CODEGEN_PYTHON_ARGS --additional-particles "${ACTS_ADDITIONAL_PARTICLES_CSV}")
endif()

acts_code_generation(
    ADD_TO_TARGET ActsCore
    PYTHON ${PROJECT_SOURCE_DIR}/Fatras/scripts/generate_particle_data_table.py
    PYTHON_ARGS ${PARTICLE_TABLE_CODEGEN_PYTHON_ARGS}
    OUTPUT ParticleDataTable.hpp
)
