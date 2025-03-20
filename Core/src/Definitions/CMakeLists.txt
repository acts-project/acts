target_sources(ActsCore PRIVATE Common.cpp Direction.cpp ParticleData.cpp)

include(ActsCodegen)

if(ACTS_USE_SYSTEM_LIBS)
    set(_codegen_offline_flags "OFFLINE")
else()
    set(_codegen_offline_flags "")
endif()

acts_code_generation(
      TARGET ParticleDataTable
      PYTHON ${CMAKE_SOURCE_DIR}/Fatras/scripts/generate_particle_data_table.py
      WITH_REQUIREMENTS ${CMAKE_SOURCE_DIR}/Fatras/scripts/requirements.txt
      OUTPUT ParticleDataTable.hpp
      ${_codegen_offline_flags}
)

unset(_codegen_offline_flags)

target_link_libraries(ActsCore PRIVATE ParticleDataTable)
