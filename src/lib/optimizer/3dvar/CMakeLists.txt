set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(3dvar STATIC)
target_sources(3dvar
  PRIVATE
  AdjointVariableMethod.cpp
  AdjointSteadyNavierStokes.cpp
  AdjointMatrixAssemblySNS.cpp
  Darcy_AdjointMatrixAssemblySNS.cpp
)

target_include_directories(3dvar PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
  ${PROJECT_SOURCE_DIR}/src/lib/fem/basefem
  ${PROJECT_SOURCE_DIR}/src/lib/fem/postprocess
  ${PROJECT_SOURCE_DIR}/src/lib/optimizer/basevar
)

target_compile_options(preprocess PUBLIC ${compile_options_tmp})