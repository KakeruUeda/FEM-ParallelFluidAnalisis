set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(stokes STATIC)
target_sources(stokes
  PRIVATE
  CreateMatrix.cpp
  CreateMatrixXFEM.cpp
  Stokes.cpp
)

target_include_directories(stokes PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
  ${PROJECT_SOURCE_DIR}/src/lib/fem
)

target_compile_options(stokes PUBLIC ${compile_options_tmp})