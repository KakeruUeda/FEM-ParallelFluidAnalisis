set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(Stokes STATIC)
target_sources(Stokes
  PRIVATE
  CreateMatrix.cpp
  Stokes.cpp
)

target_include_directories(Stokes PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
  ${PROJECT_SOURCE_DIR}/src/lib/FEM
)

target_compile_options(Stokes PUBLIC ${compile_options_tmp})