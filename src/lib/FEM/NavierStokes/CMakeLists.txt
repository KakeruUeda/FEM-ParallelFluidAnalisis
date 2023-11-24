set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(NavierStokes STATIC)
target_sources(NavierStokes
  PRIVATE
  CreateMatrix.cpp
  NavierStokes.cpp
)

target_include_directories(NavierStokes PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
  ${PROJECT_SOURCE_DIR}/src/lib/FEM
)

target_compile_options(NavierStokes PUBLIC ${compile_options_tmp})