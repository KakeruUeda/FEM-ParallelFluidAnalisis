set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(solver STATIC)
target_sources(solver
  PRIVATE
  PetscSolver.cpp
)

target_include_directories(solver PRIVATE
  ${PROJECT_SOURCE_DIR}/src/lib/base
)

target_compile_options(solver PUBLIC ${compile_options_tmp})