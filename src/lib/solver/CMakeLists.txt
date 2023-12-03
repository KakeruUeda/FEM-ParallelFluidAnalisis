set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(Solver STATIC)
target_sources(Solver
  PRIVATE
  PetscSolver.cpp
)

target_include_directories(Solver PRIVATE
  ${PROJECT_SOURCE_DIR}/src/lib/base
)

target_compile_options(Solver PUBLIC ${compile_options_tmp})