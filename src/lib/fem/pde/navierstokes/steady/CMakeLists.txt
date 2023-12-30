set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(steadynavierstokes STATIC)
target_sources(steadynavierstokes
  PRIVATE
  MatrixAssemblySNS.cpp
  XFEM_MatrixAssemblySNS.cpp
  SteadyNavierStokes.cpp
)

target_include_directories(steadynavierstokes PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
  ${PROJECT_SOURCE_DIR}/src/lib/fem/basefem
)

target_compile_options(steadynavierstokes PUBLIC ${compile_options_tmp})