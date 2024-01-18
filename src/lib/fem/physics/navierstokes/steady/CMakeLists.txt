set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(steadynavierstokes STATIC)
target_sources(steadynavierstokes
  PRIVATE
  SteadyNavierStokes.cpp
  MatrixAssemblySNS.cpp
  XFEM_MatrixAssemblySNS.cpp
  Darcy_MatrixAssemblySNS.cpp
)

target_include_directories(steadynavierstokes PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/fem/basefem
  ${PROJECT_SOURCE_DIR}/src/lib/fem/postprocess
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
)

target_compile_options(steadynavierstokes PUBLIC ${compile_options_tmp})