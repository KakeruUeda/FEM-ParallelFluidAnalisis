set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(unsteadynavierstokes STATIC)
target_sources(unsteadynavierstokes
  PRIVATE
  UnsteadyNavierStokes.cpp
  MatrixAssemblyUSNS.cpp
  Darcy_MatrixAssemblyUSNS.cpp
  XFEM_MatrixAssemblyUSNS.cpp
)

target_include_directories(unsteadynavierstokes PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/fem/postprocess
  ${PROJECT_SOURCE_DIR}/src/lib/fem/basefem
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
)

target_compile_options(unsteadynavierstokes PUBLIC ${compile_options_tmp})