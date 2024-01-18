set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(stokes STATIC)
target_sources(stokes
  PRIVATE
  SteadyStokes.cpp
  MatrixAssemblySST.cpp
  XFEM_MatrixAssemblySST.cpp
  Darcy_MatrixAssemblySST.cpp
)

target_include_directories(stokes PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/fem/postprocess
  ${PROJECT_SOURCE_DIR}/src/lib/fem/basefem
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
)

target_compile_options(stokes PUBLIC ${compile_options_tmp})