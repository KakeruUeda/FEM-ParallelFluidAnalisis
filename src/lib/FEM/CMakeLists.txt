set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(FEM STATIC)
target_sources(FEM
  PRIVATE
  baseFEM.cpp
  initializeFEM.cpp
  exportFEM.cpp
  domainFEM.cpp
  mathFEM.cpp
  solutionData.cpp
  elementBase.cpp
)

target_include_directories(FEM PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
)

target_compile_options(FEM PUBLIC ${compile_options_tmp})

add_subdirectory(NavierStokes)
add_subdirectory(Stokes)