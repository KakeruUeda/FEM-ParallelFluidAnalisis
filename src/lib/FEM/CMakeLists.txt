set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(fem STATIC)
target_sources(fem
  PRIVATE
  FEM.cpp
  BaseFEM.cpp
  InitializeFEM.cpp
  DomainFEM.cpp
  MathFEM.cpp
  SolutionData.cpp
  ElementBaseFEM.cpp
  ReadlnputFile.cpp
  PrepareMatrix.cpp
)

target_include_directories(fem PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
)

target_compile_options(fem PUBLIC ${compile_options_tmp})

add_subdirectory(navierstokes)
add_subdirectory(stokes)