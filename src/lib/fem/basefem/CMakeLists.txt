set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(fem STATIC)
target_sources(fem
  PRIVATE
  FEM.cpp
  BaseFEM.cpp
  DomainFEM.cpp
  MathFEM.cpp
  SolutionData.cpp
  ElementBaseFEM.cpp
  VisualizeVariables.cpp
)

target_include_directories(fem PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
  ${PROJECT_SOURCE_DIR}/src/lib/fem/postprocess
)

target_compile_options(fem PUBLIC ${compile_options_tmp})