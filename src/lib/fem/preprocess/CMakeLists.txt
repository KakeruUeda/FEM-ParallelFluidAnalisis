set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(preprocess STATIC)
target_sources(preprocess
  PRIVATE
  InitializeFEM.cpp
  ReadlnputFile.cpp
  PrepareMatrix.cpp
  AllocationFEM.cpp
  OctreeSubDivision.cpp
  InterfacePartition.cpp
  LineSDFSearch.cpp
)

target_include_directories(preprocess PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
  ${PROJECT_SOURCE_DIR}/src/lib/fem/basefem
  ${PROJECT_SOURCE_DIR}/src/lib/fem/postprocess
)

target_compile_options(preprocess PUBLIC ${compile_options_tmp})