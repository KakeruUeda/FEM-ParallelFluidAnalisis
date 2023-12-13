set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(navierstokes STATIC)
target_sources(navierstokes
  PRIVATE
  CreateMatrix_SNS.cpp
  CreateMatrixXFEM_SNS.cpp
  NavierStokes.cpp
)

target_include_directories(navierstokes PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
  ${PROJECT_SOURCE_DIR}/src/lib/fem/basefem
)

target_compile_options(navierstokes PUBLIC ${compile_options_tmp})