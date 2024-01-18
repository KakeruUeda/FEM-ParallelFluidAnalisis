set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(postprocess STATIC)
target_sources(postprocess
  PRIVATE
  PrepareForDA.cpp
  ExportFile.cpp
)

target_include_directories(postprocess PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
  ${PROJECT_SOURCE_DIR}/src/lib/fem/basefem
  ${PROJECT_SOURCE_DIR}/src/lib/fem/postprocess
)

target_compile_options(postprocess PUBLIC ${compile_options_tmp})