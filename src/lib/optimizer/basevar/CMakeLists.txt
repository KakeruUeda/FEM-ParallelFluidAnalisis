set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(basevar STATIC)
target_sources(basevar
  PRIVATE
  VariationalDataAssimilation.cpp
)

target_include_directories(basevar PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/solver
  ${PROJECT_SOURCE_DIR}/src/lib/fem/basefem
  ${PROJECT_SOURCE_DIR}/src/lib/fem/postprocess
)

target_compile_options(fem PUBLIC ${compile_options_tmp})