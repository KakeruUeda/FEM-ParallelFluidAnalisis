set(main_target "SteadyStokesSolver")

set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_executable(${main_target} main.cpp)

target_sources(${main_target}
  PRIVATE
    main.cpp
)

target_compile_options(${main_target} PUBLIC ${compile_options_tmp})


target_link_libraries( ${main_target} PRIVATE
  TPmpi
  petsc
  metis
  base
  fem
  solver
  stokes
  preprocess
  postprocess
  omp
)

target_include_directories(${main_target} PRIVATE
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/fem/basefem
  ${PROJECT_SOURCE_DIR}/src/lib/solver
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})