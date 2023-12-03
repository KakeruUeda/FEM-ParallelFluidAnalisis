set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(base STATIC)
target_sources(base PRIVATE
    BasicFunctions.cpp
)

target_compile_options(base PUBLIC ${compile_options_tmp})