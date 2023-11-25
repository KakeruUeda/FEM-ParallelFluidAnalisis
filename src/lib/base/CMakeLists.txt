set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(Base STATIC)
target_sources(Base PRIVATE
    basicFunctions.cpp
)


target_compile_options(Base PUBLIC ${compile_options_tmp})