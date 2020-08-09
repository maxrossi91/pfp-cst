
add_executable(pfp_cst_benchmarks pfp_cst_benchmarks.cpp)
target_link_libraries(pfp_cst_benchmarks common strdup pfp gsacak sdsl divsufsort divsufsort64 malloc_count benchmark pthread)
target_include_directories(pfp_cst_benchmarks PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# target_compile_options(ds_building_wt_sdsl_264 PUBLIC -DM64)