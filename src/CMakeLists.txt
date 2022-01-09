add_executable(pfp-cst_build pfp-cst_build.cpp)
target_link_libraries(pfp-cst_build common strdup pfp gsacak sdsl malloc_count)

add_executable(pfp-cst_build64 pfp-cst_build.cpp)
target_link_libraries(pfp-cst_build64 common strdup pfp gsacak64 sdsl malloc_count)
target_compile_options(pfp-cst_build64 PUBLIC -DM64)