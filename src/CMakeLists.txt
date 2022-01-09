add_executable(ds_building ds_building.cpp)
target_link_libraries(ds_building common strdup pfp gsacak sdsl malloc_count)

add_executable(ds_building64 ds_building.cpp)
target_link_libraries(ds_building64 common strdup pfp gsacak64 sdsl malloc_count)
target_compile_options(ds_building64 PUBLIC -DM64)