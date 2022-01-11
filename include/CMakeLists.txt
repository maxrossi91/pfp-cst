add_subdirectory(common)
target_include_directories(common PUBLIC common)

add_subdirectory(pfp)
target_include_directories(pfp PUBLIC pfp)

add_library(pfp_cst OBJECT pfp/pfp_cst.hpp)
target_link_libraries(pfp_cst pfp common gsacak malloc_count)
target_include_directories(pfp_cst PUBLIC pfp)