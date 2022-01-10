CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

# Set a default build type if none was specified
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")

# About this project
# ------------------------------------------------------------------------------
project(pfp-cst)
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "1")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

option(BUILD_TESTS "Set OFF to not compile the tests" ON)
option(BUILD_BENCHMARKS "Set OFF to not compile the benchmarks" ON)

# Set environment
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install)


find_package(Git)
if(GIT_FOUND)
    message("git found: ${GIT_EXECUTABLE}")
else()
    message(WARNING "git not found. Cloning of submodules will not work.")
endif()



# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include") # This is an hack because include(GUIInstallDirs) doesn't work

add_subdirectory(thirdparty)

# Configure the compiler with the appropriate flags
# ------------------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  include(ConfigureCompilerClang)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  include(ConfigureCompilerGcc)
else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message(FATAL_ERROR "Only the compiler gcc and clang are supported")
endif()


add_subdirectory(include)
add_subdirectory(src)

if(${BUILD_TESTS})
  message(STATUS "Building tests enabled.")
  add_subdirectory(test)
endif()

if(${BUILD_BENCHMARKS})
  message(STATUS "Building benchmarks enabled.")
  add_subdirectory(benchmarks/src)
endif()

# Configure pipeline for build folder
set(USE_INSTALL_PATH False)
set(PFP_CST_VERSION False)
configure_file(${PROJECT_SOURCE_DIR}/pipeline/pfp-cst.in ${PROJECT_BINARY_DIR}/pfp-cst @ONLY)

# Configure pipeline for install folder
set(USE_INSTALL_PATH True)
set(PFP_CST_VERSION "0.1.0")
configure_file(${PROJECT_SOURCE_DIR}/pipeline/pfp-cst.in ${PROJECT_BINARY_DIR}/pfp-cst.install @ONLY)


install(TARGETS pfp-cst_build pfp-cst_build64 TYPE RUNTIME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/pfp-cst.install RENAME pfp-cst TYPE BIN)

# # Installin the library https://dominikberner.ch/cmake-interface-lib/
# install(TARGETS pfp_cst pfp common gsacak malloc_count
#         EXPORT ${PROJECT_NAME}_Targets
#         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# include(CMakePackageConfigHelpers)
# write_basic_package_version_file("pfp-cstConfigVersion.cmake"
#                                  VERSION ${VERSION}
#                                  COMPATIBILITY SameMajorVersion)

# install(EXPORT ${PROJECT_NAME}_Targets
#         FILE ${PROJECT_NAME}Targets.cmake
#         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

# install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
#               "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
#         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

# install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/pfp 
#                   ${PROJECT_SOURCE_DIR}/include/common 
#                   DESTINATION include)
# Configure cpack variables
# ------------------------------------------------------------------------------

set(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set(CPACK_PACKAGE_VERSION "${VERSION}")

include(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "STGZ;TGZ;DEB")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "University of Florida")
set(CPACK_PACKAGE_CONTACT "rossi.m@ufl.edu")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "PFP-CST")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-sources")

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Massimiliano Rossi")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT) 
set(CPACK_COMPONENTS_GROUPING ALL_COMPONENTS_IN_ONE) # Groupp all components
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEB_COMPONENT_INSTALL YES)
include(CPack)