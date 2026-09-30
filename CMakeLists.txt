cmake_minimum_required(VERSION 3.25)
project(dune-sycl-test
        VERSION 0.0.1
        LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_COLOR_DIAGNOSTICS ON)

set(CMAKE_FIND_PACKAGE_TARGETS_GLOBAL ON)
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR}/extern)

add_subdirectory(extern/dune-common         SYSTEM)
add_subdirectory(extern/dune-geometry       SYSTEM)
add_subdirectory(extern/dune-grid           SYSTEM)
add_subdirectory(extern/dune-localfunctions SYSTEM)
add_subdirectory(extern/dune-istl           SYSTEM)

find_package(dune-common         REQUIRED)
find_package(dune-geometry       REQUIRED)
find_package(dune-grid           REQUIRED)
find_package(dune-localfunctions REQUIRED)
find_package(dune-istl           REQUIRED)

find_package(AdaptiveCpp         REQUIRED)

add_executable(main main.cc)
target_link_libraries(main PRIVATE Dune::ISTL Dune::LocalFunctions Dune::Grid Dune::Geometry Dune::Common)
add_sycl_to_target(TARGET main SOURCES main.cc)
