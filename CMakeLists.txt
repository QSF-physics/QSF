cmake_minimum_required(VERSION 3.8)

# :::::::::::::::::::::::::::::::::: Project ::::::::::::::::::::::::::::::::::
# CMake 3.16 or higher provides native support for precompiled headers
# CMake 3.11 or higher provides FetchContent for downloading libraries from git
project(QSF VERSION 0.9.0 
        DESCRIPTION "Quantum Strong Field computations under one common framework"
        LANGUAGES C CXX)

set(QSF_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
# ::::::::::::::::::::::::::::: Define imports ::::::::::::::::::::::::::::::::
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
message(STATUS "CMAKE_MODULE_PATH= ${CMAKE_MODULE_PATH}")

find_package(MPI REQUIRED)
if (NOT MPI_FOUND) 
    message(STATUS "MPI not found :(") 
else()
    message(STATUS "MPI: ${MPI_INCLUDE_PATH}") 
endif()
if(MPI_FOUND)
  include(CMakePrintHelpers)
  cmake_print_properties(TARGETS MPI::MPI_C PROPERTIES INTERFACE_LINK_LIBRARIES)
endif()
find_package(FFTW REQUIRED DOUBLE_MPI_LIB DOUBLE_LIB)
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${QSF_INCLUDE_DIR}> $<INSTALL_INTERFACE:include>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
# special workaround for gcc version < 9 to link with std::filesystem libraries
target_link_libraries(${PROJECT_NAME} INTERFACE $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>)
target_link_libraries(${PROJECT_NAME} INTERFACE MPI::MPI_C)
message(STATUS "FFTW_DOUBLE_LIB_FOUND: ${FFTW_DOUBLE_LIB_FOUND}") 
message(STATUS "FFTW_DOUBLE_MPI_LIB_FOUND: ${FFTW_DOUBLE_MPI_LIB_FOUND}") 
target_link_libraries(${PROJECT_NAME} INTERFACE FFTW::Double)
target_link_libraries(${PROJECT_NAME} INTERFACE FFTW::DoubleMPI)

add_subdirectory ("docs")