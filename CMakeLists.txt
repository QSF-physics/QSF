cmake_minimum_required(VERSION 3.8)

# :::::::::::::::::::::::::::::::::: Project ::::::::::::::::::::::::::::::::::
# CMake 3.16 or higher provides native support for precompiled headers
# CMake 3.11 or higher provides FetchContent for downloading libraries from git
project(QSF VERSION 0.9.0 
        DESCRIPTION "Quantum Strong Field computations under one common framework"
        LANGUAGES C CXX)

set(QSF_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# ::::::::::::::::::::::::::::: Define imports ::::::::::::::::::::::::::::::::
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/findFFTW")
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


add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/extern/inipp")
# add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/extern/HighFive")

add_library(${PROJECT_NAME} INTERFACE)
# add_library(QSF::QSF ALIAS ${PROJECT_NAME})

# consider: http://videocortex.io/cmake_snippets/
# target_sources(${PROJECT_NAME} INTERFACE qsf.h) # Adds header files to target, will appear in IDEs
target_include_directories(${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${QSF_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>)

# maybe move to project
# set(CMAKE_CXX_STANDARD 17) # C++17 (for all targets)
# set_property(TARGET my_app PROPERTY CXX_STANDARD 17) # C++17
# target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
# special workaround for gcc version < 9 to link with std::filesystem libraries
target_link_libraries(${PROJECT_NAME} INTERFACE $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>)

target_link_libraries(${PROJECT_NAME} INTERFACE MPI::MPI_C)
# target_link_libraries(${PROJECT_NAME} INTERFACE MPI::MPI_CXX)
message(STATUS "FFTW_DOUBLE_LIB_FOUND: ${FFTW_DOUBLE_LIB_FOUND}") 
message(STATUS "FFTW_DOUBLE_MPI_LIB_FOUND: ${FFTW_DOUBLE_MPI_LIB_FOUND}") 
target_link_libraries(${PROJECT_NAME} INTERFACE FFTW::Double)
target_link_libraries(${PROJECT_NAME} INTERFACE FFTW::DoubleMPI)
# target_link_libraries(${PROJECT_NAME} INTERFACE fmt::fmt-header-only)
target_link_libraries(${PROJECT_NAME} INTERFACE inipp)
# target_link_libraries(${PROJECT_NAME} INTERFACE HighFive)







# find_package(xtensor REQUIRED)
# if(xtensor_FOUND)
#     message(STATUS "Found xtensor: ${xtensor_INCLUDE_DIRS}")
# else()
#     message(STATUS "xtensor not found :(")
# endif()

# find_package(xtensor-fftw REQUIRED)
# if(xtensor-fftw_FOUND)
#     message(STATUS "Found xtensor-fftw: ${xtensor-fftw_INCLUDE_DIRS}")
# else()
#     message(STATUS "xtensor-fftw not found :(")
# endif()
#find_package(xtensor-fftw REQUIRED)
#find_package(xtensor-io REQUIRED)

#--------------------------------------- cmake modules

# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/findFFTW/")
# add_library(${PROJECT_NAME} INTERFACE)
# find_package(MPI REQUIRED)
# if (NOT MPI_FOUND) 
#     message(STATUS "MPI not found :(") 
# else()
#     message(STATUS "MPI: ${MPI_INCLUDE_PATH}") 
# endif()
# target_link_libraries(${PROJECT_NAME} INTERFACE MPI::MPI_C)
# target_link_libraries(${PROJECT_NAME} INTERFACE MPI::MPI_CXX)


# find_package(PkgConfig REQUIRED)
# pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
# target_include_directories(${PROJECT_NAME} PkgConfig::FFTW)
# target_compile_options(${PROJECT_NAME} INTERFACE -lm -lmpi -lfftw3_mpi -lfftw3)
# target_link_libraries(${PROJECT_NAME} INTERFACE PkgConfig::FFTW)

# target_link_libraries(${PROJECT_NAME} INTERFACE fftw3-mpi)
# target_link_libraries(${PROJECT_NAME} INTERFACE fftw3)
# Using cmake (probably not good)
# find_package(FFTW3 REQUIRED)
# if (NOT FFTW3_FOUND) 
#     message(STATUS "FFTW3 not found :(") 
# else()
#     message(STATUS "FFTW3: ${FFTW3_INCLUDE_DIRS}") 
# endif()

# target_include_directories(${PROJECT_NAME} INTERFACE ${FFTW3_INCLUDE_DIRS})
# target_link_directories(${PROJECT_NAME} INTERFACE ${FFTW3_LIBRARY_DIRS})
# target_compile_options(${PROJECT_NAME} INTERFACE -lfftw3 -lm)

# target_link_libraries(${PROJECT_NAME} INTERFACE ${FFTW3_LIBRARY_DIRS})
# set(CMAKE_OSX_SYSROOT /opt/homebrew/)
# set(CMAKE_EXE_LINKER_FLAGS -Wl -undefined dynamic_lookup)
# set(CMAKE_SHARED_LINKER_FLAGS -Wl -undefined dynamic_lookup)
# if (TARGET MPI_C)
    # message(STATUS "MPI ${MPI_VERSION} found with MPI_C bindings")
    
    # target_include_directories(${PROJECT_NAME} INTERFACE ${MPI_INCLUDE_PATH})
# endif()
# if (TARGET MPI_CXX)
    # message(STATUS "MPI ${MPI_VERSION} found with MPI_CXX bindings")
    # target_include_directories(${PROJECT_NAME} INTERFACE ${MPI_INCLUDE_PATH})
# endif()
# message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")


# set(CMAKE_SHARED_LINKER_FLAGS -lfftw3_mpi,-lfftw3)

# find_package(FFTW3 REQUIRED)
# if(FFTW3_FOUND)
#     message(STATUS "Found FFTW3: ${FFTW3_INCLUDE_DIRS}")
#     target_link_libraries(${PROJECT_NAME} INTERFACE FFTW3)
# else()
#     message(STATUS "FFTW3 not found :(")
# endif()

# find_package(MPI REQUIRED)
# target_link_libraries(${PROJECT_NAME} INTERFACE MPI::MPI_CXX)


# message(STATUS "CMAKE_SKIP_BUILD_RPATH= ${CMAKE_SKIP_BUILD_RPATH}")
# message(STATUS "CMAKE_BUILD_WITH_INSTALL_RPATH= ${CMAKE_BUILD_WITH_INSTALL_RPATH}")
# message(STATUS "CMAKE_INSTALL_RPATH= ${CMAKE_INSTALL_RPATH}")
# message(STATUS "CMAKE_INSTALL_RPATH_USE_LINK_PATH= ${CMAKE_INSTALL_RPATH_USE_LINK_PATH}")
# message(STATUS "isSystemDir= ${isSystemDir}")

# include(GNUInstallDirs)
# set(CMAKE_SKIP_BUILD_RPATH FALSE)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
# if("${isSystemDir}" STREQUAL "-1")
#     set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# endif("${isSystemDir}" STREQUAL "-1")

# message(STATUS "CMAKE_SKIP_BUILD_RPATH= ${CMAKE_SKIP_BUILD_RPATH}")
# message(STATUS "CMAKE_BUILD_WITH_INSTALL_RPATH= ${CMAKE_BUILD_WITH_INSTALL_RPATH}")
# message(STATUS "CMAKE_INSTALL_RPATH= ${CMAKE_INSTALL_RPATH}")
# message(STATUS "CMAKE_INSTALL_RPATH_USE_LINK_PATH= ${CMAKE_INSTALL_RPATH_USE_LINK_PATH}")
# message(STATUS "isSystemDir= ${isSystemDir}")



# target_link_libraries(${PROJECT_NAME} INTERFACE xtensor)



# ::::::::::::::::::::::::::::Include guards ::::::::::::::::::::::::::::::::::

# if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
# message(FATAL_ERROR
#     "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
# )
# endif()

# if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.4.1" )
# message( STATUS "This compiler version might cause problems" )
# endif()

# if( APPLE )
#   message( STATUS "Running under MacOS X" )
# # Watch out, for this check also is TRUE under MacOS X because it
# # falls under the category of Unix-like.
# elseif( UNIX )
#   message( STATUS "Running under Unix or a Unix-like OS" )
# # Despite what you might think given this name, the variable is also
# # true for 64bit versions of Windows.
# elseif( WIN32 )
#   message( STATUS "Running under Windows (either 32bit or 64bit)" )
# endif()

# if( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
#   message( STATUS "g++ for the win!" )
#   message( STATUS ${CMAKE_CXX_COMPILER_VERSION} )
# endif()

# if( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
#   message( STATUS "LLVM, yeah!" )
# endif()


# if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
#     include(CTest)
# endif()

#cmake -S . -B build
#cmake --build build --target test

# add_library(QSF STATIC)

# find_package(MPI REQUIRED)
# message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
# target_link_libraries(MyTarget PUBLIC MPI::MPI_CXX)

#target_include_directories(QSF INTERFACE $<BUILD_INTERFACE:include> $<INSTALL_INTERFACE:include>)
#target_compile_features(QSF INTERFACE cxx_std_14)
#target_link_libraries(QSF INTERFACE xtl)

# add_library(QSF STATIC include/test.cpp)
# target_include_directories(QSF PUBLIC include)
# target_compile_features(QSF PUBLIC cxx_std_14)

#add_executable(calc projects/calc.cpp)
#target_link_libraries(calc PUBLIC QSF)