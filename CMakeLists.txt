cmake_minimum_required(VERSION 3.6)
project(GEDc++)


set(CMAKE_CXX_STANDARD 17)
option(CXX "enable C++ compilation" ON)
if (CXX)
    enable_language(CXX)
endif ()


# a CMake module named "FindGUROBI.cmake" is available in cmake/modules/
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

find_package(GUROBI REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS})


add_subdirectory(src/utils)
add_subdirectory(src/executables)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


if (DEFINED ENV{GEDLIB_ROOT})
    message(STATUS "GEDLIB_ROOT in ENV I see is: $ENV{GEDLIB_ROOT}")
endif ()
if (NOT DEFINED ENV{GEDLIB_ROOT})
    message(STATUS "!!! GEDLIB_ROOT not set, please add GEDLIB_ROOT=~/gedlib to your cmake settings (e.g. in CLion), or via -DGEDLIB_ROOT=~/gedlib !!!")
endif ()

if (DEFINED ENV{LIBLSAP_ROOT})
    message(STATUS "LIBLSAP_ROOT in ENV I see is: $ENV{LIBLSAP_ROOT}")
endif ()
if (NOT DEFINED ENV{LIBLSAP_ROOT})
    message(STATUS "!!! LIBLSAP_ROOT not set, please add LIBLSAP_ROOT=~/liblsap to your cmake settings (e.g. in CLion), or via -DLIBLSAP_ROOT=~/liblsap !!!")
endif ()


add_executable(prepro_verification_aids
        src/executables/prepro_gurobi_aids_similarity.cpp
        src/utils/GXLGraphReader.cpp
        src/utils/io.cpp
)

target_include_directories(prepro_verification_aids
        PRIVATE
        include
        $ENV{LIBLSAP_ROOT}/cpp/include
        $ENV{GEDLIB_ROOT}
        $ENV{GEDLIB_ROOT}/ext/boost.1.69.0
        $ENV{GEDLIB_ROOT}/ext/eigen.3.3.4/Eigen
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/src
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/ext/sgtelib/src
        $ENV{GEDLIB_ROOT}/ext/lsape.5/include
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/include
)
target_link_directories(prepro_verification_aids
        PRIVATE
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/lib
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/lib
        $ENV{GEDLIB_ROOT}/lib
)
target_link_libraries(prepro_verification_aids
        PRIVATE
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        libdoublefann.so.2
        libsvm.so
        libnomad.so
        libgxlgedlib.so
)

add_executable(prepro_verification_prot
        src/executables/prepro_gurobi_prot_similarity.cpp
        src/utils/GXLGraphReader.cpp
        src/utils/io.cpp
)

target_include_directories(prepro_verification_prot
        PRIVATE
        include
        $ENV{LIBLSAP_ROOT}/cpp/include
        $ENV{GEDLIB_ROOT}
        $ENV{GEDLIB_ROOT}/ext/boost.1.69.0
        $ENV{GEDLIB_ROOT}/ext/eigen.3.3.4/Eigen
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/src
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/ext/sgtelib/src
        $ENV{GEDLIB_ROOT}/ext/lsape.5/include
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/include
)
target_link_directories(prepro_verification_prot
        PRIVATE
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/lib
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/lib
        $ENV{GEDLIB_ROOT}/lib
)
target_link_libraries(prepro_verification_prot
        PRIVATE
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        libdoublefann.so.2
        libsvm.so
        libnomad.so
        libgxlgedlib.so
)

add_executable(prepro_verification_muta
        src/executables/prepro_gurobi_muta_similarity.cpp
        src/utils/GXLGraphReader.cpp
        src/utils/io.cpp
)

target_include_directories(prepro_verification_muta
        PRIVATE
        include
        $ENV{LIBLSAP_ROOT}/cpp/include
        $ENV{GEDLIB_ROOT}
        $ENV{GEDLIB_ROOT}/ext/boost.1.69.0
        $ENV{GEDLIB_ROOT}/ext/eigen.3.3.4/Eigen
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/src
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/ext/sgtelib/src
        $ENV{GEDLIB_ROOT}/ext/lsape.5/include
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/include
)
target_link_directories(prepro_verification_muta
        PRIVATE
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/lib
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/lib
        $ENV{GEDLIB_ROOT}/lib
)
target_link_libraries(prepro_verification_muta
        PRIVATE
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        libdoublefann.so.2
        libsvm.so
        libnomad.so
        libgxlgedlib.so
)

add_executable(lb_heuristics
        src/executables/lb_heuristics.cpp
        src/utils/GXLGraphReader.cpp
        src/utils/io.cpp
)

target_include_directories(lb_heuristics
        PRIVATE
        include
        $ENV{LIBLSAP_ROOT}/cpp/include
        $ENV{GEDLIB_ROOT}
        $ENV{GEDLIB_ROOT}/ext/boost.1.69.0
        $ENV{GEDLIB_ROOT}/ext/eigen.3.3.4/Eigen
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/src
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/ext/sgtelib/src
        $ENV{GEDLIB_ROOT}/ext/lsape.5/include
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/include
)
target_link_directories(lb_heuristics
        PRIVATE
        $ENV{GEDLIB_ROOT}/ext/nomad.3.8.1/lib
        $ENV{GEDLIB_ROOT}/ext/libsvm.3.22
        $ENV{GEDLIB_ROOT}/ext/fann.2.2.0/lib
        $ENV{GEDLIB_ROOT}/lib
)
target_link_libraries(lb_heuristics
        PRIVATE
        optimized ${GUROBI_CXX_LIBRARY}
        debug ${GUROBI_CXX_DEBUG_LIBRARY}
        ${GUROBI_LIBRARY}
        libdoublefann.so.2
        libsvm.so
        libnomad.so
        libgxlgedlib.so
)

set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)
install(TARGETS prepro_verification_aids
        RUNTIME
        DESTINATION bin)
install(TARGETS prepro_verification_muta
        RUNTIME
        DESTINATION bin)
install(TARGETS lb_heuristics
        RUNTIME
        DESTINATION bin)
# overwrite install() command with a dummy macro that is a nop
macro(install)
endmacro()

