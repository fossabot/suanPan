cmake_minimum_required(VERSION 3.0)

project(feast Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB FEAST_SRC *.f90)

set(SPIKE_SRC
    spike-src/cpivul.f
    spike-src/cspike_smp.f90
    spike-src/dpivul.f
    spike-src/dspike_smp.f90
    spike-src/dzlbprim.f90
    spike-src/f90_functions_wrapper.f90
    spike-src/sclbprim.f90
    spike-src/spike_smp_utilities.f90
    spike-src/spivul.f
    spike-src/sspike_smp.f90
    spike-src/zpivul.f
    spike-src/zspike_smp.f90
)

if(USE_MKL)
    add_definitions(-DMKL)
    include_directories(${MKLROOT}/mkl/include)
    if((COMPILER_IDENTIFIER MATCHES "gcc-linux") OR (COMPILER_IDENTIFIER MATCHES "gcc-mac"))
        link_directories(${MKLROOT}/mkl/lib/intel64)
        link_directories(${MKLROOT}/lib/intel64)
        link_libraries(mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 m)
    elseif(COMPILER_IDENTIFIER MATCHES "gcc-win")
        # directly link dll files
        link_directories(${MKLROOT}/redist/intel64/mkl)
        link_libraries(mkl_rt)
    elseif(COMPILER_IDENTIFIER MATCHES "vs")
        link_directories(${MKLROOT}/mkl/lib/intel64)
        link_directories(${MKLROOT}/compiler/lib/intel64)
        set(DSUFIX "_dll")
        link_libraries(mkl_core${DSUFIX} mkl_intel_lp64${DSUFIX} mkl_intel_thread${DSUFIX} libiomp5md)
    endif()
else()
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
        link_libraries(pthread gfortran quadmath openblas)
    else() # MSVC COMPILER
        link_libraries(libopenblas)
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffixed-line-length-none")
endif()

if(BUILD_SHARED)
    add_library(${PROJECT_NAME} SHARED ${FEAST_SRC} ${SPIKE_SRC})
else()
    add_library(${PROJECT_NAME} STATIC ${FEAST_SRC} ${SPIKE_SRC})
endif()

message("FEAST Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")