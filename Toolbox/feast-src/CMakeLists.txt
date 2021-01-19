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
endif()

if(BUILD_SHARED)
    if(USE_MKL)
        link_directories(${MKLROOT}/lib/intel64)
        if(MKLROOT MATCHES "(oneapi|oneAPI)")
            if((COMPILER_IDENTIFIER MATCHES "gcc-linux") OR (COMPILER_IDENTIFIER MATCHES "gcc-mac"))
                find_library(IOMPPATH iomp5 ${MKLROOT}/../../compiler/latest/linux/compiler/lib/intel64_lin)
                get_filename_component(IOMPPATH ${IOMPPATH} DIRECTORY)
                link_directories(${IOMPPATH})
                link_libraries(mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 m)
            elseif(COMPILER_IDENTIFIER MATCHES "gcc-win")
                link_libraries(mkl_rt)
            elseif(COMPILER_IDENTIFIER MATCHES "vs")
                find_library(IOMPPATH libiomp5md ${MKLROOT}/../../compiler/latest/windows/compiler/lib/intel64_win)
                get_filename_component(IOMPPATH ${IOMPPATH} DIRECTORY)
                link_directories(${IOMPPATH})
                set(DSUFIX "_dll")
                link_libraries(mkl_intel_lp64${DSUFIX} mkl_intel_thread${DSUFIX} mkl_core${DSUFIX} libiomp5md)
            endif()
        else()
            if((COMPILER_IDENTIFIER MATCHES "gcc-linux") OR (COMPILER_IDENTIFIER MATCHES "gcc-mac"))
                find_library(IOMPPATH iomp5 ${MKLROOT}/../lib/intel64)
                get_filename_component(IOMPPATH ${IOMPPATH} DIRECTORY)
                link_directories(${IOMPPATH})
                link_libraries(mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 m)
            elseif(COMPILER_IDENTIFIER MATCHES "gcc-win")
                link_libraries(mkl_rt)
            elseif(COMPILER_IDENTIFIER MATCHES "vs")
                find_library(IOMPPATH libiomp5md ${MKLROOT}/../compiler/lib/intel64)
                get_filename_component(IOMPPATH ${IOMPPATH} DIRECTORY)
                link_directories(${IOMPPATH})
                set(DSUFIX "_dll")
                link_libraries(mkl_intel_lp64${DSUFIX} mkl_intel_thread${DSUFIX} mkl_core${DSUFIX} libiomp5md)
            endif()
        endif()
    else()
        if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
            link_libraries(pthread gfortran quadmath openblas)
        else() # MSVC COMPILER
            link_libraries(libopenblas)
        endif()
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