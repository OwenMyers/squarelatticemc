include_directories(include)
include_directories(${squarelatticemc_SOURCE_DIR}/include/squarelatticemc)

#mac
#set(BOOST_INCLUDE "/opt/local/include")
#centos
set(BOOST_INCLUDE "/usr/include")
#set(FFT_INCLUDE "/Users/owenmyers/TestFFTInstalLoc/include")
#set(FFT_INCLUDE "/usr/local/include")

find_path(FFTW_INCLUDES fftw3.h)
#include (FindPackageHandleStandardArgs)
#find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)
    
set(COVERAGE_LINK_FLAGS "-lfftw3")

#include_directories(${FFT_INCLUDE})
include_directories(${FFT_INCLUDES})
include_directories(${BOOST_INCLUDE})

set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}" )



# With templated classes the implimentation must be in the header file so this list does not contain
# all the estimator implementation files. Look in the corresponding include directory for everything
# templated.
SET(estimator_SRCS src/estimator.cpp
                   src/dimer_dimer_direction_vert.cpp
                   src/dimer_density.cpp
                   src/dimer_horz_struc_fac.cpp
                   src/z3_struc_fac.cpp
                   )

ADD_LIBRARY(estimators STATIC ${estimator_SRCS})

