include_directories(include)
include_directories(${squarelatticemc_SOURCE_DIR}/include/squarelatticemc)

set(BOOST_INCLUDE "/opt/local/include")
set(FFT_INCLUDE "/Users/owenmyers/TestFFTInstalLoc/include")
set(COVERAGE_LINK_FLAGS "-lfftw3")
include_directories(${BOOST_INCLUDE})
include_directories(${FFT_INCLUDE})

set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}" )



# With templated classes the implimentation must be in the header file so this list does not contain
# all the estimator implementation files. Look in the corresponding include directory for everything
# templated.
SET(estimator_SRCS src/estimator.cpp
                   src/dimer_dimer_direction_vert.cpp
                   src/dimer_density.cpp
                   src/z3_struc_fac.cpp
                   src/dimer_horz_struc_fac.cpp)

ADD_LIBRARY(estimators STATIC ${estimator_SRCS})

