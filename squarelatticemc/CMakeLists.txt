include_directories(${squarelatticemc_SOURCE_DIR}/estimatorlib/include)
include_directories(${squarelatticemc_SOURCE_DIR}/include/squarelatticemc)

#mac
#set(BOOST_INCLUDE "/opt/local/include")
#centos
set(BOOST_INCLUDE "/usr/include")
include_directories(${BOOST_INCLUDE})


find_path (FFTW_INCLUDES fftw3.h)
find_library (FFTW_LIBRARIES NAMES fftw3)

include_directories(${FFT_INCLUDES})


set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}" )

#file(GLOB SOURCES loops.cpp ${qdpm_SOURCE_DIR}/src/qdpm/*.cpp)

add_executable(squarelatticemc square_lattice_mc.cpp 
                    ${squarelatticemc_SOURCE_DIR}/src/squarelatticemc/lattice.cpp
                    ${squarelatticemc_SOURCE_DIR}/src/squarelatticemc/model.cpp
                    ${squarelatticemc_SOURCE_DIR}/src/squarelatticemc/random.cpp
                    ${squarelatticemc_SOURCE_DIR}/src/squarelatticemc/vertex.cpp
                    )

target_link_libraries(squarelatticemc estimators ${FFTW_LIBRARIES})



