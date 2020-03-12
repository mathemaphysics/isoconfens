cmake_minimum_required(VERSION 3.0.0)
project(ICEPropensity VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(ICEPropensity ICEPropensity.cc)

find_package(Boost REQUIRED program_options)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    target_link_libraries(ICEPropensity ${Boost_LIBRARIES})
endif (Boost_FOUND)

find_library(GROMACS_LIBRARIES gromacs)
if(NOT ${GROMACS_LIBRARIES} EQUAL GROMACS_LIBRARIES-NOTFOUND)
    target_link_libraries(ICEPropensity ${GROMACS_LIBRARIES})
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
