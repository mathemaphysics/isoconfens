cmake_minimum_required(VERSION 3.0.0)
project(ICEPropensity VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(ICEPropensity ICEPropensity.cc)
add_executable(GMXTemplateICE GMXTemplateICE.cc)

find_package(Boost REQUIRED program_options filesystem regex)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    target_link_libraries(ICEPropensity ${Boost_LIBRARIES})
    target_link_libraries(GMXTemplateICE ${Boost_LIBRARIES})
endif (Boost_FOUND)

find_package(OpenMP)

find_package(GROMACS)
if(GROMACS_FOUND)
    message("Directory is ${GROMACS_INCLUDE_DIRS}")
    message("Library is ${GROMACS_LIBRARIES}")
    include_directories(${GROMACS_INCLUDE_DIRS})
    target_link_libraries(ICEPropensity libgromacs.so)
    target_link_libraries(GMXTemplateICE libgromacs.so)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
