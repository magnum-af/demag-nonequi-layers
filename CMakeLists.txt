project(demag_nonequidistant)
cmake_minimum_required(VERSION 3.15)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # for IDE support

# set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wpedantic")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion")

find_package(PythonLibs REQUIRED)
find_package(Boost COMPONENTS python numpy REQUIRED)

# Without this, any build libraries automatically have names "lib{x}.so"
set(CMAKE_SHARED_MODULE_PREFIX "")

add_library(demag_nonequidistant MODULE demag_nonequidistant/demag_nonequidistant.cpp)

target_link_libraries(demag_nonequidistant ${Boost_LIBRARIES})
target_include_directories(demag_nonequidistant PRIVATE ${PYTHON_INCLUDE_DIRS})

include(GNUInstallDirs)
install(TARGETS demag_nonequidistant LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR})