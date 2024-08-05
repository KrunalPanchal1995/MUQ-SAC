cmake_minimum_required(VERSION 3.12)
project(yamlwriter)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Find Python
find_package(Python REQUIRED COMPONENTS Interpreter Development)
include_directories(${Python_INCLUDE_DIRS})

# Find pybind11 (adjust path if necessary)
find_package(pybind11 REQUIRED)
include_directories(${pybind11_INCLUDE_DIRS})

# Find yaml-cpp (adjust path if necessary)
find_package(yaml-cpp REQUIRED)
include_directories(${yaml-cpp_INCLUDE_DIRS})

# Define your C++ module
pybind11_add_module(yamlwriter yamlwriter.cpp)

# Link yaml-cpp library
target_link_libraries(yamlwriter PRIVATE yaml-cpp)

# Example installation path for the shared library
install(TARGETS yamlwriter LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

