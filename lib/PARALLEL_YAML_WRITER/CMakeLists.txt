cmake_minimum_required(VERSION 3.12)
project(parallel_yaml_writer)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Find Python
find_package(Python REQUIRED COMPONENTS Interpreter Development)

# Find pybind11 (adjust path if necessary)
find_package(pybind11 REQUIRED)

# Find yaml-cpp (adjust path if necessary)
find_package(yaml-cpp REQUIRED)

# Define your C++ module
pybind11_add_module(parallel_yaml_writer parallel_yaml_writer.cpp)

# Link yaml-cpp library
target_link_libraries(parallel_yaml_writer PRIVATE yaml-cpp)

# Include Python headers
target_include_directories(parallel_yaml_writer PRIVATE ${Python_INCLUDE_DIRS})

# Add pthread support
find_package(Threads REQUIRED)
target_link_libraries(parallel_yaml_writer PRIVATE Threads::Threads)

# Set C++ standard and pthread flag
set_target_properties(parallel_yaml_writer PROPERTIES 
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    COMPILE_OPTIONS "-pthread"
)

# Example installation path for the shared library
install(TARGETS parallel_yaml_writer LIBRARY DESTINATION /path/to/install/directory)

