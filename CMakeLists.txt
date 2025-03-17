cmake_minimum_required(VERSION 3.10)
project(trajectory_optimization)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find Ceres
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

# Add executable
add_executable(trajectory_optimization trajectory_optimization.cpp)

# Example explicit library linking. Replace with your actual paths.
link_directories("/usr/local/lib") # Replace with your Ceres library path
target_link_libraries(trajectory_optimization PUBLIC ceres glog pthread m)
