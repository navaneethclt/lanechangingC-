# Trajectory Optimization for lane changing with Ceres Solver

This project implements a trajectory optimization problem using the **Ceres Solver** to optimize a set of coefficients for vehicle trajectory based on acceleration, jerk, and time penalties. The code models the motion of a vehicle in both the `x` and `y` directions and performs optimization on the time parameter `T_opt`.

## Project Overview

The goal of this project is to minimize the acceleration, jerk, and time taken for a vehicle to travel a specified distance, given certain parameters like initial velocity, acceleration, and position. The optimization is formulated as a least-squares problem and solved using the **Ceres Solver**.

### Cost Function

The cost function is a sum of penalties for:
1. **Acceleration** (`lambda1`)
2. **Time taken** (`lambda2`)
3. **Jerk** (`lambda3`)

These penalties are used to compute the optimal trajectory by adjusting the time parameter `T_opt`.

### Key Components
- **CostFunction**: The main cost function that computes the residuals based on the provided parameters. It uses time-dependent equations to model the vehicle's position, velocity, acceleration, and jerk.
- **generateTrajectory**: After optimization, this function generates and saves the trajectory data (position, velocity, acceleration, and jerk) into a file named `trajectory.txt`.
- **Optimization**: The **Ceres Solver** is used to minimize the cost function with respect to `T_opt`, the optimal time for the trajectory.

## Installation

### Prerequisites

Before running the code, you will need to install the following dependencies:

1. **Ceres Solver**: For solving optimization problems.
   - Follow the installation instructions [here](https://ceres-solver.org/installation.html).
   
2. **Eigen**: A C++ library for linear algebra (already included in this repository as `eigen3/Eigen/Dense`).

3. **C++ Compiler**: Ensure you have a C++ compiler that supports C++11 or later.

4. **CMake**: For building the project.

### Steps to Build

   ```bash
   git clone https://github.com/your-username/trajectory-optimization.git
   cd trajectory-optimization
   cd build
   ./trajectory_optimization


