#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "ceres/ceres.h"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

struct CostFunction {
    CostFunction(double u_x_in, double ax_in, double u_in, double v_in, double y0_in, double u_y_in, double ay_in, double W_in, double lambda1_in, double lambda2_in, double lambda3_in)
        : u_x(u_x_in), ax(ax_in), u(u_in), v(v_in), y0(y0_in), u_y(u_y_in), ay(ay_in), W(W_in), lambda1(lambda1_in), lambda2(lambda2_in), lambda3(lambda3_in) {}

    template <typename T>
    bool operator()(const T* T_opt, T* residual) const {
        // Define time and square of time
        T T_sq = T_opt[0] * T_opt[0];

        // Solve for coefficients a3, a4, a5, b3, b4, b5
        T a3 = (u - u_x - ax * T_opt[0]) / (T(3.0) * T_sq);
        T a4 = -a3 / (T(4.0) * T_opt[0]);
        T a5 = -a3 / (T(5.0) * T_sq);
        T b3 = (W - y0 - u_y * T_opt[0] - 0.5 * ay * T_sq) / (T(3.0) * T_sq);
        T b4 = -b3 / (T(4.0) * T_opt[0]);
        T b5 = -b3 / (T(5.0) * T_sq);

        // Calculate the position x(t) and y(t)
        T x_t = u_x * T_opt[0] + 0.5 * ax * T_sq + a3 * T_sq * T_opt[0] + a4 * T_sq * T_sq + a5 * T_sq * T_sq * T_opt[0];
        T y_t = y0 + u_y * T_opt[0] + 0.5 * ay * T_sq + b3 * T_sq * T_opt[0] + b4 * T_sq * T_sq + b5 * T_sq * T_sq * T_opt[0];

        // Calculate the velocity and acceleration
        T x_dot = u_x + ax * T_opt[0] + T(3.0) * a3 * T_sq + T(4.0) * a4 * T_sq * T_opt[0] + T(5.0) * a5 * T_sq * T_sq;
        T y_dot = u_y + ay * T_opt[0] + T(3.0) * b3 * T_sq + T(4.0) * b4 * T_sq * T_opt[0] + T(5.0) * b5 * T_sq * T_sq;

        // Calculate the acceleration
        T x_ddot = ax + T(6.0) * a3 * T_opt[0] + T(12.0) * a4 * T_sq + T(20.0) * a5 * T_sq * T_opt[0];
        T y_ddot = ay + T(6.0) * b3 * T_opt[0] + T(12.0) * b4 * T_sq + T(20.0) * b5 * T_sq * T_opt[0];

        // Calculate the jerk (rate of change of acceleration)
        T x_jerk = T(6.0) * a3 + T(24.0) * a4 * T_opt[0] + T(60.0) * a5 * T_sq;
        T y_jerk = T(6.0) * b3 + T(24.0) * b4 * T_opt[0] + T(60.0) * b5 * T_sq;

        // Calculate max acceleration and jerk using Euclidean norm
        T a_max = ceres::sqrt(x_ddot * x_ddot + y_ddot * y_ddot);
        T j_max = ceres::sqrt(x_jerk * x_jerk + y_jerk * y_jerk);

        // Cost function: penalties for acceleration, jerk, and time
        T cost_accel = lambda1 * a_max * a_max;
        T cost_time = lambda2 * T_opt[0] * T_opt[0];
        T cost_jerk = lambda3 * j_max * j_max;

        // Residual as the sum of the penalties
        residual[0] = cost_accel + cost_time + cost_jerk;
        return true;
    }

private:
    double u_x, ax, u, v, y0, u_y, ay, W, lambda1, lambda2, lambda3;
};
void generateTrajectory(double u_x, double ax, double u, double v, double y0, double u_y, double ay, double W, double T_opt, const string& filename) {
    ofstream outputFile(filename);
    if (outputFile.is_open()) {
        outputFile << fixed << setprecision(10);

        double T_sq = T_opt * T_opt;
        double a3 = (u - u_x - ax * T_opt) / (3.0 * T_sq);
        double a4 = -a3 / (4.0 * T_opt);
        double a5 = -a3 / (5.0 * T_sq);
        double b3 = (W - y0 - u_y * T_opt - 0.5 * ay * T_sq) / (3.0 * T_sq);
        double b4 = -b3 / (4.0 * T_opt);
        double b5 = -b3 / (5.0 * T_sq);

        int numPoints = 100;
        for (int i = 0; i <= numPoints; ++i) {
            double t = (double)i / numPoints * T_opt;
            double t_sq = t * t;

            double x_t = u_x * t + 0.5 * ax * t_sq + a3 * t_sq * t + a4 * t_sq * t_sq + a5 * t_sq * t_sq * t;
            double y_t = y0 + u_y * t + 0.5 * ay * t_sq + b3 * t_sq * t + b4 * t_sq * t_sq + b5 * t_sq * t_sq * t;

            double x_dot = u_x + ax * t + 3.0 * a3 * t_sq + 4.0 * a4 * t_sq * t + 5.0 * a5 * t_sq * t_sq;
            double y_dot = u_y + ay * t + 3.0 * b3 * t_sq + 4.0 * b4 * t_sq * t + 5.0 * b5 * t_sq * t_sq;

            double x_ddot = ax + 6.0 * a3 * t + 12.0 * a4 * t_sq + 20.0 * a5 * t_sq * t;
            double y_ddot = ay + 6.0 * b3 * t + 12.0 * b4 * t_sq + 20.0 * b5 * t_sq * t;

            outputFile << t << " " << x_t << " " << y_t << " " << x_dot << " " << y_dot << " " << x_ddot << " " << y_ddot << endl;
        }
        outputFile.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
int main() {
    // Input parameters
    double u_x = 1.0;
    double ax = 0.5;
    double u = 2.0;
    double v = 3.0;
    double y0 = 0.0;
    double u_y = 1.0;
    double ay = 0.5;
    double W = 5.0;
    double lambda1 = 1.0;
    double lambda2 = 1.0;
    double lambda3 = 0.1;

    // Initial guess for the optimization variable T_opt
    double T_opt = 1.0;

    // Create a Ceres problem
    ceres::Problem problem;

    // Add the cost function to the problem
    problem.AddResidualBlock(
        new ceres::AutoDiffCostFunction<CostFunction, 1, 1>(
            new CostFunction(u_x, ax, u, v, y0, u_y, ay, W, lambda1, lambda2, lambda3)),
        nullptr, &T_opt);

    // Configure the solver
    ceres::Solver::Options options;
    options.max_num_iterations = 100;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;

    // Solve the problem
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // Output the results
    std::cout << summary.BriefReport() << "\n";
    std::cout << "Optimal time: " << T_opt << "\n";

// Generate trajectory after optimization
    generateTrajectory(u_x, ax, u, v, y0, u_y, ay, W, T_opt, "trajectory.txt");

    std::cout << "Trajectory generated and saved to trajectory.txt" << std::endl;

    return 0;
}
