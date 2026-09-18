#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// Simple structure to represent a point inside the AdS3 bulk (Poincaré patch)
struct PointAdS {
    double x; // Boundary spatial coordinate
    double z; // Holographic direction (z -> 0 corresponds to the UV boundary)
};

class RyuTakayanagiSimulator {
private:
    double R_ads;        // Radius of the Anti-de Sitter space
    double G_newton;     // 3D Newton's gravitational constant
    double cutoff_uv;    // UV regulator cutoff to prevent divergence at z = 0
    double L_boundary;   // Length of the boundary subsystem interval
    int n_points;        // Number of points to discretize the curve

    std::vector<PointAdS> curve;

public:
    RyuTakayanagiSimulator(double r, double g, double eps, double l, int pts)
        : R_ads(r), G_newton(g), cutoff_uv(eps), L_boundary(l), n_points(pts) {
        initialize_straight_curve();
    }

    // Start with a naive straight-line configuration anchored at the UV cutoff
    void initialize_straight_curve() {
        curve.clear();
        double x_start = -L_boundary / 2.0;
        double dx = L_boundary / (n_points - 1);

        for (int i = 0; i < n_points; ++i) {
            curve.push_back({x_start + i * dx, cutoff_uv});
        }
    }

    // Calculates the total proper length using the hyperbolic AdS3 metric tensor
    double calculate_ads_length() {
        double total_length = 0.0;
        for (size_t i = 0; i < curve.size() - 1; ++i) {
            double dx = curve[i+1].x - curve[i].x;
            double dz = curve[i+1].z - curve[i].z;
            double z_avg = (curve[i].z + curve[i+1].z) / 2.0;

            // AdS3 Metric: ds^2 = (R^2 / z^2) * (dx^2 + dz^2)
            double ds = (R_ads / z_avg) * std::sqrt(dx * dx + dz * dz);
            total_length += ds;
        }
        return total_length;
    }

    // Relaxation algorithm using numerical gradient descent to find the minimal surface
    void optimize_minimal_surface(int steps, double learning_rate) {
        // Boundary endpoints remain strictly fixed
        for (int step = 0; step < steps; ++step) {
            std::vector<PointAdS> new_curve = curve;

            for (int i = 1; i < n_points - 1; ++i) {
                // Local numerical evaluation of geometric spatial forces (gradient)
                double original_z = curve[i].z;
                
                // Small variations to calculate the numerical derivative of length
                curve[i].z = original_z + 1e-5;
                double len_plus = calculate_ads_length();

                curve[i].z = original_z - 1e-5;
                double len_minus = calculate_ads_length();

                curve[i].z = original_z; // Reset

                double gradient_z = (len_plus - len_minus) / 2e-5;

                // Move interior coordinates deeper into the bulk to minimize path area
                new_curve[i].z = original_z - learning_rate * gradient_z;
                
                // Physical safeguard boundary check
                if (new_curve[i].z < cutoff_uv) {
                    new_curve[i].z = cutoff_uv;
                }
            }
            curve = new_curve;
        }
    }

    // Explicit execution of the Ryu-Takayanagi Formula: S = Area / (4 * G_N)
    double calculate_holographic_entropy() {
        double minimal_area = calculate_ads_length();
        return minimal_area / (4.0 * G_newton);
    }

    void display_profile() {
        std::cout << "\n--- Minimal Surface Profile (AdS Geodesic) ---" << std::endl;
        std::cout << std::setw(12) << "X" << std::setw(15) << "Z (Bulk)" << std::endl;
        // Step interval selection to prevent console output flooding
        int stride = n_points / 10;
        for (int i = 0; i < n_points; i += stride) {
            std::cout << std::setw(12) << curve[i].x << std::setw(15) << curve[i].z << std::endl;
        }
    }
};

int main() {
    std::cout << "Starting Ryu-Takayanagi minimal surface simulation..." << std::endl;

    // Physical parameter space variables
    double R = 1.0;          // AdS space radius curvature
    double G = 0.05;         // 3D bulk Newton gravitational constant
    double epsilon = 0.01;   // UV regularization bounds limit (cutoff)
    double L = 2.0;          // Subregion spatial interval length boundary
    int points = 100;        // Discretization array grid density resolution

    RyuTakayanagiSimulator sim(R, G, epsilon, L, points);

    std::cout << "Initial length (straight line): " << sim.calculate_ads_length() << " - Non-minimal baseline." << std::endl;

    // Execute iterative curve minimization execution path
    sim.optimize_minimal_surface(5000, 0.002);

    std::cout << "Final length (after relaxation): " << sim.calculate_ads_length() << std::endl;
    
    sim.display_profile();

    // Final boundary system entropy analysis evaluation output
    double entropy = sim.calculate_holographic_entropy();
    std::cout << "\n=======================================================" << std::endl;
    std::cout << " Holographic Entanglement Entropy (S_A): " << entropy << std::endl;
    std::cout << "=======================================================" << std::endl;

    return 0;
}
