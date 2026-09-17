#include <iostream>
#include <cmath>
#include <functional>

// Function to evaluate: f(x) = e^x
double f(double x) {
    return std::exp(x);
}

// Function to calculate the left Riemann sum approximation of a definite integral
double leftRiemannSum(const std::function<double(double)>& func, double a, double b, int n) {
    double area = 0.0;
    double delta_x = (b - a) / n;

    // The loop iterates from 0 up to n-1 for the left endpoints
    for (int i = 0; i < n; ++i) {
        area += func(a + i * delta_x) * delta_x;
    }

    return area;
}

int main() {
    // Define the boundaries and number of subintervals
    double lower_bound = 0.0;
    double upper_bound = 1.0;
    int subintervals = 1000;

    std::cout << "Calculating Left Riemann Sum..." << std::endl;
    std::cout << "Interval: [" << lower_bound << ", " << upper_bound << "]" << std::endl;
    std::cout << "Subintervals (n): " << subintervals << std::endl;

    // Compute the approximation
    double result = leftRiemannSum(f, lower_bound, upper_bound, subintervals);

    std::cout << "Approximated Area: " << result << std::endl;

    return 0;
}
