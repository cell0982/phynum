// you are my sunset and my sunrise 
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

// Constants
const double R_EARTH = 6371000.0; // Earth radius in meters
const double N0 = 1.000293;       // Refractive index at sea level
const double H_SCALE = 8500.0;    // Scale height (approx 8.5km)

// Atmospheric Index: n(y) = 1 + (n0-1) * exp(-alt / H)
double get_n(double x, double y) {
    double alt = std::sqrt(x*x + y*y) - R_EARTH;
    if (alt < 0) return N0;
    return 1.0 + (N0 - 1.0) * std::exp(-alt / H_SCALE);
}

// Gradient Calculation: dn/dx and dn/dy
void get_grad_n(double x, double y, double& dndx, double& dndy) {
    double r = std::sqrt(x*x + y*y);
    double alt = r - R_EARTH;
    
    // dn/d_alt = -(n0-1)/H * exp(-alt/H)
    double dn_dalt = -(N0 - 1.0) / H_SCALE * std::exp(-alt / H_SCALE);
    
    // Convert radial gradient to Cartesian components
    dndx = dn_dalt * (x / r);
    dndy = dn_dalt * (y / r);
}

// Ray Derivative Function (Hamiltonian Form)
void ray_deriv(int n_dims, double s, double state[4], double dst_ds[4]) {
    double x = state[0];
    double y = state[1];
    double ux = state[2];
    double uy = state[3];

    double n = get_n(x, y);
    double dndx, dndy;
    get_grad_n(x, y, dndx, dndy);

    dst_ds[0] = ux; 
    dst_ds[1] = uy; 

    double product = (ux * dndx + uy * dndy);
    // Ray equation: d(nu)/ds = grad(n) -> du/ds = (grad(n) - u(u.grad(n))) / n
    dst_ds[2] = (dndx - ux * product) / n;
    dst_ds[3] = (dndy - uy * product) / n;
}

// RK4 Solver
void rk4(int n, double& s, double state[], double ds) {
    double d1[4], d2[4], d3[4], d4[4], temp[4];
    
    ray_deriv(n, s, state, d1);
    for(int i=0; i<n; i++) temp[i] = state[i] + d1[i]*ds/2.0;
    
    ray_deriv(n, s+ds/2.0, temp, d2);
    for(int i=0; i<n; i++) temp[i] = state[i] + d2[i]*ds/2.0;
    
    ray_deriv(n, s+ds/2.0, temp, d3);
    for(int i=0; i<n; i++) temp[i] = state[i] + d3[i]*ds;
    
    ray_deriv(n, s+ds, temp, d4);
    for(int i=0; i<n; i++) {
        state[i] += (ds/6.0) * (d1[i] + 2.0*d2[i] + 2.0*d3[i] + d4[i]);
    }
    s += ds;
}

int main() {
    // 1. Initial State: Observer at (0, R_EARTH), looking HORIZONTAL (ux=1, uy=0)
    double state[4] = {0.0, R_EARTH, 1.0, 0.0}; 
    double s = 0;
    double ds = 500.0; // 500 meter steps
    
    std::ofstream outFile("sunset_data.csv");
    outFile << "x,y\n";

    // 2. Run simulation until ray leaves the dense atmosphere (~100km up)
    for(int i = 0; i < 2000; i++) {
        outFile << state[0] << "," << state[1] << "\n";
        
        rk4(4, s, state, ds);
        
        double alt = std::sqrt(state[0]*state[0] + state[1]*state[1]) - R_EARTH;
        if (alt > 100000.0) break; // Ray has reached space
    }
    outFile.close();

    // 3. Calculate Refraction Angle
    // Final ray direction (ux, uy)
    double final_angle_rad = std::atan2(state[3], state[2]);
    double refraction_deg = std::abs(final_angle_rad * 180.0 / M_PI);

    // 4. Calculate Time Delay
    // Earth rotates 360 degrees in 24 hours (1440 mins). 
    // Rotation speed = 0.25 degrees per minute.
    double extra_minutes = refraction_deg * 4.0;

    std::cout << "--- Sunset Simulation Results ---" << std::endl;
    std::cout << "Total Bending (Refraction): " << refraction_deg << " degrees" << std::endl;
    std::cout << "Extra Daylight Time: " << extra_minutes << " minutes" << std::endl;
    std::cout << "Data saved to sunset_data.csv" << std::endl;

    return 0;
}
