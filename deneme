#include <iostream>
#include <cmath>
#include <vector>

// m(r) is a "donnée du problème" (data of the problem) [cite: 5, 7]
double get_m(double x, double y) {
    return 1.0 + 0.01 * x; // Example: index varying with x
}

// Gradient of the refractive index 
void get_grad_m(double x, double y, double& dmdx, double& dmdy) {
    dmdx = 0.01; // Example partial derivative
    dmdy = 0.0;
}

// This matches the signature required by your rk4.cpp
void ray_derivatives(int n, double s, double y_vec[], double dy_ds[]) {
    double x = y_vec[0];
    double y = y_vec[1];
    double ux = y_vec[2];
    double uy = y_vec[3];

    double m = get_m(x, y);
    double dmdx, dmdy;
    get_grad_m(x, y, dmdx, dmdy);

    // From the document: d(r)/ds = u [cite: 9, 13]
    dy_ds[0] = ux; 
    dy_ds[1] = uy;

    // From the document: du/ds = (1/m) * grad(m) - [u . grad(m)] * (u/m) [cite: 13]
    double u_dot_gradm = (ux * dmdx + uy * dmdy);
    
    dy_ds[2] = (1.0 / m) * dmdx - (ux / m) * u_dot_gradm;
    dy_ds[3] = (1.0 / m) * dmdy - (uy / m) * u_dot_gradm;
}

int main() {
    const int n = 4;        // x, y, ux, uy 
    double s = 0.0;         // current path length
    double ds = 0.1;        // step size (pas)
    double state[4] = {0.0, 0.0, 1.0, 0.0}; // Initial position and direction

    for (int step = 0; step < 100; step++) {
        // Call your uploaded rk4 function
        rk4(n, s, state, ds, ray_derivatives);
        
        s += ds;
        std::cout << "s: " << s << " x: " << state[0] << " y: " << state[1] << std::endl;
    }

    return 0;
}
