#include <iostream>
#include <cmath>

// Define which mirage to simulate
enum MirageType { INFERIOR, SUPERIOR };
MirageType current_mode = INFERIOR; 

double get_m(double y) {
    double m0 = 1.00027; // Standard air index
    double alpha = 0.0001; 
    // Inferior: Index increases with height (m = m0 + alpha * y)
    // Superior: Index decreases with height (m = m0 - alpha * y)
    return (current_mode == INFERIOR) ? (m0 + alpha * y) : (m0 - alpha * y);
}

void get_grad_m(double y, double& dmdx, double& dmdy) {
    double alpha = 0.0001;
    dmdx = 0.0;
    dmdy = (current_mode == INFERIOR) ? alpha : -alpha; [cite: 13]
}

// The derivative function f(X(s)) [cite: 25, 26]
void mirage_deriv(int n, double s, double y_vec[], double dy_ds[]) {
    double y = y_vec[1];
    double ux = y_vec[2];
    double uy = y_vec[3];

    double m = get_m(y);
    double dmdx, dmdy;
    get_grad_m(y, dmdx, dmdy);

    dy_ds[0] = ux; // dx/ds 
    dy_ds[1] = uy; // dy/ds 

    double u_dot_gradm = (ux * dmdx + uy * dmdy);
    
    // du/ds calculation from the document 
    dy_ds[2] = (1.0 / m) * dmdx - (ux / m) * u_dot_gradm;
    dy_ds[3] = (1.0 / m) * dmdy - (uy / m) * u_dot_gradm;
}

// 3. How they appear "Closer" or "Farther"The distance at which the image appears depends on the curvature of the ray path.Inferior Mirage (The "Close" Mirage): Because the rays bend upward rapidly near the hot surface, the "virtual image" appears on the ground quite close to the observer. You see the sky reflected on the road, making the "water" look like it's just a few hundred meters away. Superior Mirage (The "Far" Mirage): These rays bend downward, following the curvature of the Earth. This allows you to see objects that are actually below the horizon (far away). The object appears to "loom" or float high in the sky. 4. Running the SimulationTo see the difference, you can initialize the state vector $X$  and run your rk4 loop twice:Set current_mode = INFERIOR: Start at $y=2$ and point the ray slightly down ($u_y = -0.05$). You will see it "bounce" off the air near $y=0$. Set current_mode = SUPERIOR: Start at $y=10$ and point the ray slightly up ($u_y = 0.05$). You will see it curve back toward the ground. Would you like me to add a "ground collision" detection to stop the simulation when the ray hits the $y=0$ surface?


#include <fstream>

int main() {
    const int n = 4; // x, y, ux, uy 
    double s = 0.0;
    double ds = 0.1; // Step size [cite: 24]
    double state[4] = {0.0, 2.0, 1.0, -0.02}; // Start at height 2, pointing slightly down

    std::ofstream outFile("mirage_data.csv");
    outFile << "x,y\n"; // Header

    for (int step = 0; step < 500; step++) {
        outFile << state[0] << "," << state[1] << "\n";
        
        // Use your RK4 to solve dx = f(x(s)) 
        rk4(n, s, state, ds, mirage_derivatives);
        s += ds;

        if (state[1] < 0 || state[1] > 10) break; // Stop if it hits ground or goes too high
    }
    outFile.close();
    return 0;
}
