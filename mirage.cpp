#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <eigen3/Eigen/Dense>
#include <chrono>
#include <cmath>


using namespace std::chrono;
using namespace std;
using namespace Eigen;


//modele de n (mirage par example qui depend que de hauter et qui est constant par rapport a x cf notes)

//double something_real(double y, double T_y0){
//return T_y0 + 0.5 * y ;}


//mirage to simulate (mirage inferieur)
double mirage(double y){
double n0 =   1.000293; // air
double alpha = 0.0001; //increment of dependance of y in terms of refractive index
return n0 + alpha * y; //for superior mirage echqnge + with -
}

void grad_mirage(double y, double& dndx , double& dndy ){
double alpha = 0.2;
dndx = 0.0;
dndy = -alpha;

}

//fonction derivative f(X(s))

//je tracte my ray avec un state vector
void mirage_deriv(int m , double s, double state_vec[4], double dy_ds[4] ){


double x = state_vec[0];
double y = state_vec[1];
double ux = state_vec[2];
double uy = state_vec[3];


double n = mirage(y);
double dndx;
double dndy;
grad_mirage(y, dndx, dndy);

dy_ds[0] = ux; // dx/ds
dy_ds[1] = uy; //dy/ds

double produit_deriv = ( ux* dndx + uy * dndy); // definition sur ipad
dy_ds[2] = (1/n) * dndx - (ux / n) * produit_deriv;
dy_ds[3] = (1/n) * dndy - (uy / n) * produit_deriv;
}



void rk4(int n, double x, double y[] ,double dx,
               void deriv(int , double ,double [], double [])){
               //
/*-----------------------------------------
 sous programme de resolution d'equations
 differentielles du premier ordre par
 la methode de Runge-Kutta d'ordre 4
 x = abscisse
 y = valeurs des fonctions
 dx = pas
 n = nombre d'equations differentielles
 deriv = variable contenant le nom du
 sous-programme qui calcule les derivees
 ----------------------------------------*/
int i ;
double ddx ;
/* d1, d2, d3, d4 = estimations des derivees
   yp = estimations intermediaires des fonctions */
double d1[n], d2[n], d3[n], d4[n], yp[n];

ddx = dx/2;                /* demi-pas */

deriv(n,x,y,d1) ;          /* 1ere estimation */

for( i = 0; i< n; i++){ yp[i] = y[i] + d1[i]*ddx ; }
deriv(n,x+ddx,yp,d2) ;     /* 2eme estimat. (1/2 pas) */

for( i = 0; i < n; i++){ yp[i] = y[i] + d2[i]*ddx ; }
mirage_deriv(n,x+ddx,yp,d3) ; /* 3eme estimat. (1/2 pas) */

for( i = 0; i< n; i++){ yp[i] = y[i] + d3[i]*dx ;}
deriv(n,x+dx,yp,d4) ;      /* 4eme estimat. (1 pas) */
/* estimation de y pour le pas suivant en utilisant
  une moyenne pondérée des dérivées en remarquant
  que : 1/6 + 1/3 + 1/3 + 1/6 = 1 */

for( i = 0; i < n ; i++)
 { y[i] = y[i] + dx*( d1[i] + 2*d2[i] + 2*d3[i] + d4[i] )/6 ; }

}



int main(){

const int m = 4; // x, y , ux , uy
double s = 0;
double ds = 0.1 ; // taille step incriment
double state[4] = {0.0, 1.0, 0.7071, -0.7071};
    std::cout << "Starting simulation ..." << std::endl;
    std::ofstream outFile("mirage_data.csv");
    if(!outFile) return 1;

    outFile << "x,y\n";

    for(int i = 0; i < 100; i++) {


        outFile << state[0] << "," << state[1] << "\n";
      rk4(m, s, state, ds, mirage_deriv);
       //printf("iteration :%d",i);
      //fflush(stdout);

        s += ds;

        if((state[1]< 0) or! (state[1]) > 10) break; //si on est par terre ou trop haut
    }

    outFile.close();
    std::cout << "Simulation complete. Data saved to mirage_data.csv" << std::endl;
    return 0;
}
