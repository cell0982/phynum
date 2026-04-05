#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>


using namespace std;

const double R_earth = 6371000;
const double N0 = 1.000321918;  // refractive index at 20 deg sea level
double H_scale = 50000; //scale height

//atmospheic index

double atm_index( double x, double y ){
    double alt = sqrt(x*x + y*y) - R_earth;
    if (alt < 0)
    return N0;
return 1 + (N0 - 1) * exp(- alt  / H_scale );
}


//gradiant n cette fois dndx nest pas 0

void grad_index(double  x, double y , double& dndx , double& dndy){
double r = sqrt(x*x + y*y);
double alt = r - R_earth;

// dn/d_alt = - ((n0 - 1) / hscale )* exp(...)

double dndalt = - ((N0 - 1)/ H_scale ) * exp(- alt / H_scale);


//conversion r derivative on x and y

dndx = dndalt * (x / r);
dndy = dndalt * (y / r);

}
// derive function
void deriv_sunset(int n_dimension, double s, double state_vec[4], double dst_ds[4] ){


double x = state_vec[0];
double y = state_vec[1];
double ux = state_vec[2];
double uy = state_vec[3];


double n = atm_index(x,y);
double dndx;
double dndy;
grad_index(x, y, dndx, dndy);

dst_ds[0] = ux; // dx/ds
dst_ds[1] = uy; //dy/ds

double produit_deriv = ( ux* dndx + uy * dndy); // definition sur ipad
dst_ds[2] = (1/n) * dndx - (ux / n) * produit_deriv;
dst_ds[3] = (1/n) * dndy - (uy / n) * produit_deriv;
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
deriv(n,x+ddx,yp,d3) ; /* 3eme estimat. (1/2 pas) */

for( i = 0; i< n; i++){ yp[i] = y[i] + d3[i]*dx ;}
deriv(n,x+dx,yp,d4) ;      /* 4eme estimat. (1 pas) */
/* estimation de y pour le pas suivant en utilisant
  une moyenne pondérée des dérivées en remarquant
  que : 1/6 + 1/3 + 1/3 + 1/6 = 1 */

for( i = 0; i < n ; i++)
 { y[i] = y[i] + dx*( d1[i] + 2*d2[i] + 2*d3[i] + d4[i] )/6 ; }

}



int main(){

int const n_dim = 4;
double s = 0;
double ds = 500 ; //  step around 500 meters
double state[4] = {0, R_earth, 1 , 0}; // probleme dunitairite

std::cout << "Starting simulation ..." << std::endl;
    std::ofstream outFile("sunset.csv");
    if(!outFile) return 1;

outFile << "x,y\n";
for(int i = 0 ; i < 2000; i++){
outFile << state[0]<< "," << state[1] << "\n";

rk4(n_dim, s, state, ds, deriv_sunset);

double alt = sqrt(state[0]*state[0] + state[1]*state[1]) - R_earth;
if( alt > 100000) break;
}
outFile.close();

double final_angle = acos(state[2]);
double refraction = abs( final_angle * 180 / M_PI);
double time_extra = refraction * 4;
cout << "Déviation totale des rayons du soleil : " << refraction << " deg " << endl;
cout << "Il vous reste encore " << time_extra << " minutes d'ensoleillement après le coucher du soleil." << endl;


return 0;
}







