// fibre en 3d pour la visualisation

// fibre optique 2d

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>


using namespace std;

double r_0 = 50 ; //micrometre
double n_0 = 1.454;
double n_core = 1.484 ; //prise de RP Photonics

double delta_n =  (- n_0 + n_core) / n_core;

double indice_optique(double x,  double y){
double r = sqrt(y*y + x*x) ;

return n_0 +  delta_n*(1.0 - r*r / (r_0*r_0) );
}

void grad_index(double  x, double y , double z, double& dndx , double& dndy, double& dndz){
double r = sqrt(y*y + x*x);

double dndr =  - (delta_n/ r_0*r_0)*2 * r;

//conversion r deriv en x et y

dndy = dndr * (y / r);
dndx = dndr * (x / r);
dndz = 0;
}


void deriv_fibre(int n_dimension, double s, double state_vec[4], double dst_ds[4] ){


double x = state_vec[0];
double y = state_vec[1];
double z = state_vec[2];
double ux = state_vec[3];
double uy = state_vec[4];
double uz = state_vec[5];


double n = indice_optique(x,y);
double dndx;
double dndy;
double dndz;
grad_index(x, y, z, dndx, dndy, dndz);

dst_ds[0] = ux; // dx/ds
dst_ds[1] = uy; //dy/ds
dst_ds[2] = uz;

double produit_deriv = ( ux* dndx + uy * dndy + uz * dndz); // definition sur ipad
dst_ds[3] = (1/n) * dndx - (ux / n) * produit_deriv;
dst_ds[4] = (1/n) * dndy - (uy / n) * produit_deriv;
dst_ds[5] = (1/n) * dndz - (uz / n) * produit_deriv;
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

int const n_dim = 6;
double s = 0;
double ds = 0.001; //
double theta = 10;
double state[6] = {0.0, 0.1, 0.0, sin(M_PI/ 180 * theta), 0, cos(M_PI/ 180 * theta)}; // probleme dunitairite

std::cout << "Starting simulation ..." << std::endl;
    std::ofstream outFile("fibre3d.csv");
    if(!outFile) return 1;

outFile << "x,y,z\n";
for(int i = 0 ; i < 500000; i++){
outFile << state[0]<< "," << state[1] << "," << state[2]<< "\n";

rk4(n_dim, s, state, ds, deriv_fibre);

if((state[0]) > 2 * 10e6) break;

}
outFile.close();

double conservation = state[5]*state[5] + state[4]*state[4] + state[3]*state[3];

cout << "unitairy sum " << conservation <<endl;
cout << "...simulation done look at python  " << endl;


return 0;
}

