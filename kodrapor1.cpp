
double fonction_n( double x, double y ){
// une fonciton qui donne les dependances de n spatiallement
}

void gradiant_fonction_n(double y, double& dndx, double& dndy ){
//la derivée de la fonciton qui donne les dependances de n spatiallement 

}


//je trace my ray avec un state vector
void fonction_deriv(int m , double s, double state_vec[4], double dy_ds[4] ){

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

double produit_deriv = ( ux* dndx + uy * dndy); 
dy_ds[2] = (1/n) * dndx - (ux / n) * produit_deriv;
dy_ds[3] = (1/n) * dndy - (uy / n) * produit_deriv;
}

