// Variational Monte Carlo for hydrogen atom using Metropolis sampling to calculate the ground state energy and variance


#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>


using namespace std;


double WaveFunction(double xi, double a){
	return  exp(-a*xi);
}
 
 double LocalEnergy(double xi, double a){
	return -(1./xi)-a*(a-(2/xi))/2.;
}

double alea(){
	return (double)rand()/RAND_MAX;}
	
/* void condition(double &x0, double &xi){
    if (x0 < 0) x0 = jumpsize*(alea() .); 
    if (xi < 0) xi = jumpsize*(alea()); 
}*/

	
	int main(){
	
		ofstream f1;
		ofstream f2;
		
f1.open("H_metropolis_energy.data");
f2.open("H_metropolis_variance.data");
	
	
srand(time(0));
double energy, energy2, variance, error;
double x0, xi,  wfold, wfnew;
double alphaStep=0.05,jumpsize=1. ;
int mcCycles = 20000000;

for(double alpha=0.4; alpha<=1.5; alpha+=alphaStep){
        
        // initialization
	energy = 0.; energy2 = 0.;
	x0 = jumpsize * alea(); // random position
	wfold = WaveFunction(x0, alpha);
	
	for (int  mc=1; mc<=mcCycles; mc++){
	
	    // new trial positions/wave functions
	    xi = x0 + jumpsize * (alea() - 0.5);
	    if(xi>=0){
	    wfnew = WaveFunction(xi, alpha);

	    //Metropolis test
	    double ratio = xi*xi*wfnew*wfnew/(x0*x0*wfold*wfold);
	    double r = alea();
	    if( r<= ratio){
		x0= xi;
		wfold = wfnew;}} 
   
	
            energy += LocalEnergy(x0,alpha);
	    energy2 += LocalEnergy(x0,alpha) * LocalEnergy(x0,alpha);	
	}
	energy /= mcCycles;
	energy2 /= mcCycles;
	variance = energy2 - energy*energy;
	error = sqrt(variance/mcCycles);
	
	
	f1<< alpha << "\t" << energy << endl;
	f2<< alpha << "\t" << variance <<endl;
	
	cout <<"\t" << alpha <<"\t" << energy <<"\t" << variance <<"\t"<< error <<endl;
}



return 0;}
