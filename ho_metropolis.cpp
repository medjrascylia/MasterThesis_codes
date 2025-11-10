#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <fstream>

using namespace std;


double WaveFunction(double r, double alpha){
	return  exp(-0.5*alpha*alpha*r*r);
}

double LocalEnergy(double r, double alpha){
	return 0.5*r*r*(1-pow(alpha,4))+0.5*alpha*alpha;
}

double alea(){
	return (double)rand()/RAND_MAX;
}

double exact (double alpha){
    return 0.5 * alpha * alpha + (1 - pow(alpha, 4)) / (4 * alpha * alpha);
}




int main(){


double energy, energy2, variance, error;
double positionOld, positionNew,  wfold, wfnew;
double alphaStep=0.05, jumpSize=1.;
int mcCycles = 1000000;

ofstream outfile;

outfile.open("ho.data");

srand(time(0));


for(double alpha=0.4; alpha<=2; alpha+=alphaStep){
        
        // initialization
	energy = 0.; energy2 = 0.;
	positionOld = jumpSize * (alea() - 0.5); // random position
	//wfold = WaveFunction(positionOld, alpha);
	
	for (int  mc=1; mc<=mcCycles; mc++){
	wfold = WaveFunction(positionOld, alpha);
	    // new trial positions/wave functions
	    positionNew = positionOld + jumpSize * (alea() - 0.5);
	    wfnew = WaveFunction(positionNew, alpha);

	    //Metropolis test
	    double ratio = wfnew*wfnew/(wfold*wfold);
	    double r = alea();
	    if( r<= ratio){
		positionOld = positionNew;
		wfold = wfnew;} 
		
	    energy += LocalEnergy(positionOld,alpha);
	    energy2 += LocalEnergy(positionOld,alpha) * LocalEnergy(positionOld,alpha);	
	}
	
	energy /= mcCycles;
	energy2 /= mcCycles;
	variance = energy2 - energy*energy;
	error = sqrt(variance/mcCycles);
	outfile << alpha << "\t" << energy << "\t" << variance  << "\t" << error << endl;
	cout << alpha << "\t" << energy << "\t" << variance  << "\t" << error << endl;
	


	
}

outfile.close();

outfile.open("exact_ho.dat");

for(double alpha=0.4; alpha<=2; alpha+=alphaStep){
outfile << alpha <<"\t"  <<exact(alpha)<<endl;
}
outfile.close();

return 0;}

