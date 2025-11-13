// // Variational Monte Carlo for harmonic oscillator using Gaussian sampling to calculate the ground state energy


#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<ctime>

using namespace std;

#define PI 4.*atan(1.)
#define MCSTEPS 200000

double loi_uniforme(){
		return (rand()/(double)RAND_MAX);}

double loi_normale(double mu, double sigma){
	double u,v,x;
	u=loi_uniforme();
	v=loi_uniforme();

	x=sqrt(-2*log(u))*cos(2*PI*v);
	
		return (sigma*x+mu);}
	
double EL(double xi,double a){
	double s;

	s=(1/2.)*(a*a+(1-pow(a,4))*xi*xi);
		return s;}
		
int main(){

srand(time(0));
int i,mc;
double a,xi,E_moy;
double apas=0.05;

	for(a=0;a<=1.5;a+=apas){
	E_moy=0;
		for(mc=1;mc<=MCSTEPS;mc++){
		
		xi=loi_normale(0,1/(a*sqrt(2.)));
		
		E_moy+=EL(xi,a);}
		
	cout<<a<<"\t"<<E_moy/MCSTEPS<<endl; }
	
	return 0;}
		
