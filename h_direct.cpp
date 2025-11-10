#include<iostream>
#include<cstdlib>
#include<cmath>
#include<iomanip>
#include<ctime>
#include<fstream>

using namespace std;

#define MCSTEPS 10000

double loi_uniforme(){
   return (rand()/(double)RAND_MAX);}

double f(double r, double alpha){
    return alpha*alpha*r*r*exp(-2*alpha*r);}

double density(double alpha){
    double rmax=100;
    
    label: double u=rmax*loi_uniforme();
    double w=sqrt(alpha)*loi_uniforme();
    
    if (w<=f(u,alpha)) return u;
    else goto label;}
    
double EL(double r, double a){
       return -(1/r)-a*(a-(2/r))/2.;}

int main(){


srand(time(0));
int i,k,mc;
double apas=0.05;
double a,xi,E_moy;

for(a=0.4;a<=1.5; a+=apas){

E_moy=0;

  for(mc=1;mc<=MCSTEPS;mc++){
 
    xi=density(a);

   E_moy+=EL(xi,a);}
       
   cout<<a<<"\t"<<E_moy/MCSTEPS<<endl;
}

return 0;}



