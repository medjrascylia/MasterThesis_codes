// Variational Monte Carlo for Helium-like atoms using Metropolis-hastings algorithm for ground energy and variance computing


#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<ctime>
#include<fstream>
#define PI 4.*atan(1.)


using namespace std;
 
/* function to allocate memory for matrices */  
void allocMatrix(double ***matrix, int nrow, int ncol){
    *matrix = (double **)malloc(sizeof(int *) * nrow);
    for(int i = 0; i < ncol; i++)
        *(*matrix + i) = (double *)malloc(sizeof(int) * ncol);}

/* function to free allocate memory */
void deallocMatrix(double **matrix, int nrow){
    for(int i = 0; i < nrow; i++) free(matrix[i]);
    free(matrix);} 
 
/* uniform law */ 
double uniform_law(){
    return (double)rand()/RAND_MAX;} 
 
/* normal law */
double normal_law(double mu,double sigma ){
    double u=uniform_law();
    double v=uniform_law();
    return (sigma*sqrt(-2*log(u))*cos(2*PI*v)+mu);}
  
/* computing the trial wave function */
double wf(double** r, int np, int d, double a){
    // computing the norm squared for each particle 
    double norm2_single_particle[np];
    for(int i=0;i<np;i++){
        norm2_single_particle[i]=0.;
        for(int j=0;j<d;j++){
          norm2_single_particle[i]+=r[i][j]*r[i][j];}}
   // computing the argument under the exp       
   double arg=0.;
   for(int i=0;i<np;i++)
      arg+=norm2_single_particle[i]; 
   return exp(-a*sqrt(arg));
}
 

// computing the components of the quantum force  
// using its  analytical express
void quantumf_components(double **r, double **qf, int np, int d, double a){
    double norm2_single_particle[np];
    for(int i=0;i<np;i++){
        norm2_single_particle[i]=0.;
        for(int j=0;j<d;j++){
          norm2_single_particle[i]+=r[i][j]*r[i][j];}}
     // computing the x,y,z force components for each particle 
     for(int i=0;i<np;i++){
        for(int j=0;j<d;j++){
          qf[i][j]=-2*a*r[i][j]/sqrt(norm2_single_particle[i]);}}
 }
 
// computing the components of the quantum force  
// using numerical derivative
void qfc_numerical(double **r, double **qf, int np, int d, double a){
    double temp0, temp1, temp2;
    double wf1, wf2;
    double h = 0.001;
    for(int i=0;i<np;i++){
        for(int j=0;j<d;j++){
          temp0=r[i][j];
          temp1=temp0-h;
          r[i][j]=temp1;
          wf1=wf(r,np,d,a);
          temp2=temp0+h;
          r[i][j]=temp2;
          wf2=wf(r,np,d,a);
          qf[i][j]=(temp2-temp1)/(2*h);
          r[i][j]=temp0;}}
}
    
// computing the local energy using its
// analytical expression
double local_energy(double **r, int np, int d, double a){
    double norm2_single_particle[np],  arg=0.;
    for(int i=0;i<np;i++){
        norm2_single_particle[i]=0.;
        for(int j=0;j<d;j++){
          norm2_single_particle[i]+=r[i][j]*r[i][j];}
         arg+=norm2_single_particle[i];}
    return -1/sqrt(arg)-.5*a*(a-2/sqrt(arg));}
       

// computing the exp of the ratio between  green functions 
double green_function(double **r, double **r0,double **qf, int np, int d, double a, double t){
    double sum1[np],  sum2[np];
    double temp1, temp2; 
    // computing the force components at position r0
    quantumf_components(r0,qf,np,d,a);
    
    for(int i=0;i<np;i++){
        sum1[i]=0.; temp1=0.;
        for(int j=0;j<d;j++){
             temp1 = r[i][j]-r0[i][j]-0.5*qf[i][j]*t;
             sum1[i]+= temp1*temp1;}}   //somme de chaque particules
    temp1=0.;
    for(int i=0;i<np;i++)
        temp1+=sum1[i];
         
    // computing the force components at position r
    quantumf_components(r,qf,np,d,a);
    
    for(int i=0;i<np;i++){
        sum2[i]=0.; temp2=0.;
        for(int j=0;j<d;j++){
             temp2 = r0[i][j]-r[i][j]-0.5*qf[i][j]*t;
             sum2[i]+= temp2*temp2;}}
    temp2=0.;
    for(int i=0;i<np;i++)
        temp2+=sum2[i];
        
     return  exp(-0.5*(temp2-temp1)/t);
}    
   
// the main program   
int main(){

ofstream f1;
ofstream f2;

f1.open("H_hastings_energy.data");
f2.open("H_hastings_variance.data");

    const int nbr_particles = 2;
    const int dimension = 3;
    double **rold=NULL, **rnew=NULL, **qforce=NULL;
    double alphaStep=.05, owf, nwf, ratio, green_ratio,variance,energy2,deltaE2,error;
    double energy, deltaE, deltaT=.005;
    int mcCycles=10000000, thermalization=1000, accept;
     
    srand(time(0));
    // allocate matrices for old positions, new positions, and quantum force components
    allocMatrix(&rold, nbr_particles, dimension);
    allocMatrix(&rnew, nbr_particles, dimension);
    allocMatrix(&qforce, nbr_particles, dimension);
      
    for(double alpha=.4; alpha<=1.5; alpha+=alphaStep){
      energy=0.;  accept=0;
      
      for(int i=0;i<nbr_particles;i++)
        for(int j=0;j<dimension;j++)
            rold[i][j] = uniform_law()-0.5;
    
      owf=wf(rold,nbr_particles, dimension, alpha);
     
      for(int mc=1; mc<=(mcCycles+thermalization); mc++)  {      
     
      // computing the quantum force components
      quantumf_components(rold, qforce, nbr_particles, dimension, alpha);
      
      //computing the new positions
      for(int i=0;i<nbr_particles;i++){
        for(int j=0;j<dimension;j++){
           rnew[i][j]=rold[i][j]+0.5*qforce[i][j]*deltaT+normal_law(0,1)*sqrt(deltaT);}}
   
      //computing the new  wave function 
      nwf=wf(rnew,nbr_particles, dimension, alpha);
     
      // Metropolis test
     
      green_ratio=green_function(rnew,rold,qforce,nbr_particles,dimension,alpha,deltaT);  
      ratio=green_ratio*(nwf*nwf)/(owf*owf);
     
      double  w=uniform_law();
     
      if (w<=ratio){
        for(int i=0;i<nbr_particles;i++){
           for(int j=0;j<dimension;j++){
                rold[i][j] =rnew[i][j];}}
        owf=nwf;
        if(mc > thermalization)  accept++;
      }

    
      if(mc > thermalization){
        deltaE = local_energy(rold, nbr_particles, dimension,alpha);
        deltaE2 = local_energy(rold, nbr_particles, dimension,alpha)*local_energy(rold, nbr_particles, dimension,alpha);
        energy+=deltaE;
        energy2 +=deltaE2;
        }} 
      
      energy=energy/mcCycles;
      energy2=energy2/mcCycles;
      variance = energy2 - energy*energy;
      error = sqrt(variance/mcCycles);
      
   //f1 << alpha << "\t" << energy << endl;
   //f2 << alpha << "\t" << variance << endl;
   
    cout << alpha << "\t" << energy<<"\t" << variance <<"\t" << error << endl;
   
  }
   
  // delete allocate memeory
  deallocMatrix(rold,nbr_particles); 
  deallocMatrix(rnew,nbr_particles);
  deallocMatrix(qforce,nbr_particles);
  
return 0;}
