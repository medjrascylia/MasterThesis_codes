// Variational Monte Carlo using Metropolis sampling for helium like atoms with 3 different trial wave functions optimized via the swarm particle optimization (PSO) method

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <fstream>
#include<sstream>

#define  PI 3.14159265359f
using namespace std;

// functions prototypes and declarations 
void   allocMatrix(double***, int, int); 
void   deallocMatrix(double**, int);
double uniform_law();
double uniform_law_min_max(double, double);
double normal_law(double, double);
double wf(double**, int , int, double, double);
void   wfderivative(double*, double**, int, int,  double, double);
double local_energy(double**, int, int, double, double, int);
void   initialization(double*, double*, double*, double&, double*, double*, double*, double&, int, int, int, int);
void   best_particle_position(double*, double*, double*, double*, int, int, int, int);
void   best_swarm_position(double*, double&, double*, double&, int, int, int, int);
double optimization(double, double, int, int, int);
void   metropolis(int, int, int, double, double&,  double&, int, int, double, double);

int main(){


    srand(time(0));
    const int nbr_particles = 2;
    const int dimension = 3;
    const int charge = 4;
    const int np = 20; // swarm number of particles 
 
    double alpha[np], valpha[np], beta[np], vbeta[np], best_alpha[np], best_beta[np];
    double best_sw_alpha, best_sw_beta;
    double omega=0.25, c1=0.5, c2=1.; // control parameters
    int niter=50; 
    double energy,variance;
 
int i;

//ofstream f;
//f.open("swarm_test.data");

  
    initialization(alpha, valpha, best_alpha, best_sw_alpha, beta, vbeta, best_beta, best_sw_beta, np, nbr_particles, dimension, charge);
   
   
   for(int iter=1; iter<=niter; iter++){
      
         ostringstream oss;
         oss << "swarm" << iter << ".data";  
         string filename = oss.str();      
         ofstream file(filename);
      
      for(i=0; i<np; i++){
            valpha[i]=omega*valpha[i]+c1*uniform_law()*(best_alpha[i]-alpha[i])+c2*uniform_law()*(best_sw_alpha-alpha[i]);
            alpha[i]= alpha[i]+valpha[i];
            vbeta[i]=omega*vbeta[i]+c1*uniform_law()*(best_beta[i]-beta[i])+c2*uniform_law()*(best_sw_beta-beta[i]);
            beta[i]= beta[i]+vbeta[i];
            
          //  metropolis(2, 3 , 2,  0.9, energy,  variance, 10000, 2000, alpha[i], beta[i]);
            
        //    cout << alpha [i] <<"\t" << beta [i] << "\t" << energy <<endl;
          //  f << alpha [i] << "\t" << beta [i] << "\t" << energy << endl;
                           
            file << alpha[i] << "\t" << beta[i] << endl;
           }
            
            
            
        // metropolis(2, 3 , 2,  0.9, energy,  variance, 10000, 2000, alpha[i], beta[i]);
       
      
    
       best_particle_position(alpha, best_alpha,  beta, best_beta, np, nbr_particles, dimension, charge);
       best_swarm_position(best_alpha,  best_sw_alpha, best_beta, best_sw_beta, np, nbr_particles, dimension, charge);

           
           
         energy= optimization(best_sw_alpha, best_sw_beta,nbr_particles,dimension,charge);
       
       
       
        cout <<  best_sw_alpha << "\t" << best_sw_beta <<"\t" << energy <<endl;
        f << best_sw_alpha << "\t" << best_sw_beta <<"\t" << energy<<endl;
         
      }
     
     
    
     
      double alpha_opt=best_sw_alpha, beta_opt=best_sw_beta; 
      int mcCycles = 2000000, equilibrium = 2000;
      double length=.25;
   //   double energy, variance, length=.25;
     
    // metropolis(nbr_particles, dimension, charge,  length, energy,  variance, mcCycles, equilibrium, alpha_opt, beta_opt);
     
     cout << endl;
     cout << "Energy of the ground state : " << endl;
     cout <<  "energy [Ha] = " << energy << "   "  << "energy [eV] = " << energy*27.2 << "   " << "variance = " << variance << endl;
    
return 0;}
 
// Here are the definitions of the above different declared functions 
// computing the trail wave function 
double wf(double** r, int np, int d, double a, double b){
 
    // computing the norm square for each particle 
    double norm2_single_particle[np];
    double r_ij;
    for(int i=0;i<np;i++){
        norm2_single_particle[i]=0.;
        for(int j=0;j<d;j++){
          norm2_single_particle[i]+=r[i][j]*r[i][j];}}
     // computing r_ij 
     double corr=0.;
     for(int i=0;i<np-1;i++){
        for(int j=i+1;j<np;j++){
          r_ij=0.;
          for(int k=0;k<d;k++){
             r_ij+=(r[i][k]-r[j][k])*(r[i][k]-r[j][k]);}
             corr+=sqrt(r_ij);
          }}
    // computing the argument under the exp       
     double arg=0.;
     for(int i=0;i<np;i++)
         arg+=sqrt(norm2_single_particle[i]); 
    
     
       //   return exp(-a*arg);
        //  return exp(-a*arg)*(1+b*corr);
       return exp(-a*arg)*exp(0.5*corr/(1+b*corr)); 
     

  
}

// computing the wave function derivatives with respect to variational parameters  
void wfderivative(double *dwf, double **r, int np, int d,  double a, double b){
    double h=1e-4;
    dwf[0]=0.5*(wf(r,np,d,a+h,b)-wf(r,np,d,a-h,b))/h;
    dwf[1]=0.5*(wf(r,np,d,a,b+h)-wf(r,np,d,a,b-h))/h;
}
    
// computing the local energy
double local_energy(double **r, int np, int d, double a, double b, int z){
    double norm2_single_particle[np];
    double r_ij;
    double h=1e-4, h2=1e-8;
    double temp, wfminus, wfplus, wfold; 
 
   
    for(int i=0;i<np;i++){
        norm2_single_particle[i]=0.;
        for(int j=0;j<d;j++){
          norm2_single_particle[i]+=r[i][j]*r[i][j];}}
   
    double ec=0.;
    for(int i=0;i<np;i++){
        for(int j=0;j<d;j++){
         temp = r[i][j];
         r[i][j] = temp+h;
         wfplus = wf(r,np,d,a,b);
         r[i][j]= temp-h;
         wfminus = wf(r,np,d,a,b);
         r[i][j]=temp;
         wfold = wf(r,np,d,a,b);
         ec -= 0.5*(wfplus+wfminus-2*wfold)/wfold/h2;}}
      
      double ep=0.;
      for(int i=0;i<np;i++){
          ep-=z*(1./sqrt(norm2_single_particle[i]));}
      
        double ep_ee=0.;
        for(int i=0;i<np-1;i++){
         for(int j=i+1;j<np;j++){
          r_ij=0.;
          for(int k=0;k<d;k++){
             r_ij+=(r[i][k]-r[j][k])*(r[i][k]-r[j][k]);}
             ep_ee+=1./sqrt(r_ij);
             }}
             
    return ec+ep+ep_ee;}
     


// initialization
void initialization(double *alpha, double *valpha, double *best_alpha, double &best_sw_alpha, 
                    double *beta, double *vbeta, double *best_beta, double &best_sw_beta, int np, int nbr, int d, int z){
      for(int i=0; i<np; i++){
          alpha[i]=  uniform_law_min_max( 0.9 , 3. ); //  0 < alpha < Z   
          valpha[i]= uniform_law_min_max( 0.9 , 3. ); 
          beta[i]=  uniform_law_min_max( 0 , 1. );
          vbeta[i]= uniform_law_min_max( 0 , 1. );}
          
          
      for(int i=0; i<np; i++){
              best_alpha[i]=alpha[i];
              best_beta[i]=beta[i];}
              
     
       best_swarm_position(best_alpha, best_sw_alpha, best_beta, best_sw_beta, np, nbr, d, z);  
}   

// best particle position
void best_particle_position(double *alpha,  double *best_alpha, double *beta, double *best_beta, int np, int nbr, int d, int z){
   
      for(int i=1; i<np; i++)
        if (optimization(alpha[i],beta[i],nbr,d,z)<optimization(best_alpha[i],best_beta[i],nbr,d,z)){
              best_alpha[i]=alpha[i];
              best_beta[i]=beta[i];}
}

    
// best swarm position
void best_swarm_position(double *best_alpha, double &best_sw_alpha, double *best_beta,  double &best_sw_beta,  int np, int nbr, int d, int z){
     
      best_sw_alpha=best_alpha[0];
      best_sw_beta=best_beta[0];
      
      for(int i=1; i<np; i++)
          if (optimization(best_alpha[i],best_beta[i],nbr,d,z)<optimization(best_sw_alpha,best_sw_beta,nbr,d,z)){
                  best_sw_alpha=best_alpha[i];
                  best_sw_beta=best_beta[i];}
}


double optimization(double a, double b, int nbr, int d, int z){
    double **rold=NULL, **rnew=NULL;
    double  owf, nwf, ratio;
    double energy, energy2, variance, deltaE, length=0.9;
    int mcCycles=10000, equilibrium=1000;
     
    // allocate memory for matrices
    allocMatrix(&rold, nbr, d);
    allocMatrix(&rnew, nbr, d);
      

        energy=0.; energy2=0.; 
      
       // random initial positions
      for(int i=0;i<nbr;i++)
       for(int j=0;j<d;j++)
             rold[i][j] = length*(uniform_law()-0.5);
            
       // computing the wave function value at initial positions    
       owf=wf(rold,nbr, d, a, b);
     
       // loop over the Monte Carlo cycles
       for(int mc=1; mc<=(mcCycles+equilibrium); mc++)  {      
     
       // computing the new positions 
       for(int i=0;i<nbr;i++){
          for(int j=0;j<d;j++){
             rnew[i][j]=rold[i][j]+length*(uniform_law()-0.5); }}

       // computing the wave function value at the new positions   
       nwf=wf(rnew,nbr, d, a, b);
     
     
      // Metropolis test 
     
      ratio = (nwf*nwf)/(owf*owf);
     
      double  w=uniform_law();
     
      if (w<=ratio){
         for(int i=0;i<nbr;i++){
            for(int j=0;j<d;j++){
              rold[i][j] =rnew[i][j];}}
         owf=nwf;
      }
    
      if(mc > equilibrium){
         deltaE=local_energy(rold, nbr, d, a, b, z);
         energy+=deltaE;
         energy2+=deltaE*deltaE;}
      } 
    
      energy /= mcCycles;
      energy2 /=mcCycles;
      variance = energy2-energy*energy;	
	

   
     // delete allocate memory 
     deallocMatrix(rold,nbr); 
     deallocMatrix(rnew,nbr);
      
  
    
    // return variance;
	return energy;    
      
}


// the metropolis algorithm 
void metropolis(int np, int d, int z,  double length, double &energy, double &variance, int mcCycles, int equilibrium, double a, double b){
    double **rold=NULL, **rnew=NULL;
    double  owf, nwf, ratio, deltaE, energy2;    
    
    // allocate memory for matrices
    allocMatrix(&rold, np, d);
    allocMatrix(&rnew, np, d);
    
     energy=0; energy2=0.;
    // random initial positions
      for(int i=0;i<np;i++)
       for(int j=0;j<d;j++)
             rold[i][j] = length*(uniform_law()-0.5);
            
       // computing the wave function value at initial positions    
       owf=wf(rold, np, d, a, b);
     
       // loop over the Monte Carlo cycles
       for(int mc=1; mc<=(mcCycles+equilibrium); mc++)  {      
     
       // computing the new positions 
       for(int i=0;i<np;i++){
          for(int j=0;j<d;j++){
             rnew[i][j]=rold[i][j]+length*(uniform_law()-0.5); }}

       // computing the wave function value at the new positions   
       nwf=wf(rnew,np, d, a, b);
     
     
      // Metropolis test 
     
      ratio = (nwf*nwf)/(owf*owf);
      double  w=uniform_law();
     
      if (w<=ratio){
         for(int i=0;i<np;i++){
            for(int j=0;j<d;j++){
              rold[i][j] =rnew[i][j];}}
              owf=nwf;} 
      if(mc > equilibrium){
      deltaE=local_energy(rold, np, d, a, b, z);
      energy+=deltaE;
      energy2+=deltaE*deltaE;
      }}
      
     energy /= mcCycles; energy2 /=mcCycles;
     variance = energy2-energy*energy;
     
    
     // delete allocate memory 
     deallocMatrix(rold,np); 
     deallocMatrix(rnew,np);
}

 
// uniform laws 
double uniform_law(){
    return (double)rand()/RAND_MAX;} 
    
double uniform_law_min_max(double min, double max){
    return (max-min)*uniform_law()+min;} 
 
// normal law 
double normal_law(double mu,double sigma ){
    double u=uniform_law(), v=uniform_law();
    return (sigma*sqrt(-2*log(u))*cos(2*PI*v)+mu);}

// function to allocate memory for matrices   
void allocMatrix(double ***matrix, int nrow, int ncol){
    *matrix = (double **)malloc(sizeof(int *) * nrow);
    for(int i = 0; i < ncol; i++)
        *(*matrix + i) = (double *)malloc(sizeof(int) * ncol);}

// function to free allocate memory 
void deallocMatrix(double **matrix, int nrow){
    for(int i = 0; i < nrow; i++) free(matrix[i]);
    free(matrix);} 
