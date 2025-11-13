// Variational Monte Carlo for Helium like atoms using 3 different trial wave functions with Metropolis sampling method to calculate ground energy, variance and error

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include<ctime>

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
      
/* computing the trail wave function */
double wf(double** r, int np, int d, double a, double b){
 
    // computing the norm square for each particle 
    double norm2_single_particle[np];
    double r_ij;
    for(int i=0;i<np;i++){
        norm2_single_particle[i]=0.;
        for(int j=0;j<d;j++){
          norm2_single_particle[i]+=r[i][j]*r[i][j];}}
     // computing rij 
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
    
       return exp(-a*arg);
  //   return exp(-a*arg)*(1+b*corr);
 //  return exp(-a*arg)*exp(0.5*corr/(1+b*corr));
    
}
 
    
/* computing the local energy using its analytical expression */
double local_energy(double **r, int np, int d, double a, double b, int z){
   
    double norm2_single_particle[np];
    double r_ij;
    double h=0.001, h2=1e-6;
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
     
int main(){
    srand(time(0));
    ofstream f;
    f.open("metropolis_be_wf1.data");
    
    const int nbr_particles = 2;
    const int dimension = 3;
    const int charge = 4;
    double **rold=NULL, **rnew=NULL;
    double alphaStep=.05, betaStep=0.05, owf, nwf, ratio;
    double energy, energy2, variance, error, deltaE, length=1.,beta;
    int mcCycles=2000000, equilibium=2000, accept;
     
    // allocate memory for matrices
    allocMatrix(&rold, nbr_particles, dimension);
    allocMatrix(&rnew, nbr_particles, dimension);
      

   // loop over the variational parameter alpha
  //  double alpha =1.82955, beta=0.401314;
  	

  for(double alpha=0.9; alpha<=4; alpha+=alphaStep) {
    for(double beta=0.1; beta<=1; beta+=betaStep) {
        energy=0.; accept=0;
      
       // random initial positions
      for(int i=0;i<nbr_particles;i++)
      
       for(int j=0;j<dimension;j++)
             rold[i][j] = length*(uniform_law()-0.5);
            
       // computing the wave function value at initial positions    
       owf=wf(rold,nbr_particles, dimension, alpha, beta);
     
       // loop over the Monte Carlo cycles
       for(int mc=1; mc<=(mcCycles+equilibium); mc++)  {      
     
       // computing the new positions 
       for(int i=0;i<nbr_particles;i++){
          for(int j=0;j<dimension;j++){
             rnew[i][j]=rold[i][j]+length*(uniform_law()-0.5); }}

       // computing the wave function value at the new positions   
       nwf=wf(rnew,nbr_particles, dimension, alpha, beta);
     
     
      // Metropolis test 
     
      ratio = (nwf*nwf)/(owf*owf);
     
      double  w=uniform_law();
     
      if (w<=ratio){
         for(int i=0;i<nbr_particles;i++){
            for(int j=0;j<dimension;j++){
              rold[i][j] =rnew[i][j];}}
         owf=nwf;
         // incrementing the acceptance value
         if(mc > equilibium)  accept++;
      }
    
      if(mc > equilibium){
         deltaE=local_energy(rold, nbr_particles, dimension,alpha, beta, charge);
         energy+=deltaE;
         energy2+=deltaE*deltaE;}
      } 
    
      energy /= mcCycles;
	energy2 /= mcCycles;
	variance = energy2 - energy*energy;
	error = sqrt(variance/mcCycles);
	
	
	cout << alpha << "\t"<< beta << "\t" << energy << "\t" << variance  << "\t" << error<< endl;
	f << alpha <<"\t" << beta << "\t" << energy << "\t" << variance  << endl;
 }
  }
   
     // delete allocate memory 
     deallocMatrix(rold,nbr_particles); 
     deallocMatrix(rnew,nbr_particles);
     
     f.close();
   
return 0;}
    
     
