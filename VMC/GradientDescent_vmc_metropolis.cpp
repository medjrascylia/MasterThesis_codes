// Variational Monte Carlo sampled using Metropolis algorithm for hellium like atoms optimized using gradient descent method 



#include<iostream>
#include<iomanip>
#include<ctime>
#include<cmath>

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
													
    
double uniform_law(){
	return rand()/(double)RAND_MAX;}
	  
double wf(double** r,int np, int d, double a,double b){
          /* computing the first term */
  double norm2_single_particle[np];
  double arg=0.;
for(int i=0; i<np ; i++){
norm2_single_particle[i]=0;
	for(int j=0; j<d ; j++){
	norm2_single_particle[i] += r[i][j]*r[i][j];}
	arg += sqrt(norm2_single_particle[i]);}
		
	/*computing the corr term */
 double corr=0., rij;	
  for(int i=0; i <np-1; i++){
  	for(int j=0; j<np ; j++){
  	rij=0;
  		for(int k=0; k<d ; k++){
  		rij += ( r[i][k] - r[j][k])*( r[i][k] - r[j][k]);}
  		corr =+ sqrt(rij);}}    
  		
  		return exp(-a*arg); 
  	   //  return exp(-a*arg)*(1+ b*corr);
  		//return exp(-a*arg)*exp(0.5*corr/(1+b*corr));
}
  
													
double local_energy(double **r, int np, int d, double a, double b, double z){
   double norm2_single_particle[np],  rij;

    for(int i=0;i<np;i++){
        norm2_single_particle[i]=0.;
        for(int j=0;j<d;j++){
          norm2_single_particle[i]+=r[i][j]*r[i][j];
          }}
       
       double ec = 0, temp;
       double h=0.0001, h2 = 1e-8;
       double wfold, wfplus, wfminus;
       
        for(int i=0;i<np;i++){
            for(int j=0;j<d;j++){
            temp = r[i][j];
            r[i][j]=temp +h;
            wfplus = wf(r,np,d,a,b);
            r[i][j]= temp-h;                        
            wfminus = wf(r,np,d,a,b);
            r[i][j] = temp;       
            wfold = wf(r,np,d,a,b);
            ec -= 0.5*(wfplus-2*wfold+wfminus)/h2/wfold;
            }}
       
       double ep = 0.;
       
       for(int i=0;i<np;i++){
       		  ep -= z*(1./sqrt(norm2_single_particle[i]));}
 
 		double ep_ee = 0.;
 		for(int i=0;i<np-1;i++){
 			for(int j=i+1;j<np;j++){
 				rij = 0.;
 				for(int k=0;k<d;k++){
 					rij += 	(r[i][k]-r[j][k])*(r[i][k]-r[j][k]);}    	           
                ep_ee += 1/sqrt(rij);}}
 
         return ec+ep+ep_ee;} 		



void wfderivative(double** r,int np,int d, double a, double b, double* dwf){
	double wfaplus,wfaminus,wfbplus,wfbminus;
	double h=0.0001;
      
	wfaplus=wf(r,np,d,a+h,b);
	wfaminus=wf(r,np,d,a-h,b); 
	dwf[0]= (log(wfaplus)-log(wfaminus))/(2*h);  //derivative w respect to alpha
	
	wfbplus=wf(r,np,d,a,b+h);	
	wfbminus=wf(r,np,d,a,b-h);
	dwf[1]= (log(wfbplus)-log(wfbminus))/(2*h); //derivative w respect to beta 	
	}
	
	
	
void metropolis(int np, int d, int z, double length, double &energy,double &variance, int mcCycles, int equilibrium,
	double a, double b, double *denergy){
	double  dwf[2],dPsi[2]={0,0}, dPsiE[2]={0,0};

    double **rold=NULL, **rnew=NULL;
	double owf, nwf, ratio, deltaE, energy2;
// allocate memory for matrices
allocMatrix(&rold, np, d);
allocMatrix(&rnew, np, d);
energy=0;
energy2=0.;

	// random initial positions
		for(int i=0;i<np;i++)
		for(int j=0;j<d;j++)
			rold[i][j] = length*(uniform_law()-0.5);
			
	// computing the wave function value at initial positions
			owf=wf(rold, np, d, a, b);
			
		// loop over the Monte Carlo cycles
		
	for(int mc=1; mc<=(mcCycles+equilibrium); mc++) {
	
	// computing the new positions
		for(int i=0;i<np;i++){
		for(int j=0;j<d;j++){
	rnew[i][j]=rold[i][j]+length*(uniform_law()-0.5); }}
	
		// computing the wave function value at the new positions
		nwf=wf(rnew,np, d, a, b);
	// Metropolis test
		ratio = (nwf*nwf)/(owf*owf);
		double w=uniform_law();
	if (w<=ratio){
		for(int i=0;i<np;i++){
		for(int j=0;j<d;j++){
			rold[i][j] =rnew[i][j];}}
			owf=nwf;}
	if(mc > equilibrium){
		deltaE=local_energy(rold, np, d, a, b, z);
	 
    	energy+=deltaE;
		energy2+=deltaE*deltaE;
   
		wfderivative (rold, np, d, a, b,dwf);

    for(int i=0; i<2; i++){
         dPsi [i] += dwf[i];
         dPsiE [i] += dwf[i]*deltaE;}
			
	}}
	
	energy2 /=mcCycles; 
	energy /= mcCycles;   
	variance = energy2-energy*energy;       

     for(int i=0; i<2; i++)
          denergy[i] = 2*(dPsiE[i] - dPsi[i]*energy)/mcCycles;

     
// delete allocate memory
deallocMatrix(rold,np);
deallocMatrix(rnew,np);

}




void gradientDescent(double &a, double &b, double a0, double b0, double precision, double learningRate, int &iter, 
                     int np, int d, int z, double length, int mcCycles, int equilibrium, double &energy, 
                     double &variance, double*denergy){
	
	double a0=0.9; 
	double b0=0;
        
        
        metropolis(np,d,z,length,energy,variance,mcCycles,equilibrium,a0,b0,denergy); 
        
        a=a0 - learningRate*denergy[0];      
        b=b0 - learningRate*denergy[1];
  
   
    while(!(abs(a-a0) < precision && abs(b-b0)<precision)){
        
	metropolis(np,d,z,length,energy,variance,mcCycles,equilibrium,a,b,denergy);
	
	a0=a;
	a=a0 - learningRate*denergy[0];

	b0=b;
	b=b0 - learningRate*denergy[1] ;
	iter++; 
    }
} 

int main(){
  srand(time(0));
  const int np = 2;
  const int d=3;
  const int z = 3;
  
  
 
  double energy,variance, denergy[2];
 
 double a,b, a0=0.9, b0=0.1, precision=1e-3, length=0.9, eta=0.1; //eta is the learningRate
 
 int mcCycles=5000, equilibrium =1000, iter=1;

  gradientDescent(a, b, a0 , b0, precision,  eta, iter, np,  d,  z,length, mcCycles, equilibrium, energy, variance, denergy);
  
  
  
   cout << "Optimal variational parameters : " << endl;
   cout << "alpha = " << a << "   " << "beta = " << b <<  "   "  << "energy [eV] = " << energy*27.2 << "   "  << "# of iterations = " << iter << endl;
   cout << "gradient_EL_alpha = " <<  denergy[0] << "   " << "gradient_EL_beta = "  <<  denergy[1] <<   endl;
   
  
  return 0 ;}
  
  
