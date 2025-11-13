//Diffusion Monte Carlo method for harmonic oscilator and morse potential sampled by simulation branching and diffusion processes

#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
#include<vector>

#define  PI 3.14159265359f

using namespace std;

double uniform_law(){
	return rand()/(double)RAND_MAX;}
	
double uniform_law_min_max(double min, double max){
	return (max-min)*uniform_law()+min;}

double normal_law(double mu, double sigma){
    double u=uniform_law();
    double v=uniform_law();
    return (sigma*sqrt(-2*log(u))*cos(2*PI*v)+mu);}
  

	
double potential(double x){
	//return 0.5*x*x; // harmonic
     return 0.5*(exp(-2*x)-2*exp(-x));//morse
}

void distribution(vector <double> &v, double *psi, double *x, const int n){
  double vmin, vmax, h;
  double a[n], b[n];
  int nbr[n];
  
  vmin=v[0];  
  vmax=v[0];
  
  for(int i=1; i<v.size(); i++){
    if(v[i]<vmin) vmin=v[i];
    if(v[i]>vmax) vmax=v[i];}
  
  h = (vmax-vmin)/n; 
  
  for(int i=0; i<n; i++){
    b[i] = vmin + h*i;
    a[i] = b[i]-h;
    nbr[i] = 0;
    for(int j=0; j<v.size(); j++)
      if((v[j]>=a[i]) && v[j]<b[i]) nbr[i]++;
  }
    
  for(int i=0; i<n; i++){
    x[i]= (a[i]+b[i])/2.; 
    psi[i] =  nbr[i]/(v.size()*h);}
    
}

int main(){
 
 
srand(time(0));
ofstream f;

f.open("morse.dat");

double dtau=0.005; 
double energy=0., energy_2=0.,  Et=0., n=20; 
double vmoy, alpha=1/dtau;
int Nt=100000, Niter=0, m, nNt, equilibrium=2.;    


vector<double>x0(Nt);
vector<double>x(Nt);
vector<double>W(Nt);
vector<double>temp;





// initial wave function
for( int i=0; i< Nt ; i++)
      x0[i] = uniform_law_min_max(-2,2); 


nNt=Nt;


for (double step = 0; step < n; step += dtau) { 
      
    temp.clear();
      
    for (int i = 0; i < nNt; i++) {
      
        x[i] = x0[i] + normal_law(0, sqrt(dtau)); 
        
        W[i] = exp( - dtau * (0.5*(potential(x0[i])+potential(x[i])) - Et) );
        
        m = floor( W[i] + uniform_law() );
     
        for(int k=0; k < m ; k++)
            temp.push_back(x[i]);
    }
    
    nNt=temp.size();
    x0.resize(nNt);
    x0=temp;
      
    
    vmoy=0.;
    for(int i=0; i< nNt; i++)
        vmoy += potential(x0[i])/nNt; 
        
    Et = vmoy + alpha*(1-(double)nNt/Nt);
    
    

    W.resize(nNt);
    x.resize(nNt);
  

  // energy update 
  if ( step > equilibrium ){
    energy += Et;
    energy_2 += vmoy;
    Niter++;
    
    }
     f << step    << "  " << Et   << "   " << vmoy <<endl;
     cout << step  << "  " <<  Et << "  " << energy/Niter  << "   " << energy_2/Niter <<endl;
  
  
}

f.close();

// plot histogram

const int nbr_boxes=100;
double psi[nbr_boxes], position[nbr_boxes];

f.open("histogram_m0.dat");

distribution(x0, psi, position, nbr_boxes);

for (int i=0; i<nbr_boxes; i++)
    f << position[i] << "  " << psi[i]<< endl;

f.close();
	
return 0;}

