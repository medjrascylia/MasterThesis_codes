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
	//return 0.5*(exp(-2*x)-2*exp(-x));//morse
	return -1/x;
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
    psi[i] +=  nbr[i]/(v.size()*h);}
    
}

int main(){
 
 
srand(time(0));
ofstream f;

f.open("div.dat");

double dtau= 0.001; 
double energy=0., energy_2=0.,  Et=0., n=25, n_eq=20; 
double vmoy, alpha=1/dtau;
int Nt=100000, Niter=0, m, nNt;    


vector<double>x0(Nt);
vector<double>x(Nt);
vector<double>W(Nt);

vector<double>tempr;
vector<double>tempx;
vector<double>tempy;
vector<double>tempz;


vector<double>r0(Nt);
vector<double>r(Nt);
vector<double>y0(Nt);
vector<double>z0(Nt);
vector<double>y(Nt);
vector<double>z(Nt);

const int nbr_boxes=100;
double psi[nbr_boxes]={0.}, position[nbr_boxes];

// initial wave function
for( int i=0; i< Nt ; i++){
      x0[i] = 0.;
      y0[i] = 0.;
      z0[i] = 1.;
      
      r0[i] = sqrt ( x0[i]*x0[i] + y0[i]*y0[i] + z0[i]*z0[i]);}

nNt=Nt;


for (double step = 0; step < n; step += dtau) { 
      
    tempx.clear();
    tempy.clear();
    tempz.clear();
    tempr.clear();
    
      
    for (int i = 0; i < nNt; i++) {
      
        x[i] = x0[i] + normal_law(0, sqrt(dtau)); 
        y[i] = y0[i] + normal_law(0, sqrt(dtau));
        z[i] = z0[i] + normal_law(0, sqrt(dtau));
        
        r[i] =  sqrt ( x[i]*x[i] + y[i]*y[i] + z[i]*z[i]); 
        
        W[i] = exp( - dtau * (0.5*(potential(r0[i])+potential(r[i])) - Et) );
        
        m = floor( W[i] + uniform_law() );
     
        for(int k=0; k < m ; k++){
            tempx.push_back(x[i]);
            tempy.push_back(y[i]);
            tempz.push_back(z[i]);
            tempr.push_back(r[i]);    
            }
    }
    
    nNt=tempx.size();
    x0.resize(nNt);
    y0.resize(nNt);
    z0.resize(nNt);
    r0.resize(nNt);
    
    x0=tempx;
    y0=tempy;
    z0=tempz;
    r0=tempr;
    
    vmoy=0.;
    for(int i=0; i< nNt; i++)
        vmoy += potential(r0[i])/nNt; 
        
    Et = vmoy + alpha*(1-(double)nNt/Nt);
    
    

    W.resize(nNt);
    x.resize(nNt);
  	y.resize(nNt);
  	z.resize(nNt);
  	r.resize(nNt);

  // energy update 
  if ( step > n_eq ){
    energy += Et;
    energy_2 += vmoy;
    Niter++;}
   


 cout << step  << "  "  << energy/nNt  << "   " << energy_2/nNt <<endl;


  //  f << step << "\t" << energy/Niter << "\t" << energy_2/Niter <<endl;
//  f << step << "\t" << Et << "\t" << vmoy << endl;

if(step>n_eq)
    distribution(r0, psi, position, nbr_boxes);

}

f.close();
  
  




f.open("histogram_h3.dat");



for (int i=0; i<nbr_boxes; i++)
    f << position[i] << "  " << psi[i]*dtau/(n-n_eq)<< endl;

f.close();


return 0;}


