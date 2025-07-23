#include   <iostream>
#include    <cmath>
#include <vector>
#include <fstream>

using namespace std ;


int main(){                    

int xdim=200;
int ydim=200;
 int time_tot=350;
// Position of the source (center of the domain)
int xsource=100;
int ysource=100;
//Courant stability factor
int S=1/0.25;
//  Parameters of free space 
int epsilon0=(1/(36*3.1415))*1e-9;
int mu0=4*3.1415*1e-7;
int c=3e+8;
// Spatial grid step length 
int delta=1;
//Temporal grid step obtained using Courant condition
int deltat=1;
// Initialization of permittivity and permeability matrices
int epsilon=epsilon0;
int mu=mu0;

   
    ofstream out1{"Ez"};
    
vector <double> Ez(xdim*ydim);
vector <double> Hx(xdim*ydim);
vector <double> Hy(xdim*ydim);
for (int i=1;1<xdim;i++){
    for (int j=1; j<ydim ; j++){
        Ez[i,j]=0;
        Hy[i,j]=0;
        Hx[i,j]=0;
    }   
}
for (int i=0; i<200; i++){
        for (int j=0; j<200; j++){
            Ez[(xdim-1)+j*ydim]=Ez[0+j*ydim];
            Ez[0+j*ydim]=Ez[(xdim-1)+j*ydim];
            Ez[i+(ydim-1)*j]=Ez[i+0*j];
            Ez[i+0*j]=Ez[i+(xdim-1)*j];

            Hx[(xdim-1)+j*ydim]=Hx[0+j*ydim];
            Hx[0+j*xdim]=Hx[(xdim-1)+j*ydim];
            Hx[i+(xdim-1)*j]=Hx[i+0*j];
            Hx[i+0*j]=Hx[i+(ydim-1)*j];

            Hy[(xdim-1)+j*ydim]=Hy[0+j*ydim];
            Hy[0+j*ydim]=Hy[(xdim-1)+j*ydim];
            Hy[i+(xdim-1)*j]=Hy[i+0*j];
            Hy[i+0*j]=Hy[i+(ydim-1)*j];

          }
    }
// Update loop begins
for (int n=1; n<time_tot; n++) {
    int i;
    int j;
     
    if (n=1){
        for (int i=1;i<xdim ; i++){
            for (int j=1;j<ydim; j++ ){
                if (i=j=100){   
                    Ez[i,j]=1;
                }
            }
        }
       
    }
    

   
    for (int i=1;i<xdim;i++){
        for (int j=1;j<ydim;j++){
            //Vector update instead of for-loop for Hy and Hx fields
            Hx[i,j]= Hx[i,j] - (deltat/mu)*(Ez[i,j+1]-Ez[i,j])/delta; 
            Hy[i,j]= Hy[i,j] + (deltat/mu)*(Ez[i+1,j]-Ez[i,j])/delta;     
            //Vector update instead of for-loop for Ez field
            Ez[i,j]= Ez[i,j] + (deltat/epsilon)*((Hy[i,j]-Hy[i-1,j])/delta - ((Hx[i,j]-Hx[i,j-1])/delta));
        }
    }   
    out1 << n << " " << Ez[i,j] << endl;
}
}
