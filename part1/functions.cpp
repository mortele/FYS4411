#include "functions.h"
#include <cmath>
using namespace std;


//calculates the wavefunction value at a point
double psi(double *coordinates)
{
    r12 = sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2));
    //psi = exp(-*alph(r1+r2)*exp((r12)/(2*(1+beta*r12))));
    return exp(-alph(sqrt(x*x+y*y+z*z)+sqrt(x2*x2+y2*y2+z2*z2))*exp((r12)/((1+beta*r12))));
}



//calculates the local energy at a point
double EL(double *coordinates,
          double psi)
{
    r12 = sqrt((coordinates[0]-coordinates[3])*(coordinates[0]-coordinates[3])+
            (coordinates[1]-coordinates[4])*(coordinates[1]-coordinates[4])+
            (coordinates[2]-coordinates[5])*(coordinates[2]-coordinates[5]));
    r1 = sqrt(coordinates[0]*coordinates[0] + coordinates[1]*coordinates[1] + coordinates[2]*coordinates[2]);
    r2 = sqrt(coordinates[3]*coordinates[3] + coordinates[4]*coordinates[4] + coordinates[5]*coordinates[5]);
    E1 = (-Z*(1/r1 +1/r2) +1/r12);
    E2 = 0;
    for(int i=0; i< nDimensions*nParticles;i++){
        coordinates[i]-=h;
        psil = psi(coordinates);
        coordinates[i]+=2*h;
        psih = psi(coordinates);
        coordinates[i]-=h;
        E2-=derivative(psil, psi,psih);
    }
    return E2/psi +E1;

}

double derivative(double psil, double psi,double psih) {
    return (psil -2*psi + psih)/h2;

}
