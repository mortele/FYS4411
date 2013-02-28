#include <armadillo>
#include "functions.h"


////calculates the wavefunction value at a point
//double psi(double *coordinates)
//{
//    r12 = sqrt((coordinates[0]-coordinates[3])*(coordinates[0]-coordinates[3])+(coordinates[1]-coordinates[4])*(coordinates[1]-coordinates[4])+(coordinates[2]-coordinates[5])*(coordinates[2]-coordinates[5]));
//    //psi = exp(-*alph(r1+r2)*exp((r12)/(2*(1+beta*r12))));
//    return exp(-alph(sqrt(coordinates[0]*coordinates[0]+coordinates[1]*coordinates[1]+coordinates[2]*coordinates[2])+sqrt(coordinates[3]*coordinates[3]+coordinates[4]*coordinates[4]+coordinates[5]*coordinates[5]))*exp((r12)/((1+beta*r12))));
//}



////calculates the local energy at a point
//double EL(mat &R,mat &r,
//          double psi)
//{
//    double r12 = R(0,1);
//    double r1 = R(0,0);
//    double r2 = R(1,1);
//    double E1 = (-Z*(1/r1 +1/r2) +1/r12);
//    double E2 = 0;

//    /*
//    for(int i=0; i< nDimensions*nParticles;i++){
//        coordinates[i]-=h;
//        psil = psi(coordinates);
//        coordinates[i]+=2*h;
//        psih = psi(coordinates);
//        coordinates[i]-=h;
//        E2-=derivative(psil, psi,psih);
//    }
//    */

//    vec oldR(nParticles-1);
//    double psil, psih;

//    for(int i = 0; i<nParticles;i++){
//        r1 = R(i,0);
//        for(int k=0; k<i; k++){
//            oldR(k) = R(i,k); //R is the matrix of distances
//        }
//        for(int k=i+1;k<nParticles;k++){
//            oldR(k) = R(i,k);
//        }
//        for(int j = 0; j<nDimesions;j++){
//            r(j,i)-=h; //r is the array of coordinates

//            updateforderivative(R,r, i);
//            psil = psi(R);

//            r(j,i)+=2*h;
//            updateforderivative(R,r, i);
//            psih = psi(R);
//            r(j,i)-=h;
//            E2-=computeDoubleDerivative(psil, psi,psih);
//        }
//    }

//    return E2/psi +E1;

//}


//void updateforderivative(mat &R, mat r, int i){
//    vec dx(nDimensions);
//    double dxx;

//    for(int k=0; k<i; k++) {
//        for(int l =0;l<nDimensions;l++){
//            dxx =r[l+i*nDimensions]-r[l+k*nDimensions];
//            dx(l) = dxx*dxx;
//        }

//        R(i,k) = sqrt(sum(dx)); //R is the matrix of distances
//    }
//    for(int k=i+1;k<nParticles;k++){
//        for(int l =0;l<nDimensions;l++){
//            dxx =r[l+i*nDimensions]-r[l+k*nDimensions];
//            dx(l) = dxx*dxx;
//        }

//        R(i,k) = sqrt(sum(dx)); //R is the matrix of distances
//    }
//    for(int l =0;l<nDimensions;l++){
//        dxx =r[l+i*nDimensions]*r[l+i*nDimensions];
//        dx(l) = dxx*dxx;
//    }
//    R(i,0) = sqrt(sum(dx));
//}



// calculates the double derivative
//double derivative(double psil, double psi,double psih) {
//    return (psil -2*psi + psih)/h2;

//}
