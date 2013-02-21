#include <armadillo>
#include <time.h>
#include<iomanip>

#include "variationalmc.h"
#include "lib.cpp"


using namespace std;
using namespace arma;

/* VMC constructor. */
VariationalMC::VariationalMC() :
    nParticles  (2),
    nDimensions (3),
    nCycles     (10000000),
    N       (nCycles / 10),
    idum    (17),
    charge  (2.0),
    h       (0.0001),
    h2      (h * h),
    alph    (1.0),
    beta    (1.0),
    Z       (2.0),
    stepSize(0.01),
    D       (1.0) {
}

/* Runs the Metropolis algorithm nCycles times. */
double VariationalMC::runMetropolis(double alpha, double beta) {
    this->alph = alpha;
    this->beta = beta;

    mat coordinatesNew  = zeros<mat>(nParticles, nDimensions);
    mat coordinatesOld  = zeros<mat>(nParticles, nDimensions);
    mat Rnew            = zeros<mat>(nParticles, nParticles);   // Matrix of distances and magnitudes.
    mat Rold            = zeros<mat>(nParticles, nParticles);
    vec quantumForceOld(nParticles, nParticles);
    vec quantumForceNew(nParticles, nParticles);

    double ecoeff          = 0.0;
    double newWaveFunction = 0.0;
    double oldWaveFunction = 0.0;

    double energy          = 0.0;
    double energy2         = 0.0;

    double energySum       = 0.0;
    double energy2Sum      = 0.0;

    double greensFunction  = 0.0;

    int    accepted        = 0;
    // int    percent         = 0;

    double randI;
    int    iRand;

    // Fill coordinates arrays with random values.
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            coordinatesNew(i,j) = coordinatesOld(i,j) = (ran0(&idum)-0.5) / (0.5*alph);
        }
    }

    // Updates distances and magnitudes in this initial state.
    updateRmatrix(coordinatesNew, Rnew);
    updateRmatrix(coordinatesOld, Rold);

    // Compute the wave function in this initial state.
    oldWaveFunction = computePsi(Rnew);

    // Metropolis loop.
    for (int k = 0; k < nCycles; k++) {
//        if (k % (nCycles / 10) == 0) {
//            cout << percent << " %" << endl;
//            percent += 10;
//        }

        // Suggest new positions for all particles, i.e. new state.

        randI = ran0(&idum) * nParticles;
        iRand = floor(randI);

        // Trekk fra gauss-random her !!!!!!!!!!!!
        for (int j = 0; j < nDimensions; j++) {
            coordinatesNew(iRand,j) += (ran0(&idum)-0.5) * stepSize;
        }


        // Compute the wavefunction in this new state.
        updateRmatrix(coordinatesNew, Rnew);
        newWaveFunction = computePsi(Rnew);

        // Regn ut quantumForce her !!!!!!!!!!!!!

        // Compute the inside of the exponential term of the difference between Greens functions.
        for (int i = 0; i < nParticles; i++) {
            for (int j = 0; j < nDimensions; j++) {
                greensFunction += 0.5 * (quantumForceOld(nDimensions*i+j) + quantumForceNew(nDimensions*i+j)) *
                                  (D * h * 0.5 * (quantumForceOld(nDimensions*i+j) - quantumForceNew(nDimensions*i+j)) -
                                   coordinatesNew(i,j) + coordinatesOld(i,j));
            }
        }

        // Compute the fraction GreensF(new) / GreensF(old).
        greensFunction = exp(greensFunction);

        // Check if the suggested move is accepted.
        ecoeff = greensFunction * newWaveFunction * newWaveFunction / (oldWaveFunction * oldWaveFunction);
        if (ecoeff > ran0(&idum)) {
            accepted++;
            coordinatesOld = coordinatesNew;

            // Energy changes from previous state.
            energy = computeEnergy(Rnew, coordinatesNew, newWaveFunction);
//            if (Rnew(0,1) < 0.01) {
//                cout << "hei" << endl;
//            }

        } else {
            coordinatesNew = coordinatesOld;

            // Energy remains unchanged.
        }

        // Add energy of this state to the energy sum.

        if (k == N) {
            energySum  = 0.0;
            energy2Sum = 0.0;
        }

        energySum  += energy;
        energy2Sum += energy * energy;
    }

    vec quantumForce(nParticles * nDimensions);
    quantumForce = computeQuantumForce(Rnew ,coordinatesNew, newWaveFunction);

    // Calculate the expected value of the energy, the energy squared, and the variance.
    energy  = energySum  / (nCycles * 0.9);
    energy2 = energy2Sum / (nCycles * 0.9);

    cout << "<E>  = " << energy << endl;
    cout << "<E²> = " << energy2 << endl;
    cout << "Variance  = " << energy2 - energy*energy << endl;
    cout << "Std. dev. = " << sqrt(energy2 - energy*energy) << endl;
    cout << "Accepted steps / total steps = " << ((double) accepted) / nCycles << endl;

    return energy;
}


/* Computes the wavefunction in a state defined by position matrix r. */
double VariationalMC::computePsi(const mat &R) {

    double returnVal = 0.0;
    returnVal = exp(-alph * (R(0,0) + R(1,1))); // * exp(R(0,1) / (2 * (1  + beta * R(0,1))));

    //    for (int i = 0; i < nParticles; i++) {
    //        for (int j = (i + 1); j < nParticles; j++) {
    //            returnVal += (1 / R(i,j)) * exp(-alph * (R(i,i) + R(j,j))) * (exp(R(i,j)) / (1 + beta * R(i,j)));
    //        }
    //    }
    return returnVal;
}




/* Computes the local energy of a state defined by position matrix r, and distance matrix R.
 * EL = 1/psi * H * psi */
double VariationalMC::computeEnergy(mat &R, mat &r, double psi)
{



    double r12 = R(0,1);
    //double r1  = R(0,0);
    //double r2  = R(1,1);
    double E2  = 0; "hei morten !!!"
    double E1  = (1/r12); // + (-Z*(1/r1 +1/r2));    // this is the commutative part of the
                                             // hamiltonian (we can just multiply it with psi.)


    /*
                              for(int i=0; i< nDimensions*nParticles;i++){
                                  coordinates[i]-=h;
                                  psil = psi(coordinates);
                                  coordinates[i]+=2*h;
                                  psih = psi(coordinates);
                                  coordinates[i]-=h;
                                  E2-=derivative(psil, psi,psih);
                              }
                              */

    //vec oldR(nParticles);
    double psil, psih;

    for(int i = 0; i<nParticles;i++){
        /*        
        r1 = R(i,0);

        for(int k=0; k<i; k++){
            oldR(k) = R(i,k); //R is the matrix of distances
        }


        for(int k = i + 1; k < nParticles; k++) { // nParticles
            oldR(k) = R(i,k);
        }
        */

        for(int j = 0; j<nDimensions;j++){
            r(i,j)-=h; //r is the array of coordinates

            updateForDerivative(R,r, i);
            psil = computePsi(R);

            r(i,j)+=2*h;
            updateForDerivative(R, r, i);
            psih = computePsi(R);
            r(i,j)-=h;
            E2-=computeDoubleDerivative(psil, psi,psih);
            updateForDerivative(R, r, i);   // set all values back to normal
        }
    }
    //cout << E2 / psi << " " << E1 << endl;
    return E1; // + E2 / (2 * psi);

}

/* Computes a numerical approximation to the double derivative of psi. */
double VariationalMC::computeDoubleDerivative(double psiLow, double psi,double psiHigh) {
    return (psiLow - 2 * psi + psiHigh) / h2;
}

double VariationalMC::computeFirstDerivative(double psiLow, double psiHigh) {
    return (psiHigh - psiLow) / (2 * h);
}

void VariationalMC::updateRmatrix(const mat &r, mat &R) {

    for (int i = 0; i < nParticles; i++) {

        // Compute magnitude of r(i,:) position vector.
        R(i,i) = 0.0;
        for (int k = 0; k < nDimensions; k ++) {
            R(i,i) += r(i,k) * r(i,k);
        }
        R(i,i) = sqrt(R(i,i));

        for (int j = (i + 1); j < nParticles; j++) {

            // Compute magnitude of r(i,:) - r(j,:) position vector.
            R(i,j) = 0.0;
            for (int k = 0; k < nDimensions; k++) {
                R(i,j) += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            R(i,j) = sqrt(R(i,j));
        }
    }


    // Pseudo-code outline of function:
    //
    //    for i=0..nParticles
    //        compute R_ii
    //
    //            for j = i+1..nParticles
    //                compute R_ij
}


/* Updates the distance matrix when we have just changed one coordinate of particle "i"
 * (like we do when computing the derivative)*/
void VariationalMC::updateForDerivative(mat &R, const mat &r, int i){
    vec dx(nDimensions);

    double dxx;

    for(int k=0; k<i; k++) {
        for(int l =0;l<nDimensions;l++){
            dxx = r(i,l) - r(k,l); // [l+i*nDimensions]-r[l+k*nDimensions];   //this may need to be changed to r[i,l] - r[k,l]
                                                          // (likewise for the next loops)
            dx(l) = dxx*dxx;
        }

        R(k,i) = sqrt(sum(dx)); //R is the matrix of distances

    }


    for(int k=i+1;k<nParticles;k++){
        for(int l =0;l<nDimensions;l++){

            dxx = r(i,l) - r(k,l);;
            dx(l) = dxx * dxx;

        }

        R(i,k) = sqrt(sum(dx)); //R is the matrix of distances

    }

    for(int l =0;l<nDimensions;l++){
        dxx = r(i,l); //r[l+i*nDimensions]*r[l+i*nDimensions];
        dx(l) = dxx*dxx;
    }

    R(i,i) = sqrt(sum(dx));
}


/* Computes the quantum force. Uses numerical differentiation to find the gradient of Psi. */
vec VariationalMC::computeQuantumForce(mat &R, mat &r, double psi) {

    vec gradient(nDimensions*nParticles);
    double psiHigh, psiLow;

    for (int k = 0; k < nParticles; k++) {
        for (int o = 0; o < nDimensions; o++) {

            // Compute psi(R - dR).
            r(k,o) -= h;
            updateForDerivative(R, r, k);
            psiLow = computePsi(R);

            // Compute psi(R + dR).
            r(k,o) += 2 * h;
            updateForDerivative(R, r, k);
            psiHigh = computePsi(R);

            // Return to original configuration, R.
            r(k,o) -= h;

            //cout << "psi= " << psi << "   psiLow=" << psiLow << "   psiHigh=" << psiHigh << endl;
            gradient(nDimensions*k+o) = computeFirstDerivative(psiLow, psiHigh);
        }
    }

    gradient *= 2 / psi;

    cout << gradient << endl;

    // greensfunction += 0.5∗(qforceold[i][j]+qforcenew[i][j]) ∗
    //                   (D∗timestep∗0.5∗(qforceold[i][j] − qforcenew[i][j])−rnew[i][j]+rold[i][j]);

    return gradient;
}


