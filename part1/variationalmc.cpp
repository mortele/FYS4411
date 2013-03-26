#include <armadillo>
#include <time.h>
#include <iomanip>
#include <cmath>

#include "variationalmc.h"
#include "wavefunctions.h"
#include "lib.cpp"


using namespace std;
using namespace arma;

/* VMC constructor. */
VariationalMC::VariationalMC() :
    nParticles  (4),
    nDimensions (3),
    nCycles     (500000),
    N       (2 * nCycles / 10),
    idum    (time(0)),
    h       (0.00001),
    h2      (h * h),
    alph    (1.0),
    alph2   (alph * alph),
    beta    (1.0),
    Z       (nParticles),
    stepSize(0.0001),
    D       (0.5),
    dt      (0.03), // 0.0007
    dx      (zeros(nDimensions)),
    spins   (zeros(nParticles,nParticles)) {
}


/* Runs the Metropolis algorithm nCycles times. */
double VariationalMC::runMetropolis(double alpha, double beta) {
    this->alph  = alpha;
    this->alph2 = alph * alph;
    this->beta  = beta;

    mat coordinatesNew  = zeros<mat>(nParticles, nDimensions);
    mat coordinatesOld  = zeros<mat>(nParticles, nDimensions);
    mat Rnew            = zeros<mat>(nParticles, nParticles);   // Matrix of distances and magnitudes.
    mat Rold            = zeros<mat>(nParticles, nParticles);
    mat slater          = zeros<mat>(nParticles/2, nParticles/2);
    mat slaterOldUp     = zeros<mat>(nParticles/2, nParticles/2);
    mat slaterOldDown   = zeros<mat>(nParticles/2, nParticles/2);
    mat slaterNewUp     = zeros<mat>(nParticles/2, nParticles/2);
    mat slaterNewDown   = zeros<mat>(nParticles/2, nParticles/2);

    vec quantumForceOld(nDimensions * nParticles);
    vec quantumForceNew(nDimensions * nParticles);


    double ecoeff          = 0.0;
    double Rsd             = 0.0;
    double R               = 0.0;
    double newWaveFunction = 0.0;
    double oldWaveFunction = 0.0;

    double correlationOld  = 0.0;
    double correlationNew  = 0.0;

    double energy          = 0.0;
    double energy2         = 0.0;

    double energySum       = 0.0;
    double energy2Sum      = 0.0;

    double energyUp        = 0.0;
    double energyDown      = 0.0;
    double energyPot       = 0.0;

    double greensFunction  = 0.0;

    int    accepted        = 0;
    int    t               = 0;

    double randI;
    int    iRand;

    // Fills the spin matrix with relevant values for Beryllium.
    fillSpinMatrix(this->spins);

    // Fill coordinates matrix with random values.
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            coordinatesNew(i,j) = (ran0(&idum)-0.5) / (0.5*alph);
            coordinatesOld(i,j) = coordinatesNew(i,j);
        }
    }

    // Updates distances and magnitudes in this initial state.
    updateRmatrix(coordinatesNew, Rnew);
    updateRmatrix(coordinatesOld, Rold);

    // Evaluate the slater determinand, and find its inverse.
    evaluateSlater(slaterOldUp, Rnew, 0);  // Up --> k=0.
    evaluateSlater(slaterOldDown, Rnew, 1); // Down --> k=1.

    // Invert the slater derterminands.
    slaterOldUp   = slaterOldUp.i();
    slaterOldDown = slaterOldDown.i();
    slaterNewUp   = slaterOldUp;
    slaterNewDown = slaterOldDown;



    // Compute the intial energy.
    energyUp   = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterOldUp, 0);
    energyDown = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterOldDown, 1);

    energyPot = computePotentialEnergyClosedForm(Rnew);
    //cout << energyUp + energyDown + energyPot<< endl;

    // Compute the correlation factor in the initial state.
    correlationOld  = computeCorrelation(Rnew);

    // Compute the wave function in this initial state.
    oldWaveFunction = computePsi(Rnew);

    // Compute the quantum force in the intial state.
    quantumForceOld = computeQuantumForce(Rnew, coordinatesNew, oldWaveFunction);

    // Metropolis loop.
    for (int k = 0; k < nCycles; k++) {

        // Suggest new positions for a single particle, i.e. new state.
        randI = ran0(&idum) * nParticles;
        iRand = floor(randI);

        for (int j = 0; j < nDimensions; j++) {
            // Brute force way:
            coordinatesNew(iRand,j) += (ran0(&idum)-0.5) * stepSize;

            // Importance sampled way:
            //coordinatesNew(iRand, j) += gaussian_deviate(&idum) * sqrt(dt) +
               //                       quantumForceOld(nDimensions*iRand+j) * dt; // * 2 * D;
        }


        // Compute the wavefunction in this new state.
        updateForDerivative(Rnew, coordinatesNew,iRand); // updates R.
        newWaveFunction = computePsi(Rnew);

        // Compute the quantum force in this new state.
//        quantumForceNew = computeQuantumForce(Rnew, coordinatesNew, newWaveFunction);

//        // Compute the inside of the exponential term of the difference between Greens functions.
//        greensFunction = 0.0;
//        for (int i = 0; i < nParticles; i++) {
//            for (int j = 0; j < nDimensions; j++) {
//                greensFunction += 0.5 * (quantumForceOld(nDimensions*i+j) +
//                                  quantumForceNew(nDimensions*i+j)) *
//                                  (D * dt * 0.5 * (quantumForceOld(nDimensions*i+j) -
//                                  quantumForceNew(nDimensions*i+j)) -
//                                  coordinatesNew(i,j) + coordinatesOld(i,j));
//            }
//        }

//        // Compute the fraction GreensF(new) / GreensF(old).
//        greensFunction = exp(greensFunction);

        // Check if the suggested move is accepted, brute force way.
        // ecoeff = newWaveFunction * newWaveFunction / (oldWaveFunction * oldWaveFunction);


        if (iRand >= (nParticles/2)) {
            // Down part.
            Rsd = computeSlaterRatio(slaterOldDown, Rnew(iRand,iRand), iRand-(nParticles/2));
        } else {
            // Up part.
            Rsd = computeSlaterRatio(slaterOldUp, Rnew(iRand,iRand), iRand);
        }

//        correlationNew = computeCorrelation(Rnew);
//        R = Rsd * correlationNew / correlationOld;
        R = Rsd;


        // Check if the suggested move is accepted, importance sampled way.
        // ecoeff =greensFunction * newWaveFunction * newWaveFunction / (oldWaveFunction * oldWaveFunction);

        ecoeff =  R * R;
        //cout << ecoeff << endl;
        if (ecoeff > ran0(&idum)) {
            // Accept new step, calculate new energy.
            accepted++;
            coordinatesOld.row(iRand) = coordinatesNew.row(iRand);
            quantumForceOld = quantumForceNew;
            oldWaveFunction = newWaveFunction;
            correlationOld  = correlationNew;
            Rold = Rnew; // Denne burde ikke være her, men hvis ikke blir resultatet dårlig
            //energy          = computeEnergyNumerical(Rnew, coordinatesNew, newWaveFunction);


            // Check which slater determinand we need to change -- up or down.
            if (iRand >= (nParticles/2)) {
                // Down part.
                updateSlaterInverse(slaterNewDown, slaterOldDown, Rnew, Rold, iRand -nParticles/2, 1, ecoeff);
//                evaluateSlater(slaterNewDown,Rnew,1);
//                cout << slaterNewDown << endl;
//                slaterNewDown = slaterNewDown.i();
//                cout << "inv=" <<slaterNewDown << endl;
//                cout << "R=" << Rnew << endl;
                slaterOldDown = slaterNewDown;

            } else {
                // Up part.
                updateSlaterInverse(slaterNewUp, slaterOldUp, Rnew, Rold, iRand, 0, ecoeff);
//                evaluateSlater(slaterNewUp, Rnew,0);
//                slaterNewUp = slaterNewUp.i();
                slaterOldUp = slaterNewUp;
            }


            // Compute the energy.
            energyUp   = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterNewUp, 0);
            energyDown = computeKineticEnergyClosedForm(Rnew,coordinatesNew,slaterNewDown, 1);
            energyPot  = computePotentialEnergyClosedForm(Rnew);
            energy = energyUp + energyDown + energyPot;
            //cout << ecoeff << endl;
            //energy = computeEnergyNumerical(Rnew, coordinatesNew, computePsi(Rnew));


        } else {
            // Reject suggested step, energy remains as before.
            coordinatesNew.row(iRand) = coordinatesOld.row(iRand);

            // Check which slater determinand we need to change -- up or down.
//            if (iRand >= (nParticles/2)) {
//                // Down part.
//               slaterNewDown = slaterOldDown;
//            } else {
//                // Up part.
//                slaterNewUp = slaterOldUp;
//            }
        }


        // Add energy of this state to the energy sum.
        energySum  += energy;
        energy2Sum += energy * energy;
        //cout << energyUp + energyDown << endl;

        // Throw away the first N samples.
        if (k == N) {
            energySum  = 0.0;
            energy2Sum = 0.0;
        }


    }

    //vec quantumForce(nParticles * nDimensions);
    //quantumForce = computeQuantumForce(Rnew ,coordinatesNew, newWaveFunction);

    // Calculate the expected value of the energy, the energy squared, and the variance.
    energy  = energySum  / (nCycles - N);
    energy2 = energy2Sum / (nCycles - N);

    cout << "<E>  = " << energy << endl;
    cout << "<E²> = " << energy2 << endl;
    cout << "Variance  = " << (energy2 - energy*energy)/sqrt(nCycles) << endl;
    cout << "Std. dev. = " << sqrt((energy2 - energy*energy)/sqrt(nCycles)) << endl;
    cout << "Accepted steps / total steps = " << ((double) accepted) / nCycles << endl;
    cout << "Total steps, nCycles = " << nCycles << endl;

    return energy;
}


/* Old version of compute psi. */
double VariationalMC::computePsi2(const mat &R) {
    double returnVal = exp(-alph * ( R(0,0) + R(1,1) ) + R(0,1) / (2 * (1  + beta * R(0,1))));
    return returnVal;
}


/* New version of compute psi. Slater determinands, bitches! */
double VariationalMC::computePsi(const mat &R) {
    // Compute the correlation part of the wave function.
    double correlation = 1.0; //computeCorrelation(R);

    // Compute the Slater determinand part of the wave function.
    double detup   = psi_s1(R(0,0)) * psi_s2(R(1,1)) - psi_s1(R(1,1)) * psi_s2(R(0,0));
    double detdown = psi_s1(R(2,2)) * psi_s2(R(3,3)) - psi_s1(R(3,3)) * psi_s2(R(2,2));

    return detup * detdown * correlation;
}


/* Computes the local energy of a state defined by position matrix r, and distance matrix R.
 * EL = 1/psi * H * psi */
double VariationalMC::computeEnergy(mat &R, mat &r, double psi) {

    double b1 = beta * R(0,1);
    double b2 = 1 + b1;
    double b3 = 1/(2 * b2 * b2);
    double prikk = r(0,0) * r(1,0) +  r(0,1) * r(1,1) + r(0,2) * r(1,2);

    double E_L1 = (alph - Z) * (1 / R(0,0)  + 1 / R(1,1)) + 1 / R(0,1) - alph2; // (alph - Z) +  + 1 / R(0,1)
    double E_L2 = E_L1 + b3 * ( (alph * (R(0,0) + R(1,1))) / (R(0,1))  * (1 - (prikk / (R(0,0) * R(1,1)))) - b3 - 2 / R(0,1) + ((2*beta) / b2));
    return E_L2;
}

double VariationalMC::computeKineticEnergyClosedForm(const mat& R, const mat& r, const mat& slater, int k) {
    double returnVal = 0.0;
    for (int i = 0; i < nParticles/2; i++) {
        for (int j = 0; j < nParticles/2; j++) {
            //cout << "i,j,k=" <<i << " " << j << " " << " " << k << " " <<psiDoubleDerivative(R(i+2*k,i+2*k), j) << endl;
            returnVal += psiDoubleDerivative(R(i+2*k,i+2*k), j) * slater(j,i);
        }
    }
    //cout << "slaterkin = " << returnVal /(-2.0) << endl;
    return returnVal / (-2.0);
}

double VariationalMC::computePotentialEnergyClosedForm(const mat& R) {
    double returnVal = 0.0;

    // Compute the commutative part of H.
    for(int i = 0; i < nParticles; i++) {
        returnVal -= Z/R(i,i);
    }
    return returnVal;
}


/* Computes the local energy, by numerical differentiation, of a state defined by position matrix r, and distance matrix R.
 * EL = 1/psi * H * psi */
double VariationalMC::computeEnergyNumerical(mat &R, mat &r, double psi) {
    double psil, psih;
    double E1  = 0;
    double E2  = 0;
    double psi2, psi3;

    // Compute the commutative part of H.
    for(int i = 0;i<nParticles; i++) {
        E1 -= Z/R(i,i);
        //for(int j = i+1; j<nParticles; j++) {
        //    E1 += 1/R(i,j);
        //}
    }

    for(int i = 0; i<nParticles;i++){
        for(int j = 0; j<nDimensions;j++){
            psi3 = computePsi(R);
            r(i,j) += h;
            updateForDerivative(R, r, i);
            psih = computePsi(R);

            r(i,j) -= 2 * h;
            updateForDerivative(R, r, i);
            psil = computePsi(R);

            // Reset R back to oringial positions.
            r(i,j) += h;
            updateForDerivative(R, r, i);

            psi2 = computePsi(R);
            E2 -= computeDoubleDerivative(psil, psi2, psih);
        }
    }
    if (fabs(E2 / (2 * psi) + E1)>80) {
        cout << "hei";
        //return -20;
    }
    return E2 / (2 * psi) + E1;
}


/* Computes a numerical approximation to the double derivative of psi. */
double VariationalMC::computeDoubleDerivative(double psiLow, double psi,double psiHigh) {
    return (psiLow - 2 * psi + psiHigh) / h2;
}

double VariationalMC::computeFirstDerivative(double psiLow, double psiHigh) {
    return (psiHigh - psiLow) / (2 * h);
}


/* Takes two uniformely distributed random numbers rand1, and rand2, and transforms
 * them into two gaussian distributed random numbers. */
double VariationalMC::gaussian_deviate(long * idum) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if ( idum < 0) iset =0;
    if (iset == 0) {
        do {
            v1 = 2.*ran0(idum) -1.0;
            v2 = 2.*ran0(idum) -1.0;
            rsq = v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.);
        fac = sqrt(-2.*log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    } else {
        iset =0;
        return gset;
    }
} // end function for gaussian deviates


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

    double dxx, sum;
    for(int k=0; k<i; k++) {
        sum = 0;
        for(int l =0;l<nDimensions;l++){
            dxx = r(i,l) - r(k,l); // [l+i*nDimensions]-r[l+k*nDimensions];   //this may need to be changed to r[i,l] - r[k,l]
                                                         // (likewise for the next loops)
            sum += dxx*dxx;
        }
        R(k,i) = sqrt(sum); //R is the matrix of distances
    }

    for(int k=i+1;k<nParticles;k++){
        sum = 0;
        for(int l =0;l<nDimensions;l++){
            dxx = r(i,l) - r(k,l);;
            sum += dxx*dxx;
        }
        R(i,k) = sqrt(sum); //R is the matrix of distances
    }

    sum = 0;
    for(int l =0;l<nDimensions;l++){
        dxx = r(i,l); //r[l+i*nDimensions]*r[l+i*nDimensions];
        sum += dxx*dxx;
    }
    R(i,i) = sqrt(sum);
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
            updateForDerivative(R, r, k);

            gradient(nDimensions*k+o) = computeFirstDerivative(psiLow, psiHigh);
        }
    }
    return gradient;
}


double VariationalMC::psi_s1(double distance){
    return exp(-alph*distance);
}

double VariationalMC::psi_s2(double distance){
    return (1-alph*distance/2.)*exp(-alph*distance/2);
}

double VariationalMC::psi_s1_derivative(double distance) {

}

double VariationalMC::psi_s2_derivative(double distance) {

}

double VariationalMC::psi_s1_doubleDerivative(double distance) {
    return (alph2 - 2 * alph / distance) * exp(-alph * distance);
}

double VariationalMC::psi_s2_doubleDerivative(double distance) {
    return (5.0/4.0 * alph2 - 2 * alph / distance - alph2*alph * distance / 8.) * exp(-alph * distance / 2.);
}

double VariationalMC::psiDerivative(double distance, int j) {

}

double VariationalMC::psiDoubleDerivative(double distance, int j) {
    if (j == 0) {
        return psi_s1_doubleDerivative(distance);
    } else if (j == 1) {
        return psi_s2_doubleDerivative(distance);
    } else {
        return 0;
        cout << "Error in psiDoubleDerivative" << j  << endl;
    }
}


double VariationalMC::computeSlaterRatio(const mat& slaterInverseOld, double distance, int i) {
    double Rsd = 0;
    for (int j = 0; j < (nParticles / 2); j++) {
        Rsd += slaterPsi(distance, j) * slaterInverseOld(j,i);
    }
    return Rsd;
}


/* Choses which wave function to calculate from the int j, then calls the corresponding psi_j. */
double VariationalMC::slaterPsi(double distance, int j) {
    if (j == 0) {
        return psi_s1(distance);
    } else if (j == 1) {
        return psi_s2(distance);
    } else {
        return 0;
        cout << "Error is slaterPsi" << j<< endl;
    }
}




/* Computes the correlation factor of the wave function. */
double VariationalMC::computeCorrelation(const mat& R) {
    // Compute the correlation part of the wave function.
    double correlation = 0.0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = (i+1); j < nParticles; j++) {
            correlation += spins(i,j) * R(i,j) / (1 + beta * R(i,j));
        }
    }
    correlation = exp(correlation);
    return correlation;
}


/* Updates the inverse of one part of the slater determinand (up, or down part). */
void VariationalMC::updateSlaterInverse(mat&        slaterNew,
                                        const mat&  slaterOld,
                                        const mat&  Rnew,
                                        const mat&  Rold,
                                        int         particle,
                                        int         k,
                                        double      R) {
    for (int i = 0;i<nParticles/2; i++) {
        for (int j = 0; j<nParticles/2; j++) {
            if(j != particle) {
                double sum = 0;
                for( int l = 0; l< nParticles/2; l++) {
                    sum += slaterOld(l,j) * slaterPsi(Rnew(particle+2*k,particle+2*k), l);     //slaternew(particle,l) * slaterold(l,j);
                }
                slaterNew(i,j) = slaterOld(i,j) - slaterOld(i,particle) * sum / R;
            }
            else {
                double sum = 0;
                for( int l = 0; l< nParticles/2; l++) {
                    sum += slaterPsi(Rold(particle+2*k,particle+2*k), l) * slaterOld(l, j);
                }
                slaterNew(i,j) = slaterOld(i,particle)*sum / R;
            }
        }
    }
}

// Evaluates the slater determinand. k=0 means the spin up determinand, k=1 is spin down;
void VariationalMC::evaluateSlater(mat& slater, mat& R, int k) {
    // Offset is N/2 for the spin down determinand, 0 for spin up.
    int offset = (k*nParticles/2);

    for (int i = 0; i < (nParticles/2); i++) {
        for (int j = 0; j < (nParticles/2); j++) {
            slater(i,j) = slaterPsi(R(i+offset,i+offset), j);
        }
    }
}




/* Fills the matrix of spin values, which is used to compute the correlation
part of the wavefunction. */
void VariationalMC::fillSpinMatrix(mat& spin) {
    // Set which electrons have spins up, and which have spins down.
    vec electronSpins(nParticles);
    electronSpins.zeros();
    for (int i = 0; i < (nParticles / 2); i++) {
        electronSpins(i) = 1;
    }

    // Set the electron-electron spin interaction matrix.
    for (int i = 0; i < nParticles; i++) {
        for (int j = i+1; j < nParticles; j++) {
            if (electronSpins(i) != electronSpins(j)) {
                spin(i,j) = 0.5;
            } else {
                spin(i,j) = 0.25;
            }
        }
    }
}









