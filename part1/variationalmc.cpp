#include <armadillo>
#include <time.h>
#include <iomanip>

#include "variationalmc.h"
#include "lib.cpp"


using namespace std;
using namespace arma;

/* VMC constructor. */
VariationalMC::VariationalMC() :
    nParticles  (2),
    nDimensions (3),
    nCycles     (1000000),
    N       (2 * nCycles / 10),
    idum    (time(0)),
    h       (0.0001),
    h2      (h * h),
    alph    (1.0),
    alph2   (alph * alph),
    beta    (1.0),
    Z       (2.0),
    stepSize(0.07),
    D       (0.5),
    dt      (0.007),
    dx      (zeros(nDimensions)) {
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
    vec quantumForceOld(nDimensions * nParticles);
    vec quantumForceNew(nDimensions * nParticles);

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
            // cout << coordinatesNew(i,j)  << endl;
        }
    }

    // Updates distances and magnitudes in this initial state.
    updateRmatrix(coordinatesNew, Rnew);
    updateRmatrix(coordinatesOld, Rold);

    // Compute the wave function in this initial state.
    oldWaveFunction = computePsi(Rnew);

    quantumForceOld = computeQuantumForce(Rnew, coordinatesNew, oldWaveFunction);
    // Metropolis loop.
    for (int k = 0; k < nCycles; k++) {
//        if (k % (nCycles / 10) == 0) {
//            cout << percent << " %" << endl;
//            percent += 10;
//        }

        // Suggest new positions for all particles, i.e. new state.

        randI = ran0(&idum) * nParticles;
        iRand = floor(randI);

        for (int j = 0; j < nDimensions; j++) {
            // Brute force way:
            coordinatesNew(iRand,j) += (ran0(&idum)-0.5) * stepSize;

            // Importance sampled way:
            //coordinatesNew(iRand, j) += gaussian_deviate(&idum) * sqrt(dt) +
              //                          quantumForceOld(nDimensions*iRand+j) * dt; // * 2 * D;
            //cout << "qforce=" << quantumForceOld(nDimensions*iRand+j) << endl;
            //cout << "normal="<<gaussian_deviate(&idum) * sqrt(dt) << endl;
            //cout << "quantum="<<quantumForceOld(nDimensions*iRand+j) * dt << endl;
            //cout << "r="<<coordinatesOld(iRand, j) << endl << endl;

        }


        // Compute the wavefunction in this new state.
        updateRmatrix(coordinatesNew, Rnew);
        newWaveFunction = computePsi(Rnew);

        // Compute the quantum force in this new state.
        quantumForceNew = computeQuantumForce(Rnew, coordinatesNew, newWaveFunction);

        greensFunction = 0.0;
        // Compute the inside of the exponential term of the difference between Greens functions.
        for (int i = 0; i < nParticles; i++) {
            for (int j = 0; j < nDimensions; j++) {
                double pX = (coordinatesNew(i,j) - coordinatesOld(i,j) - D * dt * quantumForceOld(nDimensions*i+j));
                double pY = (coordinatesOld(i,j) - coordinatesNew(i,j) - D * dt * quantumForceNew(nDimensions*i+j));
                greensFunction += -(pX*pX/(4*D*dt)) + (pY*pY / (4*D*dt));
                greensFunction += 0.5 * (quantumForceOld(nDimensions*i+j) + quantumForceNew(nDimensions*i+j)) *
                                  (D * dt * 0.5 * (quantumForceOld(nDimensions*i+j) - quantumForceNew(nDimensions*i+j)) -
                                  coordinatesNew(i,j) + coordinatesOld(i,j));
            }
        }

        // Compute the fraction GreensF(new) / GreensF(old).
        greensFunction = exp(greensFunction);

        // cout << greensFunction << endl;

        // Check if the suggested move is accepted, brute force way.
        ecoeff = newWaveFunction * newWaveFunction / (oldWaveFunction * oldWaveFunction);

        // Check if the suggested move is accepted, importance sampled way.
        //ecoeff = greensFunction * newWaveFunction * newWaveFunction / (oldWaveFunction * oldWaveFunction);

        if (ecoeff > ran0(&idum)) {
            accepted++;
            coordinatesOld  = coordinatesNew;
            quantumForceOld = quantumForceNew;
            oldWaveFunction = newWaveFunction;
            // Energy changes from previous state.
            energy = computeEnergy(Rnew, coordinatesNew, newWaveFunction);
//            cout << "r1 coordinatesx = " << coordinatesNew(0,0) << endl;
//            cout << "r1 coordinatesy = " << coordinatesNew(0,1) << endl;
//            cout << "r1 coordinatesz = " << coordinatesNew(0,2) << endl;
        } else {
            coordinatesNew = coordinatesOld;
            newWaveFunction = oldWaveFunction;
            // Energy remains unchanged.
        }

        // Add energy of this state to the energy sum.

        if (k == N) {
            energySum  = 0.0;
            energy2Sum = 0.0;
        }

        energySum  += energy;
        energy2Sum += energy * energy;
        //if (energy*energy > 1000) cout << "energy=" << energy << endl;
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

    return energy;
}


/* Computes the wavefunction in a state defined by position matrix r. */
double VariationalMC::computePsi(const mat &R) {

    double returnVal = 0.0;
    returnVal = exp(-alph * (    R(0,0) + R(1,1)   ) + R(0,1) / (2 * (1  + beta * R(0,1)))   ); // ) * exp(
    /*cout << "0,1= " << R(0,1) << endl;
    cout << "r1= " << R(0,0) << endl;
    cout << "r2= " << R(1,1) << endl;*/

    //    for (int i = 0; i < nParticles; i++) {
    //        for (int j = (i + 1); j < nParticles; j++) {
    //            returnVal += (1 / R(i,j)) * exp(-alph * (R(i,i) + R(j,j))) * (exp(R(i,j)) / (1 + beta * R(i,j)));
    //        }
    //    }
    return returnVal;
}

double VariationalMC::computePsi2(const mat &R) {
    detup = psi_s1(R(0,0))*psi_s2(R(1,1))- psi_s1(R(1,1))*psi_s2(R(0,0));
    detdown = psi_s1(R(2,2))*psi_s2(R(3,3))- psi_s1(R(3,3))*psi_s2(R(2,2));
    return detup*detdown;
}


/* Computes the local energy of a state defined by position matrix r, and distance matrix R.
 * EL = 1/psi * H * psi */
double VariationalMC::computeEnergy(mat &R, mat &r, double psi)
{

    double b1 = beta * R(0,1);
    double b2 = 1 + b1;
    double b3 = 1/(2 * b2 * b2);
    double prikk = r(0,0) * r(1,0) +  r(0,1) * r(1,1) + r(0,2) * r(1,2);

    double E_L1 = (alph - Z) * (1 / R(0,0)  + 1 / R(1,1)) - alph2 + 1 / R(0,1) - alph2;
    double E_L2 = E_L1 + b3 * ( (alph * (R(0,0) + R(1,1))) / (R(0,1))  * (1 - (prikk / (R(0,0) * R(1,1)))) - b3 - 2 / R(0,1) + ((2*beta) / b2));
    return E_L1;


    //double r1  = R(0,0);
    //double r2  = R(1,1); HALLO
//    double psil, psih;
//    double E2  = 0;
//    double E1  = 0;   // this is the commutative part of the
//                                             // hamiltonian (we can just multiply it with psi.)
//    for(int i = 0;i<nParticles; i++) {
//        E1 -= Z/R(i,i);
//        for(int j = i+1; j<nParticles; j++) {
//            E1 += 1/R(i,j);
//        }
//    }

//    for(int i = 0;i<nParticles; i++) {

//        for(int j = 0; j<nDimensions;j++) {

//            r(i,j) -= h; //r is the array of coordinates

//            updateForDerivative(R,r, i);
//            //cout << "r12=" << R(0,1)<<" for r-h" << endl;

//            psil = computePsi(R);

//            r(i,j)+=2*h;
//            updateForDerivative(R, r, i);
//            //cout << "r12=" << R(0,1)<<" for r+h" << endl;

//            psih = computePsi(R);
//            r(i,j)-=h;

//            //cout << "psil, psi, psih = " << psil << " " << psi << " " << psih << endl;
//            //cout << "ddx=" << (psil - 2 * psi + psih) / h2 << endl;
//            E2-=computeDoubleDerivative(psil, psi,psih);
//            //cout << "E2=" << E2 << endl;
//            updateForDerivative(R, r, i);   // set all values back to normal
//        }
//    }
//    //cout << "r12=" << R(0,1) << endl;
//    //cout << E2 / psi << " " << E1 << endl;
//    return E1 + E2 / (2 * psi);








//    double psil, psih;
//    double r12 = R(0,1);

//    double r1 = R(0,0);

//    double r2 = R(1,1);

//    double E1 = (-Z*(1/r1 +1/r2) +1/r12);   //this is the commutative part of the hamiltonian (we can just multiply it with psi.)

//    double E2 = 0;
//    for(int i = 0; i<nParticles;i++){
//        for(int j = 0; j<nDimensions;j++){

//            r(i,j)+=h; //r is the array of coordinates

//            psih = computePsi(R);

//            r(i,j)-=2*h;

//            psil = computePsi(R);

//            r(i,j)+=h;



//            E2-=computeDoubleDerivative(psil, psi,psih);

//            updateForDerivative(R, r, i);   // set all values back to normal

//        }
//    }


//    return E2/psi +E1;


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

            //cout << "psi= " << psi << "   psiLow=" << psiLow << "   psiHigh=" << psiHigh << endl;
            gradient(nDimensions*k+o) = computeFirstDerivative(psiLow, psiHigh);
        }
    }

    gradient *= 2 / psi;

    // cout << gradient << endl;

    // greensfunction += 0.5∗(qforceold[i][j]+qforcenew[i][j]) ∗
    //                   (D∗timestep∗0.5∗(qforceold[i][j] − qforcenew[i][j])−rnew[i][j]+rold[i][j]);

    return gradient;
}


