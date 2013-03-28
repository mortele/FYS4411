#ifndef VARIATIONALMC_H
#define VARIATIONALMC_H

#include <armadillo>
#include "wavefunctions.h"

using namespace std;
using namespace arma;


class VariationalMC {
    private:
        int     nParticles;
        int     nDimensions;
        int     nCycles;
        int     N;
        long    idum;
        double  h;
        double  h2;
        double  alph;
        double  alph2;
        double  beta;
        double  Z;
        double  stepSize;
        double  D;
        double  dt;
        vec     dx;
        mat     spins;
        WaveFunctions waveFunc;


        double  computePsi(const mat&);
        double  computePsi2(const mat&);
        double  computeEnergy(mat&, mat&, double);
        double  computeKineticEnergyClosedForm(const mat&, const mat&, const mat&, int);
        double  computePotentialEnergyClosedForm(const mat&);
        double  computeEnergyNumerical(mat&, mat&, double);
        double  computeDoubleDerivative(double, double, double);
        double  computeFirstDerivative (double, double);
        double  gaussian_deviate(long* idum);
        double  psi_s1(double);
        double  psi_s2(double);
        double  psi_s1_derivative(double, double );
        double  psi_s2_derivative(double, double );
        double  psi_s1_doubleDerivative(double);
        double  psi_s2_doubleDerivative(double);
        double  computeSlaterRatio(const mat&, double, int);
        double  slaterPsi(double, int);
        double  psiDerivative(double, double, int);
        double  psiDoubleDerivative(double, int);
        double  computeCorrelation(const mat&);
        void    updateRmatrix(const mat&, mat&);
        void    updateForDerivative(mat&, const mat&, int );
        void    fillSpinMatrix(mat&);
        void    updateSlaterInverse(mat&, const mat&, const mat&, const mat&, int, int, double);
        void    evaluateSlater(mat&, mat&, int);
        vec     computeQuantumForce(mat&, mat&, double);
        void    computeSlaterGradient(mat&, mat&, mat& , mat &, double, int);

    public:
        VariationalMC();
        double runMetropolis(double, double);
};

#endif // VARIATIONALMC_H

