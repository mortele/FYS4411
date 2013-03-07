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
        double  computeEnergyNumerical(mat&, mat&, double);
        double  computeDoubleDerivative(double, double, double);
        double  computeFirstDerivative (double, double);
        double  gaussian_deviate(long* idum);
        double  psi_s1(double);
        double  psi_s2(double);
        void    updateRmatrix(const mat&, mat&);
        void    updateForDerivative(mat&, const mat&, int );
        void    fillSpinMatrix(mat&);
        vec     computeQuantumForce(mat&, mat&, double);

    public:
        VariationalMC();
        double runMetropolis(double, double);
};

#endif // VARIATIONALMC_H
