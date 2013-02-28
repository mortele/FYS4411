#ifndef VARMC_H
#define VARMC_H

#include <armadillo>

using namespace std;
using namespace arma;


class VarMC {
    private:
        int     nParticles;
        int     nDimensions;
        int     nCycles;
        int     N;
        long    idum;
        double  charge;
        double  h;
        double  h2;
        double  alph;
        double  alph2;
        double  beta;
        double  Z;
        double  stepSize;
        vec dx;

        double computePsi(const mat &);
        double computeEnergy(mat &, mat &, double);
        double computeDoubleDerivative(double, double, double);
        void   updateRmatrix(const mat &, mat &);
        void   updateForDerivative(mat &, const mat &, int );

    public:
        VarMC();
        double runMetropolis(double, double);
};

#endif // VARIATIONALMC_H
