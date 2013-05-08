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
        double  computeKineticEnergyClosedForm(const mat&, mat &, const mat&, int);
        double  computePotentialEnergyClosedForm(const mat&);
        double  computeEnergyNumerical(mat&, mat&, double);
        double  computeDoubleDerivative(double, double, double);
        double  computeFirstDerivative (double, double);
        double  gaussian_deviate(long* idum);
        double  psi_s1(double);
        double  psi_s2(double);
        double  psi_s1_derivative(double, double coord);
        double  psi_s2_derivative(double, double coord);
        double  psi_s1_doubleDerivative(double);
        double  psi_s2_doubleDerivative(double);
        double  psi_s1_alphaDerivative(double r);
        double  psi_s2_alphaDerivative(double r);
        double  computeSlaterRatio(const mat&, double, mat&, int, int p);
        double  slaterPsi(double, mat &r, int i, int);
        double  psiDerivative(double, mat&,int,int, int);
        double  psiDoubleDerivative(double,mat &, int,  int);
        double  computeCorrelation(const mat&);
        void    updateRmatrix(const mat&, mat&);
        void    updateForDerivative(mat&, const mat&, int);
        void    fillSpinMatrix(mat&);
        void    updateSlaterInverse(mat&, const mat&, const mat&, const mat&, mat &r, int, int, double);
        void    evaluateSlater(mat&, mat&, mat &r, int);
        //vec   computeQuantumForce(mat&, mat&, double);
        void    computeSlaterGradient(mat&, mat&, mat& , mat &, double, int);
        void    computeJastrowGradient(const mat& , mat&, int);
        void    computeJastrowLaplacian(const mat& , mat&, int);
        double  computeRc( const mat& ,const mat&, int);
        void    updateCorrelationsMatrix( mat& ,const mat&,const mat&,int);
        void    fillCorrelationsMatrix( mat&,const mat&,const mat&);
        void    computeQuantumForce(mat &, mat &, mat &, mat &, mat &, double &);
        double  computeJastrowEnergy(const mat& , mat& , mat&);
        double  computeJastrowBetaDerivative(const mat&, const mat&);
        double  slaterAlphaDerivative(int, double);
        double  computeSlaterAlphaDerivative(const mat&, const mat&, const mat&);
        void    updateVariationalGradient(mat&,mat&,const mat&,const mat&,const mat&, const mat&, double, mat&, mat&);
        void    updateVariationalGradientSum(mat&,mat&,const mat&,const mat&);
        double  psi_p2z(double distance, mat &r, int i);
        double  psi_p2y(double distance, mat &r, int i);
        double  psi_p2x(double distance, mat &r, int i);
        double  psi_p2z_derivative(double distance, mat &r, int i, int j);
        double  psi_p2y_derivative(double distance, mat &r, int i, int j);
        double  psi_p2x_derivative(double distance, mat &r, int i, int j);
        double  psi_p2z_doubleDerivative(double distance, mat &r, int i);
        double  psi_p2y_doubleDerivative(double distance, mat &r, int i);
        double  psi_p2x_doubleDerivative(double distance, mat &r, int i);
    public:
        VariationalMC();
        vec runMetropolis(double, double);
};

#endif // VARIATIONALMC_H

