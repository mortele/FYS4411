#include "wavefunctions.h"

WaveFunctions::WaveFunctions()
{

}

double WaveFunctions::psi_s1(double distance){
    return exp(-alph*distance);
}

double WaveFunctions::psi_s2(double distance){
    return (1-alph*distance)*exp(-alph*distance/2);
}



