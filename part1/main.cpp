#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <armadillo>

#include "lib.h"
#include "variationalmc.h"
#include "variationalloop.h"

//Old version of variationalmc.h, using numerical double derivatives.
//#include "varmc.h"

using namespace std;
using namespace arma;



int main() {
    // Initialize instances of variational monte carlo classes.
    VariationalMC   m;
    VariationalLoop l;

    // Initialize relevant values for the Helium atom..
    l.initialize_helium();

    // Run a loop over variational parameters.
    l.run();

    return 0;
}


