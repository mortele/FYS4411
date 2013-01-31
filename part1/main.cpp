#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <armadillo>

#include "lib.h"
#include "variationalmc.h"

using namespace std;
using namespace arma;



int main() {

    VariationalMC m;
    m.runMetropolis();

    return 0;

}


