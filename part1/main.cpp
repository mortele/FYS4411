#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <armadillo>
#include <mpi.h>

#include "lib.h"
#include "variationalmc.h"
#include "variationalloop.h"



using namespace std;
using namespace arma;


int main(int argc, char* argv[]) {
    // Initialize MPI.
    int my_rank, num_procs;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    VariationalLoop l;

    // Initialize processes
    l.initialize_processes(my_rank);

    // Run a loop over variational parameters.
    l.run();

    MPI_Finalize();
    return 0;
}


