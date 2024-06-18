#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Define my value
    int* my_values = malloc(world_size*sizeof(int));
    int* buffer_recv = malloc(world_size*sizeof(int));
    for(int i = 0; i < world_size; i++)
    {
        my_values[i] = my_rank;
        buffer_recv[i] = my_rank;
    }
 
    MPI_Alltoall(my_values, 1, MPI_INT, buffer_recv, 1, MPI_INT, MPI_COMM_WORLD);

    if(my_rank == 0){
       printf("Rank 0:");
       for(int i = 0; i < world_size; i++){
         printf("%d,", buffer_recv[i]);
       }
       printf("\n");
    }

    MPI_Finalize();
    return 0;
}
