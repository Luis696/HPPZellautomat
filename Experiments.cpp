#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <cmath>
#include <assert.h>

#define N 2 // 0010
#define S 8 // 1000
#define W 1 // 0001
#define E 4 // 0100



void create_matrix(int * ***Matrix, int nrows, int ncols, int nlayers);
void destroy_matrix(int ***Matrix, int nrows);
void printGrid(int **Matrix, int GridSize);
void initializeGrid(int **Matrix, int nrows, int ncols, int my_id);
void create_2D_matrix(int ***Matrix, int nrows, int ncols);

int main(int argc, char** argv){
    // ---------------  Initialize the MPI Environment:
    MPI_Status status;
    MPI_Init(&argc, &argv);
    // Get the number of processes avaiable
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    int processorGridSize = sqrt(num_procs);
    // Get the id of the process
    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);


    // Bestimme die Zeile und Spalte des aktuellen Prozessors
    int row = my_id / processorGridSize;
    int col = my_id % processorGridSize;

    // Bestimme die Position am Rand
    int edge = 0;
    if (row == 0) {
        edge |= N; // Am nördlichsten Ende
    }
    if (row == processorGridSize - 1) {
        edge |= S; // Am südlichsten Ende
    }
    if (col == 0) {
        edge |= W; // Am westlichsten Ende
    }
    if (col == processorGridSize - 1) {
        edge |= E; // Am östlichsten Ende
    }



    // ------------------ bestimme die ID der entsprechenden Prozessoren: ----------------
    // bestimme die ID der Nachbar prozessoren
    int above = -1, below = -1;
    above = (my_id - processorGridSize + num_procs) % num_procs; // Nachbar oben
    below = (my_id + processorGridSize) % num_procs; // Nachbar unten
    if (above == -1 || below == -1) {
        perror( "Error: above or below variable was not declared correctly.\n");
        exit(EXIT_FAILURE);
    }

    int left = -1, right = -1;
    if (my_id % processorGridSize == 0) {
        left = my_id + (processorGridSize - 1);
    } else {
        left = my_id - 1;
    }

    if ((my_id + 1) % processorGridSize == 0) {
        right = my_id - (processorGridSize - 1);
    } else {
        right = my_id + 1;
    }

    if (left == -1 || right == -1) {
        perror( "Error: left or right variable was not declared correctly.\n");
        exit(EXIT_FAILURE);
    }


    int Subgridsize = 10;
    int MainMatrixsize = Subgridsize * processorGridSize;

    // // Print off a hello world message
    // printf("Hello from processor, rank %d out of %d processors\n "
    //        "--------my neighbours are ---------:\n"
    //        "above: %d\n"
    //        "right: %d\n"
    //        "below: %d\n"
    //        "left:  %d\n"
    //        "edge:  %i\n\n\n", my_id, num_procs, above, right, below, left, edge);





    if (my_id != 0) {
        int **Matrix = NULL;
        create_2D_matrix(&Matrix,Subgridsize,Subgridsize);
        initializeGrid(Matrix, Subgridsize,Subgridsize, my_id);

        // send your submatrix
        MPI_Send(&(Matrix[0][0]), Subgridsize*Subgridsize, MPI_INT, 0, 0, MPI_COMM_WORLD);
        destroy_matrix(&Matrix, Subgridsize);

    }

    if (my_id == 0) {
        // Allocate and initialize the global matrix.
        int **GlobalMatrix = NULL;
        create_2D_matrix(&GlobalMatrix, MainMatrixsize, MainMatrixsize);
        initializeGrid(GlobalMatrix, MainMatrixsize,MainMatrixsize, 0);

        // Allocate and initialize the temporary (buffer) matrix.
        int **BufferMatrix = NULL;
        create_2D_matrix(&BufferMatrix, Subgridsize, Subgridsize);

        // Allocate and initialize the ownSubMatrix.
        int **ownSubMatrix = NULL;
        create_2D_matrix(&ownSubMatrix, Subgridsize, Subgridsize);
        initializeGrid(ownSubMatrix, Subgridsize,Subgridsize, 0);

        // Copy the ownSubMatrix into the GlobalMatrix.
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                GlobalMatrix[i][j] = ownSubMatrix[i][j];
            }
        }

        for (int p = 1; p < num_procs; ++p) {
            // Determine the starting coordinates for the submatrix in the global matrix.
            int start_row = (p / processorGridSize) * 10;
            int start_col = (p % processorGridSize) * 10;

            // Receive submatrices from the other processes into the buffer matrix.
            MPI_Recv(&(BufferMatrix[0][0]), 10 * 10, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Copy the values from the buffer matrix into the global matrix.
            for (int i = 0; i < 10; ++i) {
                for (int j = 0; j < 10; ++j) {
                    GlobalMatrix[start_row+i][start_col+j] = BufferMatrix[i][j];
                }
            }
        }

        printf("recveid Matrix\n");
        printGrid(GlobalMatrix, 10 * processorGridSize);

        // Clean up the dynamically allocated matrices.
        destroy_matrix(&ownSubMatrix, 10);
        destroy_matrix(&BufferMatrix, 10);
    }

    MPI_Finalize();
}


void create_matrix(int ****Matrix, int nrows, int ncols, int nlayers)
{
    assert ( *Matrix == NULL );    // Check, if *A=NULL. Empty A is necessary!

    *Matrix = (int***) malloc(nrows * sizeof(int**));
    if(*Matrix == NULL){
        perror("Failed to allocate memory for matrix");
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < nrows; i++){
        (*Matrix)[i] = (int**) malloc(ncols * sizeof(int*));
        if((*Matrix)[i] == NULL){
            perror("Failed to allocate memory for matrix");
            exit(EXIT_FAILURE);
        }

        for(int j = 0; j < ncols; j++){
            (*Matrix)[i][j] = (int*) malloc(nlayers * sizeof(int));
            if((*Matrix)[i][j] == NULL){
                perror("Failed to allocate memory for matrix");
                exit(EXIT_FAILURE);
            }
        }
    }
}
//
//
//
void destroy_2D_matrix(int ***Matrix, int nrows) {
    // Free the data block
    free(&((*Matrix)[0][0]));
    // Free the array of pointers
    free(*Matrix);
    *Matrix = NULL;  // Set the pointer to NULL to avoid dangling pointer issues
}

void print_2D_Grid(int **Matrix, int GridSize) {
    for (int i = 0; i < GridSize; ++i) {
        for (int j = 0; j < GridSize; ++j) {
            printf("%02d ", Matrix[i][j]);
        }
        printf("\n");   // Take a new line after printing each row
    }
    printf("\n");  // Take a new line after printing entire grid
}


void initialize_2D_Grid(int **Matrix, int nrows, int ncols, int my_id) {
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            Matrix[i][j] = my_id;
        }
    }
}


void create_2D_matrix(int ***Matrix, int nrows, int ncols) {
    assert(*Matrix == NULL); // Check, if *Matrix = NULL. An empty Matrix is necessary!
    *Matrix = (int **) malloc(nrows * sizeof(int *));
    if (*Matrix == NULL) {
        perror("Failed to allocate memory for matrix");
        exit(EXIT_FAILURE);
    }
    int *data = (int *) malloc(nrows * ncols * sizeof(int));
    if(data == NULL){
        perror("Failed to allocate memory for the data block of the matrix");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i< nrows; ++i){
        (*Matrix)[i] = &(data[i*ncols]);
    }
}