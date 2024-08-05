#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <cmath>
#include <assert.h>

//#define GRID_SIZE 0 // lege die größe des grids fest
#define numb_iterations 20

// Bits für die Partikelrichtungen (siehe Zustandsübergangstabelle)
#define N 2 // 0010
#define S 8 // 1000
#define W 1 // 0001
#define E 4 // 0100
char filename[250]; // Puffer für den Dateinamen, ausreichend groß
char path[] = "Grids/";
/* mögliche Zustände der Zelle:
 * 0010  2  1001  9  1011 11
 * 0001  1  1100 12  1101 13
 * 1000  8  0110  6  1110 14
 * 0100  4  0101  5  0111  7
 * 0011  3  1010 10  1111 15
 *
 * Collosions: 5, 10
 * 5 = E & W collosion
 * 10 = N & S collosion
*/

// definiere alle Funktionen die wir brauchen
void create_matrix(int * ***Matrix, int nrows, int ncols, int nlayers);
void destroy_matrix(int * ***Matrix, int nrows, int nlayers);
void create_vector(int **vector, int SubGridSize);
void destroy_vector(int **vector);

void initializeGrid(int ***Matrix, int nrows, int , int nlayers,int initializing_value);
void moveParticles(int ***Matrix, int SubGridSize, int OLD,int NEW, int my_id, int edge);
void handleCollisions(int ***Matrix, int SubGridSize, int OLD);
void printGrid(int ***Matrix, int GridSize, int layer);
void saveGridToFile(int ***Matrix, int GridSize, int layer, const char* filename);
void share_edges(int my_id, int ***Matrix, int SubGridSize, int NEW,int OLD, int above, int below, int left, int right,
                 int *top_edge_roll_over, int *bottom_edge_roll_over, int *left_edge, int *right_edge, int edge);
void gatherSubgrids(int ***GlobalMatrix, int *** BufferMatrix, int *** ownSubMatrix, int subGridSize, int subGridLayers, int MainMatrixsize , int num_procs, int processorGridSize, const char *filename);
void flipDirections(int ***Matrix, int nrows, int ncols, int nlayers);

int main(int argc, char** argv) {
    // printf("\n\n--------- Environmet settings:---------- \n");
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

    // // -------------- Initialize the OpenMP Environment:
    // int num_of_threads = 1;
    // //num_of_threads = omp_get_max_threads();
    // printf("- OpenMP max threads %i \n",num_of_threads);
    // omp_set_num_threads(num_of_threads);
    // printf("--------------- finished -----------------\n\n");
    // // -------------------------------------------------


    // ------------------ bestimme die ID der entsprechenden Prozessoren: ----------------

    // Bestimme ob die Subgrids welche die Ränder Bilden. Weise den Subgrids entsprechend Ihrer Randposition einen Wert.
    // Ecken haben Ränder in zwei Richtungen:
    //     linke obere Ecke: N & W = 1+2 = 3
    //     rechte obere Ecke: N & E = 2+4 = 6
    //     rechte untere Ecke: O & S = 4+8 = 12
    //     linke untere Ecke: S & W = 8+1 = 9
    // Keine Randposition = 0
    //
    // Bestimme die Zeile und Spalte des aktuellen Prozessors im Gesamtgrid
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

    if(my_id == 0) {
        // Print off a hello world message from each processor
        printf("Hello from processor %s, rank %d out of %d processors\n "
               "--------my neighbours are ---------:\n"
               "above: %d\n"
               "right: %d\n"
               "below: %d\n"
               "left:  %d\n"
               "edge value: %i \n\n\n", processor_name, my_id, num_procs, above, right, below, left,edge);
        //---------------------------------------------------------------------------------
    }


    int subGridSize = 10;
    int subGridLayers = 2;
    int MainMatrixsize = subGridSize * processorGridSize;

    // allocate SubMatrices for all processors
    int*** SubMatrix = NULL;
    //create matrix buffer:
    create_matrix(&SubMatrix,subGridSize,subGridSize,subGridLayers);
    initializeGrid(SubMatrix,subGridSize,subGridSize,subGridLayers, 0);
    printf("created SubGrid processor %i ... \n",my_id);

    // Allocate the edge buffers for all subgrids of all processors
    int *top_edge_roll_over = NULL, *bottom_edge_roll_over= NULL;
    int *left_edge= NULL, *right_edge= NULL;
    create_vector(&top_edge_roll_over, subGridSize);
    create_vector(&bottom_edge_roll_over, subGridSize);
    create_vector(&left_edge, subGridSize);
    create_vector(&right_edge, subGridSize);

    // // print the initial grid of every processor to console & saving it
    // printf("Initial grid layer 0:\n");
    // printGrid(SubMatrix, subGridSize, 0);
    // printf("Initial grid layer 1:\n");
    // printGrid(SubMatrix, subGridSize, 1);
    // sprintf(filename, "../Grids/grid_start.txt");
    // sending all subgrids to process 0 so it can fill the main grid with it:
    //
    if (my_id != 0) {
        // send your submatrix
        MPI_Send(&(SubMatrix[0][0][0]), subGridSize * subGridSize * subGridLayers, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    // // MPI_Barrier(MPI_COMM_WORLD);
    //
    int ***GlobalMatrix = NULL;  //
    int ***BufferMatrix = NULL;
    if(my_id == 0) {

        // allocate the global matrix for processor 0
        create_matrix(&GlobalMatrix, MainMatrixsize, MainMatrixsize, subGridLayers);
        // initialization of GlobalMatrix
        initializeGrid(GlobalMatrix, MainMatrixsize, MainMatrixsize, subGridLayers, 0);

        // case collidinig particles:
        SubMatrix[2][2][0] = N; // setting an initial particle
        SubMatrix[4][2][0] = S; // setting an initial particle
        SubMatrix[5][4][0] = W; // setting an initial particle
        SubMatrix[5][6][0] = E; // setting an initial particle

        // set initial particles
        int particletypes[4] = {N, S,  W, E};
        for(int p = 0; p < 20; p++) {
            int x = rand() % subGridSize; // Get a random x-coordinate.
            int y = rand() % subGridSize; // Get a random y-coordinate.
            int particle = rand() % 4; // Select a random particle.
            SubMatrix[x][y][0] = particletypes[particle]; // Use the correct character for the selected particle.
        }


        // print information about Grids sizes
        printf("Number of Processors:%i \n", num_procs);
        printf("Processor Gridzize: %i x %i \n", processorGridSize, processorGridSize);
        printf("Subgridsize: %i x %i \n", subGridSize, subGridSize);
        printf("Main Matrix size: %i x %i \n\n", MainMatrixsize, MainMatrixsize);
        // Allocate and initialize the temporary (buffer) matrix.

        create_matrix(&BufferMatrix, subGridSize, subGridSize, subGridLayers);

        int step = 0;
        sprintf(filename, "%sgrid_%i.txt",path, step);
        // collecting all Subgrids from all processors
        gatherSubgrids(GlobalMatrix, BufferMatrix, SubMatrix, subGridSize,subGridLayers, MainMatrixsize, num_procs, processorGridSize, filename);
        printf("%s: \n", filename);
        printGrid(GlobalMatrix,MainMatrixsize,0);
        saveGridToFile(GlobalMatrix, MainMatrixsize, 0, filename);

    }

    //
    // // gehe sicher, dass alle Prozessoren das gemacht haben bevor Sie weiter machen
    MPI_Barrier(MPI_COMM_WORLD);

    int flipping = 1;
    int val1, val2 = 0;
    // ------------------------------------------------ start the itteration of the grid: -----------------------------------------------

    for (int step = 0; step < numb_iterations; ++step) {
        // Erzeuge den Dateinamen mit der Iterationsnummer
        int val1 = step % 2; // alternate depending on the step
        int val2 = 1 - val1; // inverse of val1
        /*  step = 0, val1 = 0, val2 = 1
            step = 1, val1 = 1, val2 = 0
            step = 2, val1 = 0, val2 = 1
            step = 3, val1 = 1, val2 = 0
            step = 4, val1 = 0, val2 = 1
         */
        // printf("value1: %i,  value2: %i \n",val1,val2);
        // CAVE: OLD grid needs to start @ layer 0 -> see initialize Grids
        handleCollisions(SubMatrix, subGridSize, val1); // check for collions FIRST and change directions of particles if needed
        // Teilen der Randwerte mit den benachbarten Prozessoren
        share_edges(my_id, SubMatrix, subGridSize, val2,val1, above, below, left, right, top_edge_roll_over, bottom_edge_roll_over,
                    left_edge, right_edge, edge);

        moveParticles(SubMatrix, subGridSize, val1, val2, my_id, edge);  // bewege alle Partikel entlang der Richtung





        // sending all subgrids to process 0 so it can fill the main grid with it:
        if (my_id == 0) {
            gatherSubgrids(GlobalMatrix, BufferMatrix, SubMatrix, subGridSize,subGridLayers, MainMatrixsize, num_procs, processorGridSize, filename);
            sprintf(filename, "%sgrid_%i.txt",path, step);
            printf("%s: \n", filename);
            printGrid(GlobalMatrix,MainMatrixsize,val2);
            saveGridToFile(GlobalMatrix, MainMatrixsize, val2, filename);
        }
        if (my_id != 0) {
            // send your submatrix
            MPI_Send(&(SubMatrix[0][0][0]), subGridSize * subGridSize * subGridLayers, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        if(step >= numb_iterations/2 && flipping == 1) {

            flipDirections(SubMatrix, subGridSize, subGridSize, 2);
            flipping = 0;
        }


    }
    // ---------------------------------------------------- finished the itteration of the grid: ----------------------------------

    // sprintf(filename, "../Grids/grid_finished.txt");
    // // sending all subgrids to process 0 after finishing itterations so it can fill the main grid with it:
    // if (my_id == 0) {
    //     // Define the filename for saving the grid
    //     gatherSubgrids(SubMatrix, subGridSize, GlobalMatrix, num_procs, processorGridSize,filename); // gather all subgrids and write it to a file
    // }
    // else {
    //     // send your submatrix
    //     MPI_Send(&(SubMatrix[0][0][val2]), subGridSize * subGridSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
    // }
    //
    //
    MPI_Barrier(MPI_COMM_WORLD);
    // Free allocated memory:
    destroy_matrix(&SubMatrix,subGridSize, 2);
    // printf("SubGrid of processor %i succesfully destroyed... \n",my_id);
    if (my_id == 0){
        destroy_matrix(&GlobalMatrix, subGridSize*processorGridSize, 2);
        // printf("Global Matrix of processor %i succesfully destroyed... \n",my_id);
    }
    destroy_vector(&top_edge_roll_over);
    destroy_vector(&bottom_edge_roll_over);
    destroy_vector(&left_edge);
    destroy_vector(&right_edge);
    // printf("vectors of processor %i succesfully destroyed... \n",my_id);

    return 0;
    MPI_Finalize();
}




void initializeGrid(int ***Matrix, int nrows, int ncols , int nlayers, int initializing_value) {
    int rand_num = 0;
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            for( int k = 0; k < nlayers; k++) {
                rand_num = rand() % 16; // Zufälliger Zustand (0 bis 15) für alle möglichen Binärzustände, siehe Zustandsübergangstabelle
                Matrix[i][j][k] = initializing_value; //rand_num;
                Matrix[i][j][k] = initializing_value; //rand_num;
            }
        }
    }


}

void moveParticles(int ***Matrix, int SubGridSize, int OLD,int NEW, int my_id, int edge) {
    // newGrid, and oldGrid is alternating -> layer 0 and 1
    // gehe durch das gesamte alte grid und wende die Zustandsübergangstabelle an
    int i, j;
    for (i = 0; i < SubGridSize; ++i) {
        for (j = 0; j < SubGridSize; ++j) {
            if (Matrix[i][j][OLD] & N) {
                if (i > 0) Matrix[i-1][j][NEW] |= N; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (i == 0 && (edge & N)) Matrix[i+1][j][NEW] |= S; // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~N; // wenn Partikel bewegt, lösche alte Position
            }

            if (Matrix[i][j][OLD] & S) {
                if (i < SubGridSize - 1) Matrix[i+1][j][NEW] |= S; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (i == SubGridSize - 1 && (edge & S)) Matrix[i-1][j][NEW] |= N; // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~S;  // wenn Partikel bewegt, lösche alte Position
            }

            if (Matrix[i][j][OLD] & W) {
                if (j > 0) Matrix[i][j-1][NEW] |= W; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (j == 0 && (edge & W)) Matrix[i][j+1][NEW] |= E;  // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~W;  // wenn Partikel bewegt, lösche alte Position
            }

            if (Matrix[i][j][OLD] & E) {
                if (j < SubGridSize - 1) Matrix[i][j+1][NEW] |= E; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (j == SubGridSize - 1 && (edge & E)) Matrix[i][j-1][NEW] |= W; // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~E;  // wenn Partikel bewegt, lösche alte Position
            }
        }
    }

}

void handleCollisions(int ***Matrix, int SubGridSize, int OLD) {
    /*handels Collisions in the new grid*/
    int i, j;
    for (i = 0; i < SubGridSize; ++i) {
        for (j = 0; j < SubGridSize; ++j) {
            int cell = Matrix[i][j][OLD];
            if ((cell & N) && (cell & S)) { // wenn in der Zelle Nord- und Südpartikel vorhanden sind
                cell &= ~(N | S); // entferne Sie die Nord- und Südpartikel
                cell |= (W | E);  // fügen West- und Ostpartikel hinzu
                //printf("!N & S Collosion!");

            }

            else if((cell & W) && (cell & E)) { // wenn in der Zelle West und Ostpartikel vorhanden sind
                cell &= ~(W | E); // entfernen Sie die West- und Ostpartikel
                cell |= (N | S);  // fügen Nord- und Südpartikel hinzu
                //printf("!W & E Collosion!");
            }
            Matrix[i][j][OLD] = cell; // setze die aktuelle Zelle im grid auf den neuen Wert
        }
    }
}



void printGrid(int ***Matrix, int GridSize, int layer) {
    for (int i = 0; i < GridSize; ++i) {
        for (int j = 0; j < GridSize; ++j) {
            printf("%02d ", Matrix[i][j][layer]);
        }
        printf("\n");  // nach jeder ferigen Zeile beginn einen neuen Absatz
    }
    printf("\n");  // nachdem das gesamte grid ausgegeben wurde mache einen Absatz
}

void saveGridToFile(int ***Matrix, int GridSize, int layer, const char* filename) {
    printf("starting to save grid to file...\n");
    FILE* file = fopen(filename, "w");  // erzeuge eine Datei und öffne sie
    if (!file) { // wenn öffnen nicht funktioniert, dann gebe mir einen Fehler aus
        perror("Failed to open file");
        return;
    }
    // sonst speichere das grid unter dem Dateinamen ab
    for (int i = 0; i < GridSize; ++i) {
        for (int j = 0; j < GridSize; ++j) {
            fprintf(file, "%02d ", Matrix[i][j][layer]);  // schreibe grid wert in die Datei
        }
        fprintf(file, "\n");  // nachdem eine ganze Zeile geschrieben wurde beende die Line
    }
    fclose(file);  // schließe die Datei wieder
    printf("saved grid to file\n");
}

// not continuouse memory allocation
// void create_matrix(int ****Matrix, int nrows, int ncols, int nlayers)
// {
//     assert ( *Matrix == NULL );    // Check, if *A=NULL. Empty A is necessary!
//
//     *Matrix = (int***) malloc(nrows * sizeof(int**));
//     if(*Matrix == NULL){
//         perror("Failed to allocate memory for matrix");
//         exit(EXIT_FAILURE);
//     }
//
//     for(int i = 0; i < nrows; i++){
//         (*Matrix)[i] = (int**) malloc(ncols * sizeof(int*));
//         if((*Matrix)[i] == NULL){
//             perror("Failed to allocate memory for matrix");
//             exit(EXIT_FAILURE);
//         }
//
//         for(int j = 0; j < ncols; j++){
//             (*Matrix)[i][j] = (int*) malloc(nlayers * sizeof(int));
//             if((*Matrix)[i][j] == NULL){
//                 perror("Failed to allocate memory for matrix");
//                 exit(EXIT_FAILURE);
//             }
//         }
//     }
// }
// void destroy_matrix(int ****Matrix, int nrows, int nlayers) {
//     // Check, if *A != NULL. non-empty A is necessary!
//     assert(*Matrix != NULL);   // If condition is FALSE, program aborts!
//
//     // free the memory
//     int i, j;
//     for (i = 0; i < nlayers; ++i) {
//         for (j = 0; j < nrows; ++j) {
//             free((*Matrix)[i][j]);
//         }
//         free((*Matrix)[i]);
//     }
//     free(*Matrix);
//
//     *Matrix = NULL;  // Set the pointer to NULL to avoid dangling pointer issues
// }


// continouse memory allocation:
void create_matrix(int ****Matrix, int nrows, int ncols, int nlayers)
{
    assert ( *Matrix == NULL );    // Check, if *A=NULL. Empty A is necessary!

    // Create a single, contiguous block of memory for your data
    int *data = (int *) malloc(nrows * ncols * nlayers * sizeof(*data));
    assert(data != NULL);

    // Allocate your pointers into this block of memory
    *Matrix = (int ***) malloc(nrows * sizeof(**Matrix));
    assert(*Matrix != NULL);

    for (int i = 0; i < nrows; ++i) {
        (*Matrix)[i] = (int **) malloc(ncols * sizeof(***Matrix));
        assert((*Matrix)[i] != NULL);

        for (int j = 0; j < ncols; ++j) {
            (*Matrix)[i][j] = &(data[(i * ncols + j) * nlayers]);
        }
    }
}
void destroy_matrix(int ****Matrix, int nrows, int ncols)
{
    if(*Matrix != NULL){
        // Free the contiguous block of memory for data.
        // It's enough to free the first element of the first column of the first row.
        free((*Matrix)[0][0]);

        // Free the memory for the array of pointers.
        for (int i = 0; i < nrows; ++i) {
            free((*Matrix)[i]);
        }

        // Free the memory for the pointer to the pointers.
        free(*Matrix);

        // Set the pointer to null preventing accidental usage after free
        *Matrix = NULL;
    }
}


//
void create_vector(int **vector, int SubGridSize) {
    // Check, if *vector == NULL. Non-empty vector is not allowed!
    assert(*vector == NULL); // If condition is FALSE, program aborts!

    // Allocate the necessary memory
    *vector = (int*) malloc(SubGridSize * sizeof(int));

    // Check if the memory allocation was successful
    if (*vector == NULL) {
        perror("Failed to allocate memory for vector");
        exit(EXIT_FAILURE);
    }
}

void destroy_vector(int **vector) {
    assert(*vector != NULL);// If condition is FALSE, program aborts!

    // Free the memory
    free(*vector);

    // Set the pointer to NULL to avoid dangling pointer issues
    *vector = NULL;
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

void destroy_2D_matrix(int ***Matrix, int nrows) {
    // Free the data block
    free(&((*Matrix)[0][0]));
    // Free the array of pointers
    free(*Matrix);
    *Matrix = NULL;  // Set the pointer to NULL to avoid dangling pointer issues
}

void share_edges(int my_id, int ***Matrix, int SubGridSize, int NEW, int OLD, int above, int below, int left, int right,
                 int *top_edge_roll_over, int *bottom_edge_roll_over, int *left_edge, int *right_edge, int edge) {

    // -------------------------------- Ränder übergeben -----------------------------
    // entnehme oberen und unteren Rand aus der Matrix
    for (int i = 0; i < SubGridSize; ++i) {
        if (Matrix[0][i][OLD] != N || (edge & N)) {
            top_edge_roll_over[i] = 0; // set array element to 0 if it is not equal to N
        } else {
            top_edge_roll_over[i] = Matrix[0][i][OLD]; // erste Zeile alle Spalten von NEW
        }

        if (Matrix[SubGridSize - 1][i][OLD] != S || (edge & S)) {
            bottom_edge_roll_over[i] = 0; // set array element to 0 if it is not equal to N
        } else {
            bottom_edge_roll_over[i] = Matrix[SubGridSize - 1][i][OLD]; // letzte Zeile alle Spalten von NEW
        }
    }

    // Alle Subgrids reichen ihren oberen rand weiter and das Subgrid über ihnen:
    // source my_id: Subgrid das seinen oberen Rand weitereicht
    // destination: Subgrid das den oberen Rand von dem Subgrid unter Ihm erhält
    MPI_Sendrecv_replace(top_edge_roll_over, SubGridSize, MPI_INT, above, 0, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Der Processor hat den oberen Rand von dem Subgrid unter Ihm empfangen und speichert ihn in die unterste Zeile
    // damit wird top_edge_roll_over ersetze durch den oberen Rand des subgrids unter uns -> unser unterer Rand
    // Kopiere empfangene Daten in die Matrix
    for (int i = 0; i < SubGridSize; ++i) {
        Matrix[SubGridSize - 1][i][NEW] |= top_edge_roll_over[i]; // unterer Rand
    }

    // Alle Subgrids reichen ihren unteren rand weiter and das Subgrid unter ihnen:
    // source my_id: Subgrid das seinen unteren Rand weiterreicht
    // destination below: Subgrid das den unteren Rand von dem Subgrid über Ihm erhält
    MPI_Sendrecv_replace(bottom_edge_roll_over, SubGridSize, MPI_INT, below, 0, above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Der Processor hat den unteren Rand von dem Subgrid über Ihm empfangen und speichert ihn in die oberste Zeile
    // damit wird bottom_edge_roll_over ersetze durch den unteren Rand des subgrids über uns -> unser neuer oberer Rand
    // Kopiere empfangene Daten in die Matrix
    for (int i = 0; i < SubGridSize; ++i) {
        Matrix[0][i][NEW] |= bottom_edge_roll_over[i]; // oberer Rand
    }


    for (int i = 0; i < SubGridSize; ++i) {
        if (Matrix[i][0][OLD] != W || (edge & W)) { // checke if partikel auch wirklich nach Westen unterwegs und ob es sich nicht um einen Rand handelt
            left_edge[i] = 0; // set array element to 0 if it is not equal to N
        } else {
            left_edge[i] = Matrix[i][0][OLD];
        }

        if (Matrix[i][SubGridSize - 1][OLD] != E || (edge & E)) {
            right_edge[i] = 0; // set array element to 0 if it is not equal to N
        } else {
            right_edge[i] = Matrix[i][SubGridSize - 1][OLD];
        }
    }

    // Sende linken Rand und empfange linken Rand von rechtem Prozessor -> damit unser neuer rechter Rand
    // Alle Subgrids reichen ihren linken Rand weiter and das Subgrid rechts neben Ihnen:
    // source my_id: Subgrid das seinen linken Rand weiterreicht
    // destination below: Subgrid das den linken Rand von dem Subgrid rechts neben Ihm erhält
    MPI_Sendrecv_replace(left_edge, SubGridSize, MPI_INT, left, 0, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // damit wird left_edge ersetze durch den linken rand des subgrids rechts neben uns -> unser neuer rechter Rand
    // Kopiere empfangene Daten in die Matrix
    for (int i = 0; i < SubGridSize; ++i) {
        Matrix[i][SubGridSize - 1][NEW] |= left_edge[i]; // rechter Rand
    }

    // Sende rechten Rand und empfange rechten Rand von linkem Prozessor -> damit unser neuer linker Rand
    // Alle Subgrids reichen ihren rechten Rand weiter and das Subgrid rechts neben Ihnen:
    // source my_id: Subgrid das seinen rechten Rand weiterreicht
    // destination below: Subgrid das den rechten Rand von dem Subgrid links neben ihm erhält
    MPI_Sendrecv_replace(right_edge, SubGridSize, MPI_INT, right, 0, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // damit wird right_edge ersetze durch den rechten rand des subgrids links neben uns -> unser neuer linker Rand
    // Kopiere empfangene Daten in die Matrix
    for (int i = 0; i < SubGridSize; ++i) {
        Matrix[i][0][NEW] |= right_edge[i]; // linker Rand
    }
}


void gatherSubgrids(int ***GlobalMatrix, int *** BufferMatrix, int *** ownSubMatrix,  int subGridSize, int subGridLayers, int MainMatrixsize , int num_procs, int processorGridSize, const char *filename) {
    // save root processors own SubMatrix in Global Matrix
    for (int i = 0; i < subGridSize; ++i) {
        for (int j = 0; j < subGridSize; ++j) {
            for(int k = 0; k<subGridLayers; k++){
                GlobalMatrix[i][j][k] = ownSubMatrix[i][j][k];
            }
        }
    }

    // save SubMatrices from all other processors
    for (int p = 1; p < num_procs; ++p) {
        // Determine the starting coordinates for the submatrix in the global matrix.
        int start_row = (p / processorGridSize) * subGridSize;
        int start_col = (p % processorGridSize) * subGridSize;

        // Receive submatrices from the other processes into the buffer matrix.
        MPI_Recv(&(BufferMatrix[0][0][0]), subGridSize * subGridSize * subGridLayers, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Copy the values from the buffer matrix into the global matrix.
        for (int i = 0; i < subGridSize; ++i) {
            for (int j = 0; j < subGridSize; ++j) {
                for(int k = 0; k<subGridLayers; k++){
                    GlobalMatrix[start_row+i][start_col+j][k] = BufferMatrix[i][j][k];
                }
            }
        }
    }
}


double** load_matrix_from_file(const char* filename, int rowCount, int colCount) {
    // Open the file
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Could not open file %s\n", filename);
        return NULL;
    }

    // Allocate matrix
    double** matrix = (double**)malloc(rowCount * sizeof(double*));
    for (int i = 0; i < rowCount; ++i)
        matrix[i] = (double*)malloc(colCount * sizeof(double));

    // Read matrix from file
    for (int i = 0; i < rowCount; ++i)
        for (int j = 0; j < colCount; ++j)
            fscanf(file, "%lf", &matrix[i][j]);

    fclose(file);

    return matrix;
}


void flipDirections(int ***Matrix, int nrows, int ncols, int nlayers) {
    // Define a lookup table to map each direction to its opposite
    int flip_map[11] = {0, E, S, 0, W, 10, 0, 0, N, 0, 5};

    for (int layer = 0; layer < nlayers; ++layer) {
        for (int row = 0; row < nrows; ++row) {
            for (int col = 0; col < ncols; ++col) {
                Matrix[row][col][layer] = flip_map[Matrix[row][col][layer]];
            }
        }
    }
}