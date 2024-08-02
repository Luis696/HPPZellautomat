#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <cmath>
#include <assert.h>

#define GRID_SIZE 5 // lege die größe des grids fest
#define numb_iterations 10

// Bits für die Partikelrichtungen (siehe Zustandsübergangstabelle)
#define N 2 // 0010
#define S 8 // 1000
#define W 1 // 0001
#define E 4 // 0100

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

void initializeGrid(int ***Matrix, int SubGridSize);
void moveParticles(int ***Matrix, int SubGridSize, int OLD,int NEW, int my_id);
void handleCollisions(int ***Matrix, int SubGridSize, int OLD);
void printGrid(int ***Matrix, int GridSize, int layer);
void saveGridToFile(int ***Matrix, int GridSize, int layer, const char* filename);
void share_edges(int my_id, int ***Matrix, int SubGridSize, int NEW, int above, int below, int left, int right,
                 int *top_edge_roll_over, int *bottom_edge_roll_over, int *left_edge, int *right_edge);


int main(int argc, char** argv) {

    printf("\n\n--------- Environmet settings:---------- \n");
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

    // -------------- Initialize the OpenMP Environment:
    int num_of_threads = 1;
    //num_of_threads = omp_get_max_threads();
    printf("- OpenMP max threads %i \n",num_of_threads);
    omp_set_num_threads(num_of_threads);
    printf("--------------- finished -----------------\n\n");
    // -------------------------------------------------

    // allocate Vectors for send recv edges:
    int subGridSize  = GRID_SIZE / num_procs;
    int*** SubMatrix = NULL;
    //create matrix buffer:
    create_matrix(&SubMatrix,subGridSize,subGridSize,2);
    printf("created SubGrid processor %s ... \n",processor_name);
    // Allocate the edge buffers
    int *top_edge_roll_over = NULL, *bottom_edge_roll_over= NULL;
    int *left_edge= NULL, *right_edge= NULL;


    create_vector(&top_edge_roll_over, subGridSize);
    create_vector(&bottom_edge_roll_over, subGridSize);
    create_vector(&left_edge, subGridSize);
    create_vector(&right_edge, subGridSize);

    // allocate the global matrix
    int ***GlobalMatrix = NULL;  //
    if (my_id == 0) {
        // Allocation and initialization of GlobalMatrix
        create_matrix(&GlobalMatrix, subGridSize*processorGridSize, subGridSize*processorGridSize, 2);
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

    // Print off a hello world message
    printf("Hello from processor %s, rank %d out of %d processors\n "
           "--------my neighbours are ---------:\n"
           "above: %d\n"
           "right: %d\n"
           "below: %d\n"
           "left:  %d\n\n\n", processor_name, my_id, num_procs, above, right, below, left);
    //---------------------------------------------------------------------------------


    initializeGrid(SubMatrix, subGridSize); // fill grid with random numbers
    SubMatrix[int(subGridSize/2)][int(subGridSize/2)][0] = W; // setting an initial particle
    SubMatrix[int(subGridSize/2)][0][0] = E; // setting an initial particle

    // print the initial grid to console & saving it
    printf("Initial grid layer 0:\n");
    printGrid(SubMatrix, subGridSize, 0);
    printf("Initial grid layer 1:\n");
    printGrid(SubMatrix, subGridSize, 1);

    // saveGridToFile(grid, "grid_0.txt");
    int val1, val2 = 0;

    // sending all subgrids to process 0 so it can fill the main grid with it:
    if (my_id == 0) {
        // schreibe mein eigenes subgrid in Global Matrix (process 0's subgrid)
        int i,j;
#pragma omp parallel for private(j)
        for (i = 0; i < subGridSize; ++i) {
            for (j = 0; j < subGridSize; ++j) {
                GlobalMatrix[i][j][0] = SubMatrix[i][j][0];
            }
        }

        for (int p = 1; p < num_procs; ++p) {
            // get the starting positions for the process in the Global Matrix to store its Submatrix
            int start_row = (p / processorGridSize) * subGridSize;
            int start_col = (p % processorGridSize) * subGridSize;
            // receiving submatrices from processes
            MPI_Recv(&(GlobalMatrix[start_row][start_col][0]), subGridSize * subGridSize, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // print resulting grid to console:
        printf("starting Grid to save: \n");
        printGrid(GlobalMatrix, subGridSize*processorGridSize, 0);  // layer 1 will be the result !!!!
        char filename[20]; // erzeuge array of chars um Dateinamen zu speichern
        sprintf(filename, "starting_grid.txt"); // generiere Dateiname und speichere ihn ab
        saveGridToFile(GlobalMatrix,subGridSize*processorGridSize,0, filename); // erzeuge eine .txt datei mit Dateiname und speichere darin das grid des aktuellen standes ab

    }
    else {
        // send your submatrix
        MPI_Send(&(SubMatrix[0][0][0]), subGridSize * subGridSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    // gehe sicher, dass alle Prozessoren das gemacht haben bevor Sie weiter machen
    MPI_Barrier(MPI_COMM_WORLD);

    // ------------------------------------------------ start the itteration of the grid: -----------------------------------------------

    for (int step = 0; step < numb_iterations; ++step) {
        int val1 = step % 2; // alternate depending on the step
        int val2 = 1 - val1; // inverse of val1
        /*  step = 0, val1 = 0, val2 = 1
            step = 1, val1 = 1, val2 = 0
            step = 2, val1 = 0, val2 = 1
            step = 3, val1 = 1, val2 = 0
            step = 4, val1 = 0, val2 = 1
         */
        printf("value1: %i,  value2: %i \n",val1,val2);
        // CAVE: OLD grid needs to start @ layer 0 -> see initialize Grid
        handleCollisions(SubMatrix, subGridSize, val1); // check for collions FIRST and change directions of particles if needed
        moveParticles(SubMatrix, subGridSize, val1, val2, my_id);  // bewege alle Partikel entlang der Richtung
        //  // Teilen der Randwerte mit den benachbarten Prozessoren
        // share_edges(my_id, SubMatrix, subGridSize, val2, above, below, left, right, top_edge_roll_over, bottom_edge_roll_over,
        //             left_edge, right_edge);
        // print resulting grid to console:
        printf("New Grid after step %d:\n", step + 1);
        printGrid(SubMatrix, subGridSize, val2);
    }

    // sending all subgrids to process 0 after finishing itterations so it can fill the main grid with it:
    if (my_id == 0) {
        // schreibe mein eigenes subgrid in Global Matrix (process 0's subgrid)
        int i,j;
#pragma omp parallel for private(j)
        for (i = 0; i < subGridSize; ++i) {
            for (j = 0; j < subGridSize; ++j) {
                GlobalMatrix[i][j][1] = SubMatrix[i][j][val2];
            }
        }

        for (int p = 1; p < num_procs; ++p) {
            // get the starting positions for the process in the Global Matrix to store its Submatrix
            int start_row = (p / processorGridSize) * subGridSize;
            int start_col = (p % processorGridSize) * subGridSize;
            // receiving submatrices from processes
            MPI_Recv(&(GlobalMatrix[start_row][start_col][1]), subGridSize * subGridSize, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // print resulting grid to console:
        printf("resulting Grid: \n");
        printGrid(GlobalMatrix, subGridSize*processorGridSize, 1);  // layer 1 will be the result !!!!
        char filename[20]; // erzeuge array of chars um Dateinamen zu speichern
        sprintf(filename, "resulting_grid.txt"); // generiere Dateiname und speichere ihn ab
        saveGridToFile(GlobalMatrix,subGridSize*processorGridSize,1, filename); // erzeuge eine .txt datei mit Dateiname und speichere darin das grid des aktuellen standes ab

    }
    else {
        // send your submatrix
        MPI_Send(&(SubMatrix[0][0][val2]), subGridSize * subGridSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }



    // Free allocated memory:
    destroy_matrix(&SubMatrix,subGridSize, 2);
    printf("SubGrid of processor %s succesfully destroyed... \n",processor_name);
    if (my_id == 0){
        destroy_matrix(&GlobalMatrix, subGridSize*processorGridSize, 2);
        printf("Global Matrix of processor %s succesfully destroyed... \n",processor_name);
    }
    destroy_vector(&top_edge_roll_over);
    destroy_vector(&bottom_edge_roll_over);
    destroy_vector(&left_edge);
    destroy_vector(&right_edge);
    printf("vectors of processor %s succesfully destroyed... \n",processor_name);

    return 0;
}








void initializeGrid(int ***Matrix, int SubGridSize) {
    int rand_num = 0;
    for (int i = 0; i < SubGridSize; ++i) {
        for (int j = 0; j < SubGridSize; ++j) {
            rand_num = rand() % 16; // Zufälliger Zustand (0 bis 15) für alle möglichen Binärzustände, siehe Zustandsübergangstabelle
            Matrix[i][j][0] = 0; //rand_num;
            Matrix[i][j][1] = 0; //rand_num;
        }
    }


}

void moveParticles(int ***Matrix, int SubGridSize, int OLD,int NEW, int my_id) {
    // newGrid, and oldGrid is alternating -> layer 0 and 1
    // gehe durch das gesamte alte grid und wende die Zustandsübergangstabelle an
    int i, j;
# pragma omp parallel for private (j)
    for (i = 0; i < SubGridSize; ++i) {
        for (j = 0; j < SubGridSize; ++j) {
            if (Matrix[i][j][OLD] & N) {
                if (i > 0) Matrix[i-1][j][NEW] |= N; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (i == 0) Matrix[i+1][j][NEW] |= S; // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~N; // wenn Partikel bewegt, lösche alte Position
            }

            if (Matrix[i][j][OLD] & S) {
                if (i < SubGridSize - 1) Matrix[i+1][j][NEW] |= S; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (i == SubGridSize - 1) Matrix[i-1][j][NEW] |= N; // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~S;  // wenn Partikel bewegt, lösche alte Position
            }

            if (Matrix[i][j][OLD] & W) {
                if (j > 0) Matrix[i][j-1][NEW] |= W; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (j == 0) Matrix[i][j+1][NEW] |= E;  // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~W;  // wenn Partikel bewegt, lösche alte Position
            }

            if (Matrix[i][j][OLD] & E) {
                if (j < SubGridSize - 1) Matrix[i][j+1][NEW] |= E; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (j == SubGridSize - 1) Matrix[i][j-1][NEW] |= W; // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~E;  // wenn Partikel bewegt, lösche alte Position
            }
        }
    }

}

void handleCollisions(int ***Matrix, int SubGridSize, int OLD) {
    /*handels Collisions in the new grid*/
    int i, j;
# pragma omp parallel for private (j)
    for (i = 0; i < SubGridSize; ++i) {
        for (j = 0; j < SubGridSize; ++j) {
            int cell = Matrix[i][j][OLD];
            if ((cell & N) && (cell & S)) { // wenn in der Zelle Nord- und Südpartikel vorhanden sind
                cell &= ~(N | S); // entferne Sie die Nord- und Südpartikel
                cell |= (W | E);  // fügen West- und Ostpartikel hinzu
                printf("!N & S Collosion!");

            }

            else if((cell & W) && (cell & E)) { // wenn in der Zelle West und Ostpartikel vorhanden sind
                cell &= ~(W | E); // entfernen Sie die West- und Ostpartikel
                cell |= (N | S);  // fügen Nord- und Südpartikel hinzu
                printf("!W & E Collosion!");
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

void create_matrix(int ****Matrix, int nrows, int ncols, int nlayers)
{
    // Check, if *A=NULL. Empty A is necessary!
    assert ( *Matrix == NULL );

    // allocate the necessary memory

    *Matrix = (int***) malloc(nrows * sizeof(int**));

    int i, j;
# pragma omp parallel for private (j)
    for(i = 0; i < nrows; i++){
        (*Matrix)[i] = (int**) malloc(ncols * sizeof(int*));
        for(j = 0; j < ncols; j++){
            (*Matrix)[i][j] = (int*) malloc(nlayers * sizeof(int));
        }
    }

    // Check if the memory allocation was successful
    if (*Matrix == NULL) {
        perror("Failed to allocate memory for matrix");
        exit(EXIT_FAILURE);
    }
    return;
}



void destroy_matrix(int ****Matrix, int nrows, int nlayers) {
    // Check, if *A != NULL. non-empty A is necessary!
    assert(*Matrix != NULL);   // If condition is FALSE, program aborts!

    // free the memory
    int i, j;
# pragma omp parallel for private (j)
    for (i = 0; i < nlayers; ++i) {
        for (j = 0; j < nrows; ++j) {
            free((*Matrix)[i][j]);
        }
        free((*Matrix)[i]);
    }
    free(*Matrix);

    *Matrix = NULL;  // Set the pointer to NULL to avoid dangling pointer issues
    return;
}

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


void share_edges(int my_id, int ***Matrix, int SubGridSize, int NEW, int above, int below, int left, int right,
                 int *top_edge_roll_over, int *bottom_edge_roll_over, int *left_edge, int *right_edge) {

    // -------------------------------- Ränder übergeben -----------------------------
    // entnehme oberen und unteren Rand aus der Matrix
# pragma omp parallel for
    for (int i = 0; i < SubGridSize; ++i) {
        top_edge_roll_over[i] = Matrix[0][i][NEW]; // erste Zeile alle Spalten von NEW
        bottom_edge_roll_over[i] = Matrix[SubGridSize - 1][i][NEW]; // letzte Zeile alle Spalten von NEW
    }

    // Alle Subgrids reichen ihren oberen rand weiter and das Subgrid über ihnen:
    // source my_id: Subgrid das seinen oberen Rand weitereicht
    // destination: Subgrid das den oberen Rand von dem Subgrid unter Ihm erhält
    MPI_Sendrecv_replace(top_edge_roll_over, SubGridSize, MPI_INT, above, 0, my_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Der Processor hat den oberen Rand von dem Subgrid unter Ihm empfangen und speichert ihn in die unterste Zeile
    // damit wird top_edge_roll_over ersetze durch den oberen Rand des subgrids unter uns -> unser unterer Rand
    // Kopiere empfangene Daten in die Matrix
# pragma omp parallel for
    for (int i = 0; i < SubGridSize; ++i) {
        Matrix[SubGridSize - 1][i][NEW] = top_edge_roll_over[i]; // unterer Rand
    }

    // Alle Subgrids reichen ihren unteren rand weiter and das Subgrid unter ihnen:
    // source my_id: Subgrid das seinen unteren Rand weiterreicht
    // destination below: Subgrid das den unteren Rand von dem Subgrid über Ihm erhält
    MPI_Sendrecv_replace(bottom_edge_roll_over, SubGridSize, MPI_INT, below, 0, my_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Der Processor hat den unteren Rand von dem Subgrid über Ihm empfangen und speichert ihn in die oberste Zeile
    // damit wird bottom_edge_roll_over ersetze durch den unteren Rand des subgrids über uns -> unser neuer oberer Rand
    // Kopiere empfangene Daten in die Matrix
# pragma omp parallel for
    for (int i = 0; i < SubGridSize; ++i) {
        Matrix[0][i][NEW] = bottom_edge_roll_over[i]; // oberer Rand
    }


# pragma omp parallel for
    for (int i = 0; i < SubGridSize; ++i) {
        left_edge[i] = Matrix[i][0][NEW];
        right_edge[i] = Matrix[i][SubGridSize - 1][NEW];
    }

    // Sende linken Rand und empfange linken Rand von rechtem Prozessor -> damit unser neuer rechter Rand
    // Alle Subgrids reichen ihren linken Rand weiter and das Subgrid rechts neben Ihnen:
    // source my_id: Subgrid das seinen linken Rand weiterreicht
    // destination below: Subgrid das den linken Rand von dem Subgrid rechts neben Ihm erhält
    MPI_Sendrecv_replace(left_edge, SubGridSize, MPI_INT, left, 0, my_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // damit wird left_edge ersetze durch den linken rand des subgrids rechts neben uns -> unser neuer rechter Rand
    // Kopiere empfangene Daten in die Matrix
# pragma omp parallel for
    for (int i = 0; i < SubGridSize; ++i) {
        Matrix[i][SubGridSize - 1][NEW] = left_edge[i]; // rechter Rand
    }

    // Sende rechten Rand und empfange rechten Rand von linkem Prozessor -> damit unser neuer linker Rand
    // Alle Subgrids reichen ihren rechten Rand weiter and das Subgrid rechts neben Ihnen:
    // source my_id: Subgrid das seinen rechten Rand weiterreicht
    // destination below: Subgrid das den rechten Rand von dem Subgrid links neben ihm erhält
    MPI_Sendrecv_replace(right_edge, SubGridSize, MPI_INT, right, 0, my_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // damit wird right_edge ersetze durch den rechten rand des subgrids links neben uns -> unser neuer linker Rand
    // Kopiere empfangene Daten in die Matrix
# pragma omp parallel for
    for (int i = 0; i < SubGridSize; ++i) {
        Matrix[i][0][NEW] = right_edge[i]; // linker Rand
    }
}
