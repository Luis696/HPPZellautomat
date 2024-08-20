#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <cmath>
#include <assert.h>
#include <string.h>
//#define GRID_SIZE 0 // lege die größe des grids fest
#define numb_iterations 43
#define collosion_on true
#define debug_mode false
// Bits für die Partikelrichtungen (siehe Zustandsübergangstabelle)
#define N 2 // 0010
#define S 8 // 1000
#define W 1 // 0001
#define E 4 // 0100
char filename[250]; // Puffer für den Dateinamen, ausreichend groß
char filename_message[250];
char path[] = "Grids/";
char image_name [] = "resizedImage.txt";
/*
* 0000 0: Kein Partikel +
* 0001 1: Partikel geht nach Westen (W) +
* 0010 2: Partikel geht nach Norden (N) +
* 0011 3: Partikel gehen nach Norden und Westen (N & W) +
* 0100 4: Partikel geht nach Osten (E) +
* 0101 5: Partikel gehen nach Osten und Westen (E & W) KOLLOSION +
* 0110 6: Partikel gehen nach Norden und Osten (N & E) +
* 0111 7: Partikel gehen nach Norden, Osten und Westen (N & E & W) + maybe also a Kollosion ?
* 1000 8: Partikel geht nach Süden (S) +
* 1001 9: Partikel gehen nach Süden und Westen (S & W) +
* 1010 10: Partikel gehen nach Süden und Norden (S & N) KOLLOSION +
* 1011 11: Partikel gehen nach Süden, Norden und Westen (S & N & W) + maybe also a Kollosion ?
* 1100 12: Partikel gehen nach Süden und Osten (S & E) +
* 1101 13: Partikel gehen nach Süden, Osten und Westen (S & E & W) + maybe also a Kollosion ?
* 1110 14: Partikel gehen nach Süden, Norden und Osten (S & N & E) + maybe also a Kollosion ?
* 1111 15: Partikel gehen in alle Richtungen (S & N & E & W) +
*
* Kollisionen: 5, 10
* 5 = Kollision zwischen Osten und Westen (E & W)
* 10 = Kollision zwischen Süden und Norden (S & N)
*/

// definiere alle Funktionen die wir brauchen
void create_matrix(int * ***Matrix, int nrows, int ncols, int nlayers);
void destroy_matrix(int * ***Matrix, int nrows, int nlayers);
void create_vector(int **vector, int SubGridSize);
void destroy_vector(int **vector);
void load_matrix_from_file(int*** Matrix, const char* filename, int rowCount, int colCount, int layer);
void initializeGrid(int ***Matrix, int nrows, int , int nlayers,int initializing_value);
void moveParticles(int ***Matrix, int SubGridSize, int OLD,int NEW, int my_id, int edge);
void handleCollisions(int ***Matrix, int SubGridSize, int OLD);
void printGrid(int ***Matrix, int GridSize, int layer);
void saveGridToFile(int ***Matrix, int GridSize, int layer, const char* filename);
void share_edges(int my_id, int ***Matrix, int SubGridSize, int NEW,int OLD, int above, int below, int left, int right,
                 int *top_edge_roll_over, int *bottom_edge_roll_over, int *left_edge, int *right_edge, int edge);
void gatherSubgrids(int ***GlobalMatrix, int *** BufferMatrix, int *** ownSubMatrix, int subGridSize, int subGridLayers, int MainMatrixsize , int num_procs, int processorGridSize, const char *filename);
void flipDirections(int ***Matrix, int nrows, int ncols, int nlayers);
void distributeSubgrids(int ***GlobalMatrix, int *** BufferMatrix, int *** ownSubMatrix,  int subGridSize, int subGridLayers, int MainMatrixsize , int num_procs, int processorGridSize, const char *filename);
void print_vector(int vector, int vectorsize);



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


    int subGridSize = 3;
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

    int ***GlobalMatrix = NULL;  //
    int ***BufferMatrix = NULL;

    if(my_id == 0) {
        // print information about Grids sizes
        printf("Number of Processors:%i \n", num_procs);
        printf("Processor Gridzize: %i x %i \n", processorGridSize, processorGridSize);
        printf("Subgridsize: %i x %i \n", subGridSize, subGridSize);
        printf("Main Matrix size: %i x %i \n\n", MainMatrixsize, MainMatrixsize);


        // allocate the global matrix for processor 0
        create_matrix(&GlobalMatrix, MainMatrixsize, MainMatrixsize, subGridLayers);
        initializeGrid(GlobalMatrix, MainMatrixsize, MainMatrixsize,subGridLayers, 0);
        // load_matrix_from_file(GlobalMatrix, image_name,MainMatrixsize,MainMatrixsize, 2);

        // // N & S & E & W @ edge = W +
        // GlobalMatrix[0][0][0] |= S; // setting an initial particle
        // GlobalMatrix[0][0][0] |= N; // setting an initial particle
        // GlobalMatrix[0][0][0] |= E; // setting an initial particle
        // GlobalMatrix[0][0][0] |= W; // setting an initial particle


        // set random initial particles
        int particletypes[4] = {N, S,  W, E};
        for(int p = 0; p < 80; p++) {
            int x = rand() % MainMatrixsize; // Get a random x-coordinate.
            int y = rand() % MainMatrixsize; // Get a random y-coordinate.
            int particle = rand() % 4; // Select a random particle.
            GlobalMatrix[x][y][0] = particletypes[particle];

        }

        // Allocate and initialize the temporary (buffer) matrix.
        create_matrix(&BufferMatrix, subGridSize, subGridSize, subGridLayers);

        // collecting all Subgrids from all processors
        distributeSubgrids(GlobalMatrix, BufferMatrix, SubMatrix, subGridSize,subGridLayers, MainMatrixsize, num_procs, processorGridSize, filename);

        int step = 0;
        sprintf(filename, "%sencrypting_grid_%i.txt",path, step);
        printf("%s: \n", filename);
        printf("---------------------- intial Grid ---------------\n");
        printGrid(GlobalMatrix,MainMatrixsize,0);
        // save grid with particles, beginning state:
        saveGridToFile(GlobalMatrix, MainMatrixsize, 0, filename);

    }

    if (my_id != 0) {
        // send your submatrix
         MPI_Recv(&(SubMatrix[0][0][0]), subGridSize * subGridSize * subGridLayers, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // gehe sicher, dass alle Prozessoren das gemacht haben bevor Sie weiter machen
    MPI_Barrier(MPI_COMM_WORLD);
    //
    int new_layer = 0;
    int old_layer = 1;
    // ------------------------------------------------ start the Encryption itteration of the grid: -----------------------------------------------
    if(my_id == 0){printf("---------------------- start Encryption ---------------\n");}
    for (int step = 1; step < numb_iterations; ++step) {
        {int swap=new_layer; new_layer = old_layer; old_layer = swap;}
        if(debug_mode) {
            printf("NEW: %i \n",new_layer);
            printf("OLD: %i \n", old_layer);
        }
        /*  step = 0, val1 = 0, val2 = 1
            step = 1, val1 = 1, val2 = 0
            step = 2, val1 = 0, val2 = 1
            step = 3, val1 = 1, val2 = 0
            step = 4, val1 = 0, val2 = 1
         */
        // printf("value1: %i,  value2: %i \n",val1,val2);
        // CAVE: OLD grid needs to start @ layer 0 -> see initialize Grids
        if(collosion_on) {
            handleCollisions(SubMatrix, subGridSize, old_layer);
        }// check for collions FIRST and change directions of particles if needed
        // Teilen der Randwerte mit den benachbarten Prozessoren
        share_edges(my_id, SubMatrix, subGridSize, new_layer,old_layer, above, below, left, right, top_edge_roll_over, bottom_edge_roll_over,
                    left_edge, right_edge, edge);

        moveParticles(SubMatrix, subGridSize, old_layer, new_layer, my_id, edge);  // bewege alle Partikel entlang der Richtung

        // sending all subgrids to process 0 so it can fill the main grid with it:
        if (my_id == 0 & debug_mode) {
            gatherSubgrids(GlobalMatrix, BufferMatrix, SubMatrix, subGridSize,subGridLayers, MainMatrixsize, num_procs, processorGridSize, filename);
            sprintf(filename, "%sencrypting_grid_%i.txt",path, step);
            printf("%s: \n", filename);
            printGrid(GlobalMatrix,MainMatrixsize,new_layer);

            // save Grid with particles:
            sprintf(filename, "%sencrypting_grid_%i.txt",path, step+1);
            saveGridToFile(GlobalMatrix, MainMatrixsize, new_layer, filename);
        }
        if (my_id != 0 & debug_mode) {
            // send your submatrix
            MPI_Send(&(SubMatrix[0][0][0]), subGridSize * subGridSize * subGridLayers, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
    // ---------------------------------------------------- finished the Encryption of the grid: ----------------------------------
    printf("---------------------- finished Encryption ---------------\n");

    // ------------------------------------------------ start the Decryption itteration of the grid: -----------------------------------------------

    int Decrypt_Map[16] = {0, 4, 8, 12, 1, 10, 9, 13, 2, 6, 5, 14, 3, 7, 11, 15};

    if(debug_mode){printf("start flipping layer \n");}
    flipDirections(SubMatrix, subGridSize, subGridSize, 2);
    if(debug_mode){printf("flipped layers\n");}


    if (my_id == 0 & debug_mode) {
        printf("---------------------- intial Grid Decryption ---------------\n");
        // removing Particle Encryption by subtraction
        gatherSubgrids(GlobalMatrix, BufferMatrix, SubMatrix, subGridSize,subGridLayers, MainMatrixsize, num_procs, processorGridSize, filename);
        sprintf(filename, "%sdecrypting_grid_%i.txt",path, 0);
        printf("%s: \n", filename);
        printGrid(GlobalMatrix,MainMatrixsize,new_layer);
        // save grid with particles:
        saveGridToFile(GlobalMatrix, MainMatrixsize, new_layer, filename);

    }

    if (my_id != 0 & debug_mode) {
        // send your submatrix
        MPI_Send(&(SubMatrix[0][0][0]), subGridSize * subGridSize * subGridLayers, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    if(my_id == 0){printf("---------------------- start Decryption ---------------\n");}
    for (int step = 1; step < numb_iterations; ++step) {
        // Erzeuge den Dateinamen mit der Iterationsnummer
        {int swap=new_layer; new_layer = old_layer; old_layer = swap;}
        if(debug_mode) {
            printf("NEW: %i \n",new_layer);
            printf("OLD: %i \n", old_layer);
        }
        /*  step = 0, val1 = 0, val2 = 1
            step = 1, val1 = 1, val2 = 0
            step = 2, val1 = 0, val2 = 1
            step = 3, val1 = 1, val2 = 0
            step = 4, val1 = 0, val2 = 1
         */
        // printf("value1: %i,  value2: %i \n",val1,val2);
        // CAVE: OLD grid needs to start @ layer 0 -> see initialize Grids
        if(collosion_on) {
            handleCollisions(SubMatrix, subGridSize, old_layer);
        }// check for collions FIRST and change directions of particles if needed
        // Teilen der Randwerte mit den benachbarten Prozessoren
        share_edges(my_id, SubMatrix, subGridSize, new_layer,old_layer, above, below, left, right, top_edge_roll_over, bottom_edge_roll_over,
                    left_edge, right_edge, edge);

        moveParticles(SubMatrix, subGridSize, old_layer, new_layer, my_id, edge);  // bewege alle Partikel entlang der Richtung


        // sending all subgrids to process 0 so it can fill the main grid with it:
        if (my_id == 0 & debug_mode) {
            gatherSubgrids(GlobalMatrix, BufferMatrix, SubMatrix, subGridSize,subGridLayers, MainMatrixsize, num_procs, processorGridSize, filename);
            sprintf(filename, "%sdecrypting_grid_%i.txt",path, step);
            printf("%s: \n", filename);
            printGrid(GlobalMatrix,MainMatrixsize,new_layer);

            // save grid with particles:
            saveGridToFile(GlobalMatrix, MainMatrixsize, new_layer, filename);
        }
        if (my_id != 0 & debug_mode) {
            // send your submatrix
            MPI_Send(&(SubMatrix[0][0][0]), subGridSize * subGridSize * subGridLayers, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
    printf("---------------------- finished Decryption ---------------\n");
    // ---------------------------------------------------- finished the Decryption of the grid: ----------------------------------


    if (my_id == 0) {
        gatherSubgrids(GlobalMatrix, BufferMatrix, SubMatrix, subGridSize,subGridLayers, MainMatrixsize, num_procs, processorGridSize, filename);
        sprintf(filename, "%sdecrypting_grid_final.txt",path);
        printf("%s: \n", filename);
        printGrid(GlobalMatrix,MainMatrixsize,new_layer);
        // save grid with particles:
        saveGridToFile(GlobalMatrix, MainMatrixsize, new_layer, filename);
    }
    if (my_id != 0) {
        // send your submatrix
        MPI_Send(&(SubMatrix[0][0][0]), subGridSize * subGridSize * subGridLayers, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }


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
    MPI_Finalize();
    return 0;

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
                if (i == 0 && (edge & N)) Matrix[i][j][NEW] |= S; // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~N; // wenn Partikel bewegt, lösche alte Position
            }

            if (Matrix[i][j][OLD] & S) {
                if (i < SubGridSize - 1) Matrix[i+1][j][NEW] |= S; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (i == SubGridSize - 1 && (edge & S)) Matrix[i][j][NEW] |= N; // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~S;  // wenn Partikel bewegt, lösche alte Position
            }

            if (Matrix[i][j][OLD] & W) {
                if (j > 0) Matrix[i][j-1][NEW] |= W; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (j == 0 && (edge & W)) Matrix[i][j][NEW] |= E;  // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
                Matrix[i][j][OLD] &= ~W;  // wenn Partikel bewegt, lösche alte Position
            }

            if (Matrix[i][j][OLD] & E) {
                if (j < SubGridSize - 1) Matrix[i][j+1][NEW] |= E; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
                if (j == SubGridSize - 1 && (edge & E)) Matrix[i][j][NEW] |= W; // if hitting a frame bounce back -> ONLY for the Submatrixes which are the edges of the main Matrix
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
            if((cell == 5) || (cell == 10)) {
                if ((cell & N) && (cell & S)) { // wenn in der Zelle Nord- und Südpartikel vorhanden sind
                    cell &= ~(N | S); // entferne die Nord- und Südpartikel
                    cell |= (W | E);  // fügen West- und Ostpartikel hinzu
                    //printf("!N & S Collosion!");

                }

                else if((cell & W) && (cell & E)) { // wenn in der Zelle West und Ostpartikel vorhanden sind
                    cell &= ~(W | E); // entferne die West- und Ostpartikel
                    cell |= (N | S);  // fügen Nord- und Südpartikel hinzu
                    //printf("!W & E Collosion!");
                }
                Matrix[i][j][OLD] = cell;
            }
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
    // printf("starting to save grid to file...\n");
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
    // printf("saved grid to file\n\n");
}


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
void create_vector(int **vector, int length) {
    // Check, if *vector == NULL. Non-empty vector is not allowed!
    assert(*vector == NULL); // If condition is FALSE, program aborts!

    // Allocate the necessary memory
    *vector = (int*) malloc(length * sizeof(int));

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


void print_vector(int *vector, int vectorsize) {
    for (int i = 0; i < vectorsize; ++i) {
            printf("%02d ", vector[i]);
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
        if (Matrix[0][i][OLD] & N && !(edge & N)) { // if particle is moving north and if subgrid has not northern edge
            top_edge_roll_over[i] = N;
            Matrix[0][i][OLD] &=~N; // wenn Partikel übergeben, entferne es aus der Matrix
        }
        else{top_edge_roll_over[i] = 0;}

        if (Matrix[SubGridSize - 1][i][OLD] & S && !(edge & S)) { // if particle is moving south and if subgrid has not southern edge
            bottom_edge_roll_over[i] = S;
            Matrix[SubGridSize - 1][i][OLD] &= ~S; // wenn Partikel übergeben, entferne es aus der Matrix
        }
        else{bottom_edge_roll_over[i] = 0;}

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

        if (Matrix[i][0][OLD] & W && !(edge & W)) { // checke if partikel auch wirklich nach Westen unterwegs und ob es sich nicht um einen Rand handelt
            left_edge[i] = W;
            Matrix[i][0][OLD] &= ~W; // wenn Partikel übergeben, entferne es aus der Matrix
        }
        else {left_edge[i] = 0;}


        if (Matrix[i][SubGridSize - 1][OLD] & E && !(edge & E)) {
            right_edge[i] = E;
            Matrix[i][SubGridSize - 1][OLD] &= ~E; // wenn Partikel übergeben, entferne es aus der Matrix
        }
        else {right_edge[i] = 0;}


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

void distributeSubgrids(int ***GlobalMatrix, int *** BufferMatrix, int *** ownSubMatrix,  int subGridSize, int subGridLayers, int MainMatrixsize , int num_procs, int processorGridSize, const char *filename) {
    // Assign root processors own SubMatrix from Global Matrix
    for (int i = 0; i < subGridSize; ++i) {
        for (int j = 0; j < subGridSize; ++j) {
            for(int k = 0; k<subGridLayers; k++){
                ownSubMatrix[i][j][k] = GlobalMatrix[i][j][k];
            }
        }
    }

    // Assign SubMatrices to all other processors
    for (int p = 1; p < num_procs; ++p) {
        // Determine the starting coordinates for the submatrix in the global matrix.
        int start_row = (p / processorGridSize) * subGridSize;
        int start_col = (p % processorGridSize) * subGridSize;

        // Assign the submatrix values from the global matrix to the buffer matrix.
        for (int i = 0; i < subGridSize; ++i) {
            for (int j = 0; j < subGridSize; ++j) {
                for(int k = 0; k<subGridLayers; k++){
                    BufferMatrix[i][j][k] = GlobalMatrix[start_row+i][start_col+j][k];
                }
            }
        }

        // Send buffer matrices to the other processes.
        MPI_Send(&(BufferMatrix[0][0][0]), subGridSize * subGridSize * subGridLayers, MPI_INT, p, 0, MPI_COMM_WORLD);
    }
}

void load_matrix_from_file(int*** Matrix, const char* filename, int rowCount, int colCount, int layer)
{
    // Open the file
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Could not open file %s\n", filename);
        return;
    }

    // Read matrix from file
    for (int i = 0; i < rowCount; ++i)
        for (int j = 0; j < colCount; ++j)
            fscanf(file, "%d", &Matrix[i][j][layer]);

    fclose(file);
}


void flipDirections(int ***Matrix, int nrows, int ncols, int nlayers) {
    /* FlipMap Lookup Table
    * Jede Partikelrichtung wird zu ihrer entgegengesetzten Richtung gewechselt.
    */
    // int flip_map[16] = {0, 4, 8, 12, 1, 10, 9, 11, 2, 6, 5, 7, 3, 7, 6, 0};
    int flip_map[16] = {0, 4, 8, 12, 1, 10, 9, 13, 2, 6, 5, 14, 3, 7, 11, 15};

    for (int layer = 0; layer < nlayers; ++layer) {
        for (int row = 0; row < nrows; ++row) {
            for (int col = 0; col < ncols; ++col) {
                Matrix[row][col][layer] = flip_map[Matrix[row][col][layer]];
            }
        }
    }
}

