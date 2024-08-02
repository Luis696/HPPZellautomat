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


int num_procs = 16;

int main(){
    int processorGridSize = sqrt(num_procs);
    for(int my_id = 0; my_id<num_procs; my_id++){


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

        // Print off a hello world message
        printf("Hello from processor, rank %d out of %d processors\n "
               "--------my neighbours are ---------:\n"
               "above: %d\n"
               "right: %d\n"
               "below: %d\n"
               "left:  %d\n"
               "edge:  %i\n\n\n", my_id, num_procs, above, right, below, left, edge);
    }
}