#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

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
*/

// definiere alle Funktionen die wir brauchen
void initializeGrid(int grid[GRID_SIZE][GRID_SIZE]);
void moveParticles(int grid[GRID_SIZE][GRID_SIZE]);
void handleCollisions(int grid[GRID_SIZE][GRID_SIZE]);
void printGrid(int grid[GRID_SIZE][GRID_SIZE]);
void saveGridToFile(int grid[GRID_SIZE][GRID_SIZE], const char* filename);

int main() {
    int num_of_threads = omp_get_max_threads();
    printf("omp max threads %i \n",num_of_threads);
    omp_set_num_threads(num_of_threads);


    int grid[GRID_SIZE][GRID_SIZE];  // assign grid
    initializeGrid(grid); // fill grid with random numbers

    // frint the initial grid to console & saving it
    printf("Initial grid:\n");
    printGrid(grid);
    saveGridToFile(grid, "grid_0.txt");

    for (int step = 0; step < numb_iterations; ++step) {
        moveParticles(grid);  // bewege alle Partikel entlang der Richtung
        handleCollisions(grid); // bearbeite kollidierende partikel N<->S, E<->W

         // print resulting grid to console:
        printf("Grids after step %d:\n", step + 1);
        printGrid(grid);
        // char filename[20]; // erzeuge array of chars um Dateinamen zu speichern
        // sprintf(filename, "grid_%d.txt", step + 1); // generiere Dateiname und speichere ihn ab
        // saveGridToFile(grid, filename); // erzeuge eine .txt datei mit Dateiname und speichere darin das grid des aktuellen standes ab
    }
    printf("resulting array: \n");
    printGrid(grid);
    return 0;
}

void initializeGrid(int grid[GRID_SIZE][GRID_SIZE]) {
    for (int i = 0; i < GRID_SIZE; ++i) {
        for (int j = 0; j < GRID_SIZE; ++j) {
            //grid[i][j] = rand() % 16; // Zufälliger Zustand (0 bis 15) für alle möglichen Binärzustände, siehe Zustandsübergangstabelle
            grid[i][j] = 0;
        }
    }
    grid[int(GRID_SIZE/2)][int(GRID_SIZE/2)]=N;
}

void moveParticles(int grid[GRID_SIZE][GRID_SIZE]) {
    int newGrid[GRID_SIZE][GRID_SIZE] = {0};  // erzeuge ein neues Grids das den nächsten schritt enthält

    // gehe durch das gesamte alte grid und wende die Zustandsübergangstabelle an
# pragma omp parallel for
    for (int i = 0; i < GRID_SIZE; ++i) {
        for (int j = 0; j < GRID_SIZE; ++j) {
            if (grid[i][j] & N && i > 0) newGrid[i-1][j] |= N; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
            if (grid[i][j] & N && i == 0) newGrid[i+1][j] |= S; // if hitting a frame flied bounce back

            if (grid[i][j] & S && i < GRID_SIZE - 1) newGrid[i+1][j] |= S; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
            if (grid[i][j] & S && i == GRID_SIZE - 1) newGrid[i-1][j] |= N; // if hitting a frame flied bounce back

            if (grid[i][j] & W && j > 0) newGrid[i][j-1] |= W; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
            if (grid[i][j] & W && j == 0) newGrid[i][j+1] |= E;  // if hitting a frame flied bounce back

            if (grid[i][j] & E && j < GRID_SIZE - 1) newGrid[i][j+1] |= E; // checke in welche Richtung das Partikel unterwegs ist und schiebe es weiter
            if (grid[i][j] & E && j == GRID_SIZE - 1) newGrid[i][j-1] |= W; // if hitting a frame flied bounce back

        }
    }
    // überschreibe das alte Grids mit den neuen Zellwerten
    for (int i = 0; i < GRID_SIZE; ++i) {
        for (int j = 0; j < GRID_SIZE; ++j) {
            grid[i][j] = newGrid[i][j];
        }
    }
}

void handleCollisions(int grid[GRID_SIZE][GRID_SIZE]) {
#pragma omp parallel for
    for (int i = 0; i < GRID_SIZE; ++i) {
        for (int j = 0; j < GRID_SIZE; ++j) {
            int cell = grid[i][j];
            if ((cell & N) && (cell & S)) { // wenn in der Zelle Nord- und Südpartikel vorhanden sind
                cell &= ~(N | S); // entferne Sie die Nord- und Südpartikel
                cell |= (W | E);  // fügen West- und Ostpartikel hinzu
            }
            if ((cell & W) && (cell & E)) { // wenn in der Zelle West und Ostpartikel vorhanden sind
                cell &= ~(W | E); // entfernen Sie die West- und Ostpartikel
                cell |= (N | S);  // fügen Nord- und Südpartikel hinzu
            }
            grid[i][j] = cell; // setze die aktuelle Zelle im grid auf den neuen Wert
        }
    }
}



void printGrid(int grid[GRID_SIZE][GRID_SIZE]) {
    for (int i = 0; i < GRID_SIZE; ++i) {
        for (int j = 0; j < GRID_SIZE; ++j) {
            printf("%02d ", grid[i][j]);
        }
        printf("\n");  // nach jeder ferigen Zeile beginn einen neuen Absatz
    }
    printf("\n");  // nachdem das gesamte grid ausgegeben wurde mache einen Absatz
}

void saveGridToFile(int grid[GRID_SIZE][GRID_SIZE], const char* filename) {
    FILE* file = fopen(filename, "w");  // erzeuge eine Datei und öffne sie
    if (!file) { // wenn öffnen nicht funktioniert, dann gebe mir einen Fehler aus
        perror("Failed to open file");
        return;
    }
    // sonst speichere das grid unter dem Dateinamen ab
    for (int i = 0; i < GRID_SIZE; ++i) {
        for (int j = 0; j < GRID_SIZE; ++j) {
            fprintf(file, "%02d ", grid[i][j]);  // schreibe grid wert in die Datei
        }
        fprintf(file, "\n");  // nachdem eine ganze Zeile geschrieben wurde beende die Line
    }
    fclose(file);  // schließe die Datei wieder
}
