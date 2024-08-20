#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Funktionsprototypen
char* generateStringFromOctalArray(int *arr, int length);
void convertStringToOctalArray(char *str, int **octalArray, int *arrayLength);
char* readStringFromFile(const char* filename);
void writeStringToFile(const char* filename, const char* content);
void print_vector(int *vector, int vectorsize);

int main(int argc, char** argv) {

    const char* inputFilename = "../Messages/original_message.txt";
    const char* outputFilename = "../Messages/decrypted_message.txt";

    // char message[]= "Hello, World!";
    // writeStringToFile("original_message.txt", message);



    char* inputMessage = readStringFromFile(inputFilename);
    if (inputMessage == NULL) {
        fprintf(stderr, "Fehler beim Lesen der Datei: %s\n", inputFilename);
        return 1;
    }

    printf("original message: %s\n", inputMessage);

    int *octal_array = NULL;
    int octal_array_length = 0;

    convertStringToOctalArray(inputMessage, &octal_array, &octal_array_length);

    printf("converted message in octal: ");
    print_vector(octal_array, octal_array_length);

    char *reconstructedMessage = generateStringFromOctalArray(octal_array, octal_array_length);
    printf("reconstructed message from octal: %s \n", reconstructedMessage);
    // Schreibe den rekonvertierten String in die Ausgabedatei
    writeStringToFile(outputFilename, reconstructedMessage);



    // Speicher freigeben
    free(inputMessage);
    free(octal_array);
    free(reconstructedMessage);

    return 0;
}

char* generateStringFromOctalArray(int *arr, int length) {
    if (length % 3 != 0) {
        printf("Das Array muss eine durch 3 teilbare Länge haben.\n");
        return NULL;
    }

    int resultLength = length / 3;
    char *result = (char*)malloc(resultLength + 1); // +1 for null terminator
    if (result == NULL) {
        printf("Speicher konnte nicht zugewiesen werden.\n");
        return NULL;
    }

    int index = 0;
    for (int i = 0; i < length; i += 3) {
        int octalValue = arr[i] * 64 + arr[i + 1] * 8 + arr[i + 2];
        result[index++] = (char)octalValue;
    }

    result[index] = '\0'; // Null terminator for the string

    return result;
}

void convertStringToOctalArray(char *str, int **octalArray, int *arrayLength) {
    int length = strlen(str);
    *octalArray = (int*)malloc(length * 3 * sizeof(int));
    if (*octalArray == NULL) {
        printf("Speicher konnte nicht zugewiesen werden.\n");
        return;
    }

    int index = 0;

    for (int i = 0; i < length; i++) {
        int asciiValue = (int)str[i];
        char octalString[4];
        sprintf(octalString, "%03o", asciiValue);

        (*octalArray)[index++] = octalString[0] - '0';
        (*octalArray)[index++] = octalString[1] - '0';
        (*octalArray)[index++] = octalString[2] - '0';
    }

    *arrayLength = index;
}

void print_vector(int *vector, int vectorsize) {
    for (int i = 0; i < vectorsize; ++i) {
        printf("%02d ", vector[i]);
    }
    printf("\n");
}

char* readStringFromFile(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Fehler beim Öffnen der Datei");
        return NULL;
    }

    fseek(file, 0, SEEK_END);
    long fileSize = ftell(file);
    fseek(file, 0, SEEK_SET);

    char* content = (char*)malloc(fileSize + 1);
    if (content == NULL) {
        perror("Fehler bei der Speicherzuweisung");
        fclose(file);
        return NULL;
    }

    fread(content, 1, fileSize, file);
    content[fileSize] = '\0';

    fclose(file);
    return content;
}

void writeStringToFile(const char* filename, const char* content) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("Fehler beim Öffnen der Datei");
        return;
    }

    fprintf(file, "%s", content);

    fclose(file);
}
