#include <iostream>
#include <termios.h>
#include <unistd.h>
#include <ctime>

// Dummy Data Initialization
void DummyDataInitialization(double* pMatrix, double* pVector, int Size) {
    for (int i = 0; i < Size; i++) {
        pVector[i] = 1;
        for (int j = 0; j < Size; j++)
            pMatrix[i * Size + j] = i;
    }
}

// Random Data Initialization
void RandomDataInitialization(double* pMatrix, double* pVector, int Size) {
    srand(unsigned(clock()));
    for (int i = 0; i < Size; i++) {
        pVector[i] = rand() / double(1000);
        for (int j = 0; j < Size; j++)
            pMatrix[i * Size + j] = rand() / double(1000);
    }
}

// Process Initialization
void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult, int Size) {
    pMatrix = new double[Size * Size];
    pVector = new double[Size];
    pResult = new double[Size];

    DummyDataInitialization(pMatrix, pVector, Size);
}

// Print Matrix
void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
    for (int i = 0; i < RowCount; i++) {
        for (int j = 0; j < ColCount; j++)
            printf("%7.4f ", pMatrix[i * RowCount + j]);
        printf("\n");
    }
}

// Print Vector
void PrintVector(double* pVector, int Size) {
    for (int i = 0; i < Size; i++)
        printf("%7.4f ", pVector[i]);
    printf("\n");
}

// Result Calculation
void ResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size) {
    for (int i = 0; i < Size; i++) {
        pResult[i] = 0;
        for (int j = 0; j < Size; j++)
            pResult[i] += pMatrix[i * Size + j] * pVector[j];
    }
}

// Process Termination
void ProcessTermination(double* pMatrix, double* pVector, double* pResult) {
    delete[] pMatrix;
    delete[] pVector;
    delete[] pResult;
}

// Get character input (for future use if needed)
int getch() {
    struct termios oldt, newt;
    int ch;
    tcgetattr(STDIN_FILENO, &oldt);
    newt = oldt;
    newt.c_lflag &= ~(ICANON | ECHO);
    tcsetattr(STDIN_FILENO, TCSANOW, &newt);
    ch = getchar();
    tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
    return ch;
}

// Test runner function (already provided)
void runTest(int testNumber, int size) {
    double* pMatrix;
    double* pVector;
    double* pResult;

    clock_t start, finish;
    double duration;

    ProcessInitialization(pMatrix, pVector, pResult, size);

    start = clock();
    ResultCalculation(pMatrix, pVector, pResult, size);
    finish = clock();
    duration = (finish - start) / double(CLOCKS_PER_SEC);

    std::cout << "Test:" << testNumber << " | Size: " << size << " | Time: " << duration << " seconds\n";

    ProcessTermination(pMatrix, pVector, pResult);
}

// Main function to run tests
int main() {
    int testSizes[] = {10, 100, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
    int numTests = sizeof(testSizes) / sizeof(testSizes[0]);

    for (int i = 0; i < numTests; ++i) {
        runTest(i + 1, testSizes[i]);
    }

    return 0;
}
