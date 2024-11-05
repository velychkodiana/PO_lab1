#include <mpi.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>

int ProcNum = 0; // Кількість доступних процесів
int ProcRank = 0; // Ранг поточного процесу

// Ініціалізація значень матриці та вектора (випадкові дані)
void RandomDataInitialization(double* pMatrix, double* pVector, int Size) {
    srand(static_cast<unsigned>(time(nullptr)));
    for (int i = 0; i < Size; ++i) {
        pVector[i] = rand() / double(RAND_MAX);
        for (int j = 0; j < Size; ++j) {
            pMatrix[i * Size + j] = rand() / double(RAND_MAX);
        }
    }
}

// Ініціалізація та розподіл пам'яті
void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult,
                           double*& pProcRows, double*& pProcResult, int Size, int& RowNum) {
    int RestRows = Size;
    for (int i = 0; i < ProcRank; i++)
        RestRows -= RestRows / (ProcNum - i);
    RowNum = RestRows / (ProcNum - ProcRank);

    pVector = new double[Size];
    pResult = new double[Size];
    pProcRows = new double[RowNum * Size];
    pProcResult = new double[RowNum];

    if (ProcRank == 0) {
        pMatrix = new double[Size * Size];
        RandomDataInitialization(pMatrix, pVector, Size);
    }
}

// Функція для розподілу даних між процесами
void DataDistribution(double* pMatrix, double* pProcRows, double* pVector, int Size, int RowNum) {
    MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int *pSendInd = new int[ProcNum], *pSendNum = new int[ProcNum];
    int RestRows = Size;
    RowNum = Size / ProcNum;
    pSendNum[0] = RowNum * Size;
    pSendInd[0] = 0;

    for (int i = 1; i < ProcNum; i++) {
        RestRows -= RowNum;
        RowNum = RestRows / (ProcNum - i);
        pSendNum[i] = RowNum * Size;
        pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
    }

    MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows, pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    delete[] pSendInd;
    delete[] pSendNum;
}

// Функція для множення частини матриці на вектор
void ParallelResultCalculation(double* pProcRows, double* pVector, double* pProcResult, int Size, int RowNum) {
    for (int i = 0; i < RowNum; ++i) {
        pProcResult[i] = 0;
        for (int j = 0; j < Size; ++j) {
            pProcResult[i] += pProcRows[i * Size + j] * pVector[j];
        }
    }
}

// Збирання результатів
void ResultReplication(double* pProcResult, double* pResult, int Size, int RowNum) {
    int *pReceiveInd = new int[ProcNum], *pReceiveNum = new int[ProcNum];
    int RestRows = Size;
    RowNum = Size / ProcNum;
    pReceiveNum[0] = RowNum;
    pReceiveInd[0] = 0;

    for (int i = 1; i < ProcNum; i++) {
        RestRows -= RowNum;
        RowNum = RestRows / (ProcNum - i);
        pReceiveNum[i] = RowNum;
        pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
    }

    MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank], MPI_DOUBLE, pResult, pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);

    delete[] pReceiveInd;
    delete[] pReceiveNum;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    std::vector<int> sizes = {10, 100, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
    double SequentialTime = 0.0;

    for (size_t test = 0; test < sizes.size(); ++test) {
        int Size = sizes[test];
        double *pMatrix = nullptr, *pVector = nullptr, *pResult = nullptr;
        double *pProcRows = nullptr, *pProcResult = nullptr;
        int RowNum = 0;

        ProcessInitialization(pMatrix, pVector, pResult, pProcRows, pProcResult, Size, RowNum);
        // Виконуємо послідовний обчислення для процесу з рангом 0
        if (ProcRank == 0) {
            double SeqStart = MPI_Wtime();
            for (int i = 0; i < Size; ++i) {
                pResult[i] = 0;
                for (int j = 0; j < Size; ++j) {
                    pResult[i] += pMatrix[i * Size + j] * pVector[j];
                }
            }
            double SeqFinish = MPI_Wtime();
            SequentialTime = SeqFinish - SeqStart;
        }

        // Розсилаємо час послідовного виконання всім процесам
        MPI_Bcast(&SequentialTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        double Start = MPI_Wtime();
        DataDistribution(pMatrix, pProcRows, pVector, Size, RowNum);
        ParallelResultCalculation(pProcRows, pVector, pProcResult, Size, RowNum);
        ResultReplication(pProcResult, pResult, Size, RowNum);
        double Finish = MPI_Wtime();
        double ParallelTime = Finish - Start;

        if (ProcRank == 0) {
            // Обчислення прискорення
            double Speedup = SequentialTime / ParallelTime;
            std::cout << "Test: " << test + 1 << ": Size = " << Size
                      << ", Time: " << SequentialTime << " seconds"
                      << ", Processes: " << ProcNum
                      << ", Speedup: " << Speedup << std::endl;
        }

        delete[] pMatrix;
        delete[] pVector;
        delete[] pResult;
        delete[] pProcRows;
        delete[] pProcResult;
    }

    MPI_Finalize();
    return 0;
}
