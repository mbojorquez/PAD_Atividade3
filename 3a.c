#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define N 2048 // Tamanho do tabuleiro
#define GENERATIONS 2000 // Número de gerações
#define SRAND_VALUE 1985
#define NUM_THREADS 16

// Função para contar vizinhos vivos
int getNeighbors(int **grid, int i, int j) {
    int count = 0;
    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            if (x == 0 && y == 0) continue; // Ignorar a própria célula
            int ni = (i + x + N) % N; // Contornar bordas
            int nj = (j + y + N) % N; // Contornar bordas
            count += grid[ni][nj];
        }
    }
    return count;
}

void initializeGrid(int **grid) {
    srand(SRAND_VALUE);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            grid[i][j] = rand() % 2; // Células aleatórias
        }
    }
}

int main() {
    // Alocação das matrizes
    int **grid = malloc(N * sizeof(int *));
    int **newgrid = malloc(N * sizeof(int *));
    for (int i = 0; i < N; i++) {
        grid[i] = malloc(N * sizeof(int));
        newgrid[i] = malloc(N * sizeof(int));
    }

    // Inicialização do tabuleiro
    initializeGrid(grid);

    // Contagem de células vivas na condição inicial
    int initial_count = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            initial_count += grid[i][j];
        }
    }
    printf("Initial live cells: %d\n", initial_count);

    // Versão paralela
    double start_time = omp_get_wtime();
    for (int gen = 0; gen < GENERATIONS; gen++) {
        omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                int alive_neighbors = getNeighbors(grid, i, j);
                if (grid[i][j] == 1) {
                    newgrid[i][j] = (alive_neighbors == 2 || alive_neighbors == 3) ? 1 : 0;
                } else {
                    newgrid[i][j] = (alive_neighbors == 3) ? 1 : 0;
                }
            }
        }
        // Troca de matrizes
        int **temp = grid;
        grid = newgrid;
        newgrid = temp;
    }
    double end_time = omp_get_wtime();
    printf("Parallel execution time: %.2f seconds\n", end_time - start_time);

    // Contagem de células vivas após todas as gerações
    int final_count_parallel = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            final_count_parallel += grid[i][j];
        }
    }
    printf("Final live cells after generations (parallel, %d threads): %d\n", NUM_THREADS, final_count_parallel);

    double elapsed_time = end_time - start_time;
    //269 seconds e o TEMPO PROMEDIO EN UM THREAD
    double speedup = 269 / elapsed_time;
    double efficiency = (speedup / NUM_THREADS) * 100;

    printf("Speedup: %.2f, Efficiency: %.0f%%\n", speedup, efficiency);

    // Liberação da memória
    for (int i = 0; i < N; i++) {
        free(grid[i]);
        free(newgrid[i]);
    }

    free(grid);
    free(newgrid);

    

    return 0;
}
