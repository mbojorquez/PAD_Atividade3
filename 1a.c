#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 2 // Definindo um tamanho pequeno para as matrizes

// Função para multiplicar matrizes usando blocagem
void multiply_matrices(float **a, float **b, float **c, int block_size) {
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < N; i += block_size) {
        for (int j = 0; j < N; j += block_size) {
            for (int k = 0; k < N; k += block_size) {
                for (int ii = i; ii < i + block_size && ii < N; ii++) {
                    for (int jj = j; jj < j + block_size && jj < N; jj++) {
                        float sum = 0.0;
                        for (int kk = k; kk < k + block_size && kk < N; kk++) {
                            sum += a[ii][kk] * b[kk][jj];
                        }
                        c[ii][jj] += sum;
                    }
                }
            }
        }
    }
}

// Função para imprimir uma matriz
void print_matrix(float **matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%.1f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main() {
    // Alocação de matrizes
    float **a = (float **)malloc(N * sizeof(float *));
    float **b = (float **)malloc(N * sizeof(float *));
    float **c = (float **)malloc(N * sizeof(float *));
    
    for (int i = 0; i < N; i++) {
        a[i] = (float *)malloc(N * sizeof(float));
        b[i] = (float *)malloc(N * sizeof(float));
        c[i] = (float *)calloc(N, sizeof(float)); // Inicializa com zeros
    }

    // Inicialização de matrizes com valores conhecidos
    a[0][0] = 1.0; a[0][1] = 2.0;
    a[1][0] = 3.0; a[1][1] = 4.0;
    
    b[0][0] = 5.0; b[0][1] = 6.0; 
    b[1][0] = 7.0; b[1][1] = 8.0; 

    printf("Matriz a:\n");
    print_matrix(a, N);
    printf("\nMatriz b:\n");
    print_matrix(b, N);

    // Teste de multiplicação com um número fixo de threads
    int num_threads[] = {1, 2, 4,6,8,16};
    int num_tests = sizeof(num_threads) / sizeof(num_threads[0]);

    for (int t = 0; t < num_tests; t++) {
        omp_set_num_threads(num_threads[t]);
        double start_time = omp_get_wtime();
        
        // Multiplicação das matrizes
        multiply_matrices(a, b, c, 2); // Tamanho do bloco de 2
        
        double end_time = omp_get_wtime();
        double elapsed_time = end_time - start_time;

        // Impressão dos resultados
        printf("\nResultado da multiplicação com %d thread(s):\n", num_threads[t]);
        print_matrix(c, N);
        printf("Tempo: %f segundos\n\n", elapsed_time);
        
        // Reseta a matriz c para a próxima multiplicação
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                c[i][j] = 0.0; // Reinicializa a matriz resultado
            }
        }
    }

    // Libera a memória
    for (int i = 0; i < N; i++) {
        free(a[i]);
        free(b[i]);
        free(c[i]);
    }
    free(a);
    free(b);
    free(c);

    return 0;
}
