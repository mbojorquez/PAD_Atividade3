#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <papi.h>

#define N 1000 // Defina como 1000 o 10000 para testes maiores

// Função para multiplicar matrizes usando blocagem

void multiply_matrices(float **a, float **b, float **c, int block_size) {
    #pragma omp parallel for
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


// Função para liberar memória das matrizes
void free_matrices(float **matrix) {
    for (int i = 0; i < N; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

int main() {
    int EventSet = PAPI_NULL;
    long long values[6]; // Array to store event counts
    long long s,e;
    int retval;

    // Alocação de matrizes
    float **a = (float **)malloc(N * sizeof(float *));
    float **b = (float **)malloc(N * sizeof(float *));
    float **c = (float **)malloc(N * sizeof(float *));
    
    for (int i = 0; i < N; i++) {
        a[i] = (float *)malloc(N * sizeof(float));
        b[i] = (float *)malloc(N * sizeof(float));
        c[i] = (float *)calloc(N, sizeof(float)); // Inicializa com zeros
    }

    // Medir tempo para caso de thread única
    omp_set_num_threads(1);
    double start_time_single = omp_get_wtime();
    multiply_matrices(a, b, c, 32);
    double end_time_single = omp_get_wtime();
    double single_thread_time = end_time_single - start_time_single;

    // Print do tempo de execução de thread única
    printf("Single thread time: %.2f seconds\n", single_thread_time);


    // Initialize PAPI library
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        printf("PAPI_library_init failed: %d\n", retval);
        return 1;
    }

    // Create an event set
    retval = PAPI_create_eventset(&EventSet);
    if (retval != PAPI_OK) {
        printf("PAPI_create_eventset failed: %d\n", retval);
        return 1;
    }


    // Add cache miss events to the event set
    if (PAPI_add_event(EventSet, PAPI_L1_DCM) != PAPI_OK) {
        printf("PAPI_L1_DCM event addition failed\n");
        return 1;
    }

    // Testa a multiplicação com diferentes números de threads
    int num_threads[] = {1, 2, 4, 6, 8, 16};
    int num_tests = sizeof(num_threads) / sizeof(num_threads[0]);

    for (int t = 0; t < num_tests; t++) {
        omp_set_num_threads(num_threads[t]);
        double start_time = omp_get_wtime();

        /*// Inicia contadores PAPI
        PAPI_start_counters(events, 1);
        */
        
        // Start counting events
        retval = PAPI_start(EventSet);
        if (retval != PAPI_OK) {
            printf("PAPI_start failed");
            return 1;
        }

        s = PAPI_get_real_usec();
        
        // Multiplica as matrizes
        multiply_matrices(a, b, c, 64);

        double end_time = omp_get_wtime();
        /*PAPI_stop_counters(values, 1);*/
        e = PAPI_get_real_usec();

        // Stop counting events
        retval = PAPI_stop(EventSet, values);
        if (retval != PAPI_OK) {
            printf("Error in PAPI_stop");
            return 1;
        }

        long long cache_misses_L1 = values[0];

        double elapsed_time = end_time - start_time;

        // Calcular MFlops
        double mflops = (2.0 * N * N * N) / (elapsed_time * 1e6); // 2*N*N*N operações

        // Calcular Speedup e Eficiência
        double speedup = single_thread_time / elapsed_time;
        double efficiency = (speedup / num_threads[t]) * 100;

        // Print resultados
        printf("Threads: %d, Time: %.2f s, MFlops: %.2f, Speedup: %.2f, Efficiency: %.0f%%, Cache Misses: %lld\n", 
               num_threads[t], elapsed_time, mflops, speedup, efficiency, cache_misses_L1);

        // Reseta a matriz c para a próxima multiplicação
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                c[i][j] = 0.0; // Reinicializa a matriz resultado
            }
        }
    }

    // Libera a memória
    free_matrices(a);
    free_matrices(b);
    free_matrices(c);

    return 0;
}
