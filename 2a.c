#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define nmax 32000
#define NUM_THREADS 16

unsigned long *create_sieve_to_number(unsigned long number);

int main() {
    unsigned long number;
    unsigned long *sieve;
    int *qtty;

    qtty = (int *) malloc(nmax * sizeof(int));
    sieve = create_sieve_to_number(nmax);

    // Usando um número fixo de threads
    //int NUM_THREADS= 2; // Ajuste conforme necessário
    omp_set_num_threads(NUM_THREADS);
    printf("Number of threads: %d\n", NUM_THREADS);


    // Sem escalonamento
    double start_time = omp_get_wtime();
    #pragma omp parallel for
    for (number = 2; number < nmax; number += 2) {
        for (unsigned long i = 2; i < number; i++) {
            if (sieve[i] == 1) {
                for (unsigned long j = i; j < number; j++) {
                    if (sieve[j] == 1) {
                        if (i + j == number) {
                            #pragma omp atomic
                            qtty[number]++;
                            break;
                        }
                    }
                }
                if (qtty[number] > 0) break;
            }
        }
    }
    double end_time = omp_get_wtime();
    printf("Without scheduling execution time: %.2f seconds\n", end_time - start_time);

    // Resetando qtty para o próximo teste
    for (int i = 0; i < nmax; i++) qtty[i] = 0;

    // STATIC schedule
    start_time = omp_get_wtime();
    #pragma omp parallel for schedule(static, 10000)
    for (number = 2; number < nmax; number += 2) {
        for (unsigned long i = 2; i < number; i++) {
            if (sieve[i] == 1) {
                for (unsigned long j = i; j < number; j++) {
                    if (sieve[j] == 1) {
                        if (i + j == number) {
                            #pragma omp atomic
                            qtty[number]++;
                            break;
                        }
                    }
                }
                if (qtty[number] > 0) break;
            }
        }
    }
    end_time = omp_get_wtime();
    printf("STATIC schedule execution time: %.2f seconds\n", end_time - start_time);

    // Resetando qtty para o próximo teste
    for (int i = 0; i < nmax; i++) qtty[i] = 0;

    // DYNAMIC schedule
    start_time = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic)
    for (number = 2; number < nmax; number += 2) {
        for (unsigned long i = 2; i < number; i++) {
            if (sieve[i] == 1) {
                for (unsigned long j = i; j < number; j++) {
                    if (sieve[j] == 1) {
                        if (i + j == number) {
                            #pragma omp atomic
                            qtty[number]++;
                            break;
                        }
                    }
                }
                if (qtty[number] > 0) break;
            }
        }
    }
    end_time = omp_get_wtime();
    printf("DYNAMIC schedule execution time: %.2f seconds\n", end_time - start_time);

    // Resetando qtty para o próximo teste
    for (int i = 0; i < nmax; i++) qtty[i] = 0;

    // GUIDED schedule
    start_time = omp_get_wtime();
    for (number = 2; number < nmax; number += 2) {
        for (unsigned long i = 2; i < number; i++) {
            if (sieve[i] == 1) {
                for (unsigned long j = i; j < number; j++) {
                    if (sieve[j] == 1) {
                        if (i + j == number) {
                            #pragma omp atomic
                            qtty[number]++;
                            break;
                        }
                    }
                }
                if (qtty[number] > 0) break;
            }
        }
    }
    end_time = omp_get_wtime();
    printf("GUIDED schedule execution time: %.2f seconds\n", end_time - start_time);

    // Limpeza
    free(qtty);
    free(sieve);

    return 0;
}

unsigned long *create_sieve_to_number(unsigned long number) {
    unsigned long *sieve = (unsigned long *)malloc(sizeof(unsigned long) * (number + 1));
    for (int i = 0; i < number; i++) {
        sieve[i] = 1;
    }
    for (unsigned long i = 2; i < number; i++) {
        if (sieve[i] == 1) {
            for (unsigned long j = i * i; j < number; j = j + i) {
                sieve[j] = 0;
            }
        }
    }
    return sieve;
}
