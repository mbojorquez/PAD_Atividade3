#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

long long int fib(int n) {
    if (n < 2) return n;

    long long int x, y;

    // Usar tareas solo para n >= 30
    if (n >= 30) {
        #pragma omp parallel
        {
            #pragma omp single
            {
                #pragma omp task shared(x)
                {
                    x = fib(n - 1);
                }
                #pragma omp task shared(y)
                {
                    y = fib(n - 2);
                }
            }
        }
    } else {
        // Iterative approach for n < 30
        long long int n1 = 0, n2 = 1, n3;
        for (int i = 2; i <= n; ++i) {
            n3 = n1 + n2;
            n1 = n2;
            n2 = n3;
        }
        return n2; // Return the nth Fibonacci number directly
    }

    return x + y; // Return the result of the recursive calls
}


int main(int argc, char* argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Uso: %s <n>\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    long long int result;

    // Testa a multiplicação com diferentes números de threads
    int num_threads[] = {1,2, 4, 6, 8, 16};
    int num_tests = sizeof(num_threads) / sizeof(num_threads[0]);
    double single_thread_time=0.0;
    for (int t = 0; t < num_tests; t++) {
        omp_set_num_threads(num_threads[t]);
        double start_time = omp_get_wtime();
        result = fib(n);
        double end_time = omp_get_wtime();
        if(t==0){
            single_thread_time = end_time - start_time;
        }
        double elapsed_time = end_time - start_time;
        double speedup = single_thread_time / elapsed_time;
        double efficiency = (speedup / num_threads[t]) * 100;

        printf("Fibonacci de %d es %lld\n", n, result);
        printf("Tiempo de ejecución: %f segundos\n", end_time - start_time);

        printf("Threads: %d, Time: %.8f s, Speedup: %.2f, Efficiency: %.0f%%\n", num_threads[t], elapsed_time, speedup, efficiency);
    }
    return 0;
}